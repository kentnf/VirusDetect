from __future__ import annotations

import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

from virusdetect import __db_release_repo__, __db_version__, __version__

DEFAULT_DB_ENV_VAR = "VIRUSDETECT_DB_DIR"
DEFAULT_DB_RELEASE_REPO_ENV_VAR = "VIRUSDETECT_DB_RELEASE_REPO"
DEFAULT_DB_RELEASE_TAG_ENV_VAR = "VIRUSDETECT_DB_RELEASE_TAG"
DEFAULT_DB_ASSET_NAME_ENV_VAR = "VIRUSDETECT_DB_ASSET_NAME"
DEFAULT_DB_SHA256_ASSET_NAME_ENV_VAR = "VIRUSDETECT_DB_SHA256_ASSET_NAME"
USER_DB_DIR = Path.home() / ".local" / "share" / "virusdetect" / "database"
MANIFEST_NAME = "manifest.json"
LEGACY_MARKER_FILES = (
    "vrl_genbank.info.gz",
    "vrl_idmapping.gz",
)
LEGACY_REFERENCE_FILES = (
    "vrl_plant",
    "vrl_plant_prot",
)


@dataclass(frozen=True)
class DatabaseLocation:
    path: str
    source: str
    kind: str


@dataclass(frozen=True)
class DatabaseBundle:
    archive_path: str
    sha256_path: str


@dataclass(frozen=True)
class DatabaseDownloadSource:
    url: str
    sha256: str | None
    sha256_url: str | None
    release_repo: str | None
    release_tag: str | None
    asset_name: str | None
    sha256_asset_name: str | None


def read_manifest(db_path: str) -> dict | None:
    manifest_path = Path(db_path) / MANIFEST_NAME
    if not manifest_path.exists():
        return None

    with manifest_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, dict):
        raise ValueError(f"Invalid database manifest: {manifest_path}")
    return data


def is_v2_database_dir(db_path: str) -> bool:
    try:
        manifest = read_manifest(db_path)
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    return manifest is not None


def is_legacy_database_dir(db_path: str) -> bool:
    candidate = Path(db_path)
    if not candidate.is_dir():
        return False

    required = LEGACY_MARKER_FILES + LEGACY_REFERENCE_FILES
    return all((candidate / file_name).exists() for file_name in required)


def detect_database_kind(db_path: str) -> str | None:
    if is_v2_database_dir(db_path):
        return "v2"
    if is_legacy_database_dir(db_path):
        return "legacy"
    return None


def iter_database_candidates():
    env_db_path = os.getenv(DEFAULT_DB_ENV_VAR)
    if env_db_path:
        yield DatabaseLocation(path=env_db_path, source=DEFAULT_DB_ENV_VAR, kind="unknown")

    cwd_database = Path.cwd() / "database"
    yield DatabaseLocation(path=str(cwd_database), source="cwd:database", kind="unknown")

    cwd_databases = Path.cwd() / "databases"
    yield DatabaseLocation(path=str(cwd_databases), source="cwd:databases", kind="unknown")

    conda_prefix = os.getenv("CONDA_PREFIX")
    if conda_prefix:
        yield DatabaseLocation(
            path=str(Path(conda_prefix) / "share" / "virusdetect" / "database"),
            source="conda",
            kind="unknown",
        )

    yield DatabaseLocation(path=str(USER_DB_DIR), source="user", kind="unknown")


def resolve_database_location() -> DatabaseLocation | None:
    for candidate in iter_database_candidates():
        db_kind = detect_database_kind(candidate.path)
        if db_kind is not None:
            return DatabaseLocation(path=candidate.path, source=candidate.source, kind=db_kind)
    return None


def resolve_database_target_dir(custom_path: str | None = None) -> str:
    if custom_path:
        return str(Path(custom_path))

    conda_prefix = os.getenv("CONDA_PREFIX")
    if conda_prefix:
        return str(Path(conda_prefix) / "share" / "virusdetect" / "database")

    return str(USER_DB_DIR)


def resolve_reference_path(reference: str, db_path: str) -> str:
    candidate = Path(reference)
    if candidate.exists():
        return str(candidate.resolve())

    db_candidate = Path(db_path) / reference
    if db_candidate.exists():
        return str(db_candidate.resolve())

    raise FileNotFoundError(f"Could not resolve reference '{reference}' in {db_path}")


def resolve_data_path(name_or_path: str, db_path: str) -> str:
    candidate = Path(name_or_path)
    if candidate.exists():
        return str(candidate.resolve())

    db_candidate = Path(db_path) / name_or_path
    if db_candidate.exists():
        return str(db_candidate.resolve())

    raise FileNotFoundError(f"Could not resolve data path '{name_or_path}' in {db_path}")


def resolve_seq_info_path(db_path: str) -> str:
    database_dir = Path(db_path)
    candidates = (
        database_dir / "vrl_genbank_info.gz",
        database_dir / "vrl_genbank.info.gz",
    )
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())
    raise FileNotFoundError(f"Could not find virus sequence info in {db_path}")


def verify_database_files(db_path: str) -> dict[str, list[str]]:
    candidate = Path(db_path)
    if not candidate.exists():
        return {"database": [f"missing directory: {candidate}"]}

    manifest = read_manifest(str(candidate))
    if manifest is not None:
        files = manifest.get("files", [])
        if not isinstance(files, list):
            return {"database": [f"invalid manifest files entry: {candidate / MANIFEST_NAME}"]}
        missing = [file_name for file_name in files if not (candidate / file_name).exists()]
        return {"database": missing} if missing else {}

    missing = [
        file_name
        for file_name in LEGACY_MARKER_FILES + LEGACY_REFERENCE_FILES
        if not (candidate / file_name).exists()
    ]
    return {"database": missing} if missing else {}


def collect_database_files(db_path: str) -> list[str]:
    candidate = Path(db_path)
    return [
        path.relative_to(candidate).as_posix()
        for path in sorted(candidate.rglob("*"))
        if path.is_file() and path.name != MANIFEST_NAME
    ]


def build_database_manifest(db_path: str, db_version: str, db_name: str = "virusdetect") -> dict:
    candidate = Path(db_path)
    if not candidate.is_dir():
        raise FileNotFoundError(f"Database source directory was not found: {candidate}")

    source_kind = detect_database_kind(str(candidate))
    if source_kind is None:
        raise ValueError(f"Database source directory is not a recognized VirusDetect layout: {candidate}")

    return {
        "format": "virusdetect-database-manifest",
        "manifest_version": 1,
        "db_name": db_name,
        "db_version": db_version,
        "source_kind": source_kind,
        "created_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "files": collect_database_files(str(candidate)),
    }


def default_database_asset_name(db_version: str) -> str:
    return f"virusdetect-db-{db_version}.tar.gz"


def resolve_database_download_source(
    url: str | None = None,
    sha256: str | None = None,
    sha256_url: str | None = None,
    release_repo: str | None = None,
    release_tag: str | None = None,
    asset_name: str | None = None,
    sha256_asset_name: str | None = None,
    db_version: str = __db_version__,
) -> DatabaseDownloadSource:
    if url:
        return DatabaseDownloadSource(
            url=url,
            sha256=sha256,
            sha256_url=sha256_url,
            release_repo=release_repo,
            release_tag=release_tag,
            asset_name=asset_name,
            sha256_asset_name=sha256_asset_name,
        )

    resolved_release_repo = release_repo or os.getenv(DEFAULT_DB_RELEASE_REPO_ENV_VAR) or __db_release_repo__
    resolved_release_tag = release_tag or os.getenv(DEFAULT_DB_RELEASE_TAG_ENV_VAR) or f"v{__version__}"
    resolved_asset_name = asset_name or os.getenv(DEFAULT_DB_ASSET_NAME_ENV_VAR) or default_database_asset_name(db_version)
    resolved_sha256_asset_name = (
        sha256_asset_name
        or os.getenv(DEFAULT_DB_SHA256_ASSET_NAME_ENV_VAR)
        or f"{resolved_asset_name}.sha256"
    )

    base_url = f"https://github.com/{resolved_release_repo}/releases/download/{resolved_release_tag}"
    return DatabaseDownloadSource(
        url=f"{base_url}/{resolved_asset_name}",
        sha256=sha256,
        sha256_url=sha256_url or f"{base_url}/{resolved_sha256_asset_name}",
        release_repo=resolved_release_repo,
        release_tag=resolved_release_tag,
        asset_name=resolved_asset_name,
        sha256_asset_name=resolved_sha256_asset_name,
    )


def download_file(url: str, destination_path: Path) -> None:
    with urllib.request.urlopen(url) as response, destination_path.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def compute_sha256(file_path: Path) -> str:
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_expected_sha256(sha256_file: Path) -> str:
    with sha256_file.open("r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    return first_line.split()[0]


def write_sha256_file(file_path: Path, sha256: str) -> Path:
    sha256_path = Path(f"{file_path}.sha256")
    sha256_path.write_text(f"{sha256}  {file_path.name}\n", encoding="utf-8")
    return sha256_path


def find_database_source_dir(extract_root: Path) -> Path:
    candidates = [extract_root]
    candidates.extend(path for path in extract_root.rglob("database") if path.is_dir())
    candidates.extend(path for path in extract_root.rglob("databases") if path.is_dir())

    for candidate in candidates:
        if detect_database_kind(str(candidate)) is not None:
            return candidate

    raise FileNotFoundError("Could not find a valid VirusDetect database directory in the extracted archive")


def extract_archive(archive_path: Path, destination_dir: Path) -> None:
    destination_root = destination_dir.resolve()

    with tarfile.open(archive_path, "r:gz") as archive_handle:
        for member in archive_handle.getmembers():
            if member.issym() or member.islnk():
                raise ValueError(f"Archive member is not supported: {member.name}")

            member_path = (destination_root / member.name).resolve()
            if os.path.commonpath([str(destination_root), str(member_path)]) != str(destination_root):
                raise ValueError(f"Archive member escapes destination directory: {member.name}")

        archive_handle.extractall(destination_root)


def install_database_archive(url: str, destination_dir: str, sha256: str | None = None, sha256_url: str | None = None) -> str:
    if not url:
        raise ValueError("A database archive URL is required")

    target_dir = Path(destination_dir).resolve()
    target_dir.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="virusdetect-db-") as tmp_dir:
        tmp_root = Path(tmp_dir)
        archive_path = tmp_root / "database.tar.gz"
        download_file(url, archive_path)

        expected_sha256 = sha256
        if expected_sha256 is None and sha256_url:
            sha256_path = tmp_root / "database.sha256"
            download_file(sha256_url, sha256_path)
            expected_sha256 = read_expected_sha256(sha256_path)

        if expected_sha256:
            actual_sha256 = compute_sha256(archive_path)
            if actual_sha256 != expected_sha256:
                raise ValueError(
                    f"SHA256 mismatch for {url}: expected {expected_sha256}, got {actual_sha256}"
                )

        extract_root = tmp_root / "extract"
        extract_root.mkdir()
        extract_archive(archive_path, extract_root)
        source_dir = find_database_source_dir(extract_root)

        if target_dir.exists():
            shutil.rmtree(target_dir)

        shutil.copytree(source_dir, target_dir)

    return str(target_dir)


def bundle_database_archive(
    source_dir: str,
    destination_archive: str,
    db_version: str,
    db_name: str = "virusdetect",
) -> DatabaseBundle:
    source_path = Path(source_dir).resolve()
    if detect_database_kind(str(source_path)) is None:
        raise ValueError(f"Database source directory is not a recognized VirusDetect layout: {source_path}")

    archive_path = Path(destination_archive).resolve()
    archive_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="virusdetect-db-bundle-") as tmp_dir:
        tmp_root = Path(tmp_dir)
        bundle_root = tmp_root / "database"
        shutil.copytree(source_path, bundle_root)

        manifest = build_database_manifest(str(bundle_root), db_version=db_version, db_name=db_name)
        (bundle_root / MANIFEST_NAME).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

        with tarfile.open(archive_path, "w:gz") as archive_handle:
            archive_handle.add(bundle_root, arcname="database")

    sha256 = compute_sha256(archive_path)
    sha256_path = write_sha256_file(archive_path, sha256)
    return DatabaseBundle(archive_path=str(archive_path), sha256_path=str(sha256_path))
