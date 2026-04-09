import json
import tarfile
import tempfile
import unittest
from pathlib import Path

from virusdetect.db import (
    LEGACY_MARKER_FILES,
    LEGACY_REFERENCE_FILES,
    build_database_manifest,
    bundle_database_archive,
    verify_database_files,
)


class DatabaseBundleTests(unittest.TestCase):
    def populate_legacy_database(self, root: Path) -> None:
        root.mkdir(parents=True, exist_ok=True)
        for file_name in LEGACY_MARKER_FILES + LEGACY_REFERENCE_FILES:
            (root / file_name).write_text("x\n", encoding="utf-8")

    def test_build_database_manifest_records_legacy_layout(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            source = Path(tmp_dir) / "databases"
            self.populate_legacy_database(source)

            manifest = build_database_manifest(str(source), db_version="2026.04")

        self.assertEqual(manifest["db_version"], "2026.04")
        self.assertEqual(manifest["source_kind"], "legacy")
        self.assertIn("vrl_plant", manifest["files"])
        self.assertIn("vrl_plant_prot", manifest["files"])
        self.assertNotIn("manifest.json", manifest["files"])

    def test_bundle_database_archive_creates_installable_v2_bundle(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp = Path(tmp_dir)
            source = tmp / "databases"
            self.populate_legacy_database(source)

            bundle = bundle_database_archive(
                source_dir=str(source),
                destination_archive=str(tmp / "dist" / "virusdetect-db-test.tar.gz"),
                db_version="2026.04",
                db_name="VirusDetect plant references",
            )

            archive_path = Path(bundle.archive_path)
            sha256_path = Path(bundle.sha256_path)
            self.assertTrue(archive_path.exists())
            self.assertTrue(sha256_path.exists())
            self.assertIn(archive_path.name, sha256_path.read_text(encoding="utf-8"))

            extract_root = tmp / "extract"
            extract_root.mkdir()
            with tarfile.open(archive_path, "r:gz") as archive_handle:
                archive_handle.extractall(extract_root)

            installed_db = extract_root / "database"
            manifest = json.loads((installed_db / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(verify_database_files(str(installed_db)), {})
            self.assertEqual(manifest["db_name"], "VirusDetect plant references")
            self.assertEqual(manifest["db_version"], "2026.04")
            self.assertIn("vrl_genbank.info.gz", manifest["files"])


if __name__ == "__main__":
    unittest.main()
