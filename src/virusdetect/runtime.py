from __future__ import annotations

import os
import subprocess
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ResolvedTool:
    name: str
    required: bool
    path: str | None
    source: str | None

    @property
    def found(self) -> bool:
        return self.path is not None


TOOL_SPECS = (
    ("bwa", True),
    ("samtools", True),
    ("blastn", True),
    ("blastx", True),
    ("makeblastdb", True),
    ("hisat2", False),
    ("hisat2-build", False),
    ("spades.py", False),
)


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def is_runnable_binary(path: str) -> bool:
    if not os.access(path, os.X_OK):
        return False

    try:
        completed = subprocess.run(
            [path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return False

    return completed.returncode >= 0


def resolve_tool(name: str) -> tuple[str | None, str | None]:
    system_path = shutil.which(name)
    if system_path and is_runnable_binary(system_path):
        return system_path, "PATH"

    return None, None


def resolve_tools():
    resolved = []
    for name, required in TOOL_SPECS:
        path, source = resolve_tool(name)
        resolved.append(ResolvedTool(name=name, required=required, path=path, source=source))
    return resolved


def tool_path_map() -> dict[str, str]:
    return {status.name: status.path for status in resolve_tools() if status.path is not None}


def missing_required_tools(statuses: list[ResolvedTool]) -> list[str]:
    return [status.name for status in statuses if status.required and not status.found]
