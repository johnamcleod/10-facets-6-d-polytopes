"""Stage manifests: JSON records of pipeline stage runs."""

import json
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def git_hash():
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True
        )
        return result.stdout.strip()
    except Exception:
        return "unknown"


def python_version():
    return sys.version


def write_manifest(stage_dir, stage_name, params, counts, extra=None):
    """Write a manifest.json to stage_dir."""
    stage_dir = Path(stage_dir)
    stage_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "stage": stage_name,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_hash": git_hash(),
        "python_version": python_version(),
        "params": params,
        "counts": counts,
    }
    if extra:
        manifest.update(extra)
    path = stage_dir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return manifest


def read_manifest(stage_dir):
    path = Path(stage_dir) / "manifest.json"
    return json.loads(path.read_text())
