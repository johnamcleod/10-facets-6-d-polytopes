#!/usr/bin/env python3
"""Build the data archive (the Zenodo deposit) from a local tree.

Reads `MANIFEST` (one path or glob per line, relative to the repository root;
`#` starts a comment), checks that every entry matches at least one file, and
writes

    dist/coxeter-6-10-data.tar.gz     the files, with their paths
    dist/SHA256SUMS                   one line per file, for `shasum -a 256 -c`

Run:  python3 pack_data.py [--list]
"""
import hashlib
import sys
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIST = ROOT / "dist"


def files():
    out, missing = [], []
    for line in (ROOT / "MANIFEST").read_text().splitlines():
        pat = line.split("#", 1)[0].strip()
        if not pat:
            continue
        hits = sorted(p for p in ROOT.glob(pat) if p.is_file())
        if not hits:
            missing.append(pat)
        out += hits
    if missing:
        raise SystemExit(f"MANIFEST entries with no file: {missing}")
    return sorted(set(out))


def main():
    fs = files()
    if "--list" in sys.argv:
        for f in fs:
            print(f.relative_to(ROOT))
        return 0
    DIST.mkdir(exist_ok=True)
    sums = []
    with tarfile.open(DIST / "coxeter-6-10-data.tar.gz", "w:gz") as tar:
        for f in fs:
            rel = f.relative_to(ROOT)
            tar.add(f, arcname=str(rel))
            sums.append(f"{hashlib.sha256(f.read_bytes()).hexdigest()}  {rel}")
    (DIST / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    size = (DIST / "coxeter-6-10-data.tar.gz").stat().st_size
    print(f"{len(fs)} files -> dist/coxeter-6-10-data.tar.gz ({size / 2**20:.0f} MB), "
          f"dist/SHA256SUMS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
