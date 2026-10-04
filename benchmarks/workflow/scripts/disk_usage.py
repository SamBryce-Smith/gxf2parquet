#!/usr/bin/env python3
"""Record the on-disk size of an artifact (a file or a partitioned directory)."""

import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--path", required=True, type=Path)
    p.add_argument("--annotation", required=True)
    p.add_argument("--artifact", required=True, help="gtf, gtf_gz or a build name")
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    path = args.path.resolve()  # follow a staged symlink to the real file
    if path.is_dir():
        files = [f for f in path.rglob("*") if f.is_file()]
    elif path.is_file():
        files = [path]
    else:
        raise SystemExit(f"ERROR: not found: {args.path}")

    record = {
        "annotation": args.annotation,
        "artifact": args.artifact,
        "bytes": sum(f.stat().st_size for f in files),
        "n_files": len(files),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
