#!/usr/bin/env python3
"""Put the page cache into a known state before a timed run (untimed).

Linux keeps recently read file data in RAM (the page cache). Without a fixed
policy, whichever engine/replicate happens to run first pays for reading its
input from disk and later runs get it for free, so results depend on run
order. Modes:

  warm   read every input file once, so each timed run starts from RAM
         (measures parse/filter cost; portable, no privileges)
  evict  posix_fadvise(POSIX_FADV_DONTNEED) on every input file, so the timed
         run reads them cold from disk (Linux only, no root needed; only the
         files under test are evicted)
  drop   `sync; echo 3 > /proc/sys/vm/drop_caches` via `sudo -n` (riker's
         approach: whole-system cold cache; needs passwordless sudo)

Directories (partitioned Parquet datasets) are expanded recursively.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

CHUNK = 8 * 1024 * 1024


def iter_files(paths: list[Path]):
    for p in paths:
        if p.is_dir():
            yield from sorted(f for f in p.rglob("*") if f.is_file())
        elif p.is_file():
            yield p
        else:
            raise SystemExit(f"ERROR: input not found: {p}")


def warm(files) -> None:
    for f in files:
        with open(f, "rb") as fh:
            while fh.read(CHUNK):
                pass


def evict(files) -> None:
    if not hasattr(os, "posix_fadvise"):
        raise SystemExit("ERROR: cache_mode=evict needs posix_fadvise (Linux)")
    os.sync()  # dirty pages can't be dropped; flush them first
    for f in files:
        fd = os.open(f, os.O_RDONLY)
        try:
            os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def drop() -> None:
    subprocess.run(["sync"], check=True)
    # -n: fail immediately instead of hanging on a password prompt mid-run.
    result = subprocess.run(
        ["sudo", "-n", "tee", "/proc/sys/vm/drop_caches"],
        input="3\n",
        text=True,
        stdout=subprocess.DEVNULL,
    )
    if result.returncode != 0:
        raise SystemExit(
            "ERROR: cache_mode=drop needs passwordless sudo for "
            "`tee /proc/sys/vm/drop_caches` (see benchmarks/README.md)"
        )


def main() -> int:
    p = argparse.ArgumentParser(
        description="Prepare the page cache before a timed run."
    )
    p.add_argument("--mode", required=True, choices=["warm", "evict", "drop"])
    p.add_argument("paths", nargs="+", type=Path, help="input files or directories")
    args = p.parse_args()

    files = list(iter_files(args.paths))
    if args.mode == "warm":
        warm(files)
    elif args.mode == "evict":
        evict(files)
    else:
        drop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
