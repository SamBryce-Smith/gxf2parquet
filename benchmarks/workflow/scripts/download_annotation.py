#!/usr/bin/env python3
"""Fetch an annotation listed in annotations.yaml and validate its md5.

Used both by `pixi run download <name>` and by the pipeline's staging rule,
so there is a single code path for getting a verified annotation on disk.

For `url` entries the file is downloaded with curl into
<outdir>/<name>/<basename>, checked against the expected md5 (the pinned
`md5`, else the matching line of `md5sums_url`) with `md5sum`, checked with
`gzip -t`, and only then moved into place next to a `<basename>.md5`
sidecar. An existing file whose md5 still matches is not re-downloaded.

For `path` entries nothing is downloaded; the md5 is checked if pinned.

With --print-path the resolved local file path is printed to stdout (all
diagnostics go to stderr), which is how the staging rule consumes it.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


def log(msg: str) -> None:
    print(f"[download_annotation] {msg}", file=sys.stderr)


def md5sum(path: Path) -> str:
    out = subprocess.run(
        ["md5sum", str(path)], check=True, capture_output=True, text=True
    ).stdout
    return out.split()[0]


def curl(url: str, dest: Path | None = None) -> str:
    """Fetch url to dest (or return the body when dest is None)."""
    argv = [
        "curl",
        "-sSfL",
        "--retry",
        "10",
        "--retry-all-errors",
        "--retry-delay",
        "10",
    ]
    argv += [url] if dest is None else [url, "-o", str(dest)]
    return (
        subprocess.run(argv, check=True, capture_output=dest is None, text=True).stdout
        or ""
    )


def lookup_md5(md5sums_url: str, basename: str) -> str | None:
    """Find basename's checksum in an MD5SUMS-style file ("<md5>  <file>")."""
    for line in curl(md5sums_url).splitlines():
        parts = line.split()
        if len(parts) >= 2 and Path(parts[-1].lstrip("*")).name == basename:
            return parts[0].lower()
    return None


def sidecar_md5(path: Path) -> str | None:
    sidecar = path.with_name(path.name + ".md5")
    if sidecar.exists():
        text = sidecar.read_text().split()
        return text[0].lower() if text else None
    return None


def fetch_url_entry(name: str, entry: dict, outdir: Path, force: bool) -> Path:
    url = entry["url"]
    basename = url.rstrip("/").rsplit("/", 1)[-1]
    dest = outdir / name / basename
    pinned = (entry.get("md5") or "").strip().lower() or None

    # Re-use an existing download if it still matches a known-good checksum
    # (the pin, or the sidecar written after a previous verified download).
    known = pinned or sidecar_md5(dest)
    if dest.exists() and not force and known:
        if md5sum(dest) == known:
            log(f"{dest} present and md5 verified ({known}); skipping download")
            return dest
        log(f"{dest} exists but md5 does not match {known}; re-downloading")

    expected = pinned
    if expected is None:
        md5sums_url = entry.get("md5sums_url")
        if not md5sums_url:
            raise SystemExit(
                f"ERROR: annotation {name!r} has neither a pinned md5 nor an md5sums_url"
            )
        log(f"looking up checksum for {basename} in {md5sums_url}")
        expected = lookup_md5(md5sums_url, basename)
        if expected is None:
            raise SystemExit(f"ERROR: {basename} not listed in {md5sums_url}")

    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(basename + ".partial")
    partial.unlink(missing_ok=True)
    try:
        log(f"downloading {url}")
        curl(url, partial)
        actual = md5sum(partial)
        if actual != expected:
            raise SystemExit(
                f"ERROR: md5 mismatch for {basename}\n  expected: {expected}\n  actual:   {actual}"
            )
        if basename.endswith(".gz"):
            # Guards against e.g. an HTML error page served with a 200.
            subprocess.run(["gzip", "-t", str(partial)], check=True)
    except BaseException:
        partial.unlink(missing_ok=True)
        raise

    partial.replace(dest)
    dest.with_name(basename + ".md5").write_text(f"{expected}  {basename}\n")
    log(f"verified {dest} (md5 {expected})")
    if pinned is None:
        log(
            f"to pin this checksum, set `md5: {expected}` for {name!r} in annotations.yaml"
        )
    return dest


def check_path_entry(name: str, entry: dict, base_dir: Path) -> Path:
    path = Path(entry["path"])
    if not path.is_absolute():
        path = (base_dir / path).resolve()
    if not path.exists():
        raise SystemExit(f"ERROR: annotation {name!r}: {path} not found")
    pinned = (entry.get("md5") or "").strip().lower()
    if pinned:
        actual = md5sum(path)
        if actual != pinned:
            raise SystemExit(
                f"ERROR: md5 mismatch for {path}\n  expected: {pinned}\n  actual:   {actual}"
            )
        log(f"{path} md5 verified ({pinned})")
    return path


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--config", required=True, type=Path, help="annotations.yaml")
    p.add_argument("--name", required=True, help="annotation key in --config")
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("data"),
        help="download root (default: data)",
    )
    p.add_argument(
        "--base-dir",
        type=Path,
        default=Path.cwd(),
        help="directory that relative `path` entries are resolved "
        "against (default: current directory)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="re-download even if a verified copy exists",
    )
    p.add_argument(
        "--print-path",
        action="store_true",
        help="print the resolved local file path to stdout",
    )
    args = p.parse_args()

    with open(args.config) as fh:
        annotations = yaml.safe_load(fh) or {}
    if args.name not in annotations:
        raise SystemExit(
            f"ERROR: {args.name!r} not in {args.config} (known: {', '.join(annotations)})"
        )
    entry = annotations[args.name]

    if entry.get("scale_from"):
        raise SystemExit(
            f"ERROR: {args.name!r} is derived from {entry['scale_from']!r} "
            "(scale_from); the pipeline builds it, there is nothing to download"
        )
    try:
        if entry.get("url"):
            path = fetch_url_entry(args.name, entry, args.outdir, args.force)
        elif entry.get("path"):
            path = check_path_entry(args.name, entry, args.base_dir)
        else:
            raise SystemExit(
                f"ERROR: annotation {args.name!r} needs a `url` or a `path`"
            )
    except subprocess.CalledProcessError as e:
        raise SystemExit(
            f"ERROR: `{' '.join(map(str, e.cmd))}` failed with exit status {e.returncode}"
        ) from None

    if args.print_path:
        print(path.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
