#!/usr/bin/env bash
# Benchmark two versions of gxf2parquet against each other and write a report.
#
#   ./compare.sh                                   # main vs this checkout, smoke config
#   ./compare.sh --before v0.1.0 --after my-branch
#   ./compare.sh --config config/performance.config.yaml --replicates 5
#
# Options (defaults in brackets):
#   --before REF      git ref for the baseline [main, else origin/main]
#   --after REF       git ref to compare; omit to use this checkout's src/ as is,
#                     including uncommitted changes [this checkout]
#   --config FILE     pipeline config [config/smoke.config.yaml]
#   --replicates N    replicates per timed cell [3]
#   --cores N         snakemake --cores [all]
#   --out DIR         output directory [results-compare]
#
# Both sides run the CURRENT pipeline in the same pixi environment; only the
# gxf2parquet source tree differs (via the gxf2parquet_src config key). The
# sides run one after the other, never concurrently. Result: <out>/report.md.
set -euo pipefail

cd "$(dirname "$0")"
BENCH_DIR="$PWD"
REPO="$(git rev-parse --show-toplevel)"

BEFORE=""
AFTER=""
CONFIG="config/smoke.config.yaml"
REPLICATES=3
CORES="all"
OUT="results-compare"
COMMAND="./compare.sh $*"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --before)     BEFORE="$2"; shift 2 ;;
    --after)      AFTER="$2"; shift 2 ;;
    --config)     CONFIG="$2"; shift 2 ;;
    --replicates) REPLICATES="$2"; shift 2 ;;
    --cores)      CORES="$2"; shift 2 ;;
    --out)        OUT="$2"; shift 2 ;;
    -h|--help)    sed -n '2,21p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

resolve_ref() {
  git -C "$REPO" rev-parse --verify --quiet "$1^{commit}" || {
    echo "unknown git ref: $1" >&2
    exit 2
  }
}

if [[ -z "$BEFORE" ]]; then
  if git -C "$REPO" rev-parse --verify --quiet "main^{commit}" > /dev/null; then
    BEFORE="main"
  else
    BEFORE="origin/main"
  fi
fi
resolve_ref "$BEFORE" > /dev/null
[[ -z "$AFTER" ]] || resolve_ref "$AFTER" > /dev/null
[[ -f "$CONFIG" ]] || { echo "config not found: $CONFIG" >&2; exit 1; }

mkdir -p "$OUT"
OUT_ABS="$(cd "$OUT" && pwd)"
TREES="$(mktemp -d "${TMPDIR:-/tmp}/gxf2parquet-compare.XXXXXX")"
cleanup() {
  for tree in "$TREES"/*/; do
    [[ -d "$tree" ]] && git -C "$REPO" worktree remove --force "$tree" > /dev/null 2>&1 || true
  done
  rm -rf "$TREES"
}
trap cleanup EXIT

checkout() {  # checkout <ref> <name> -> prints the tree's src/ path
  git -C "$REPO" worktree add --detach "$TREES/$2" "$1" > /dev/null 2>&1
  echo "$TREES/$2/src"
}

BEFORE_SRC="$(checkout "$BEFORE" before)"
if [[ -n "$AFTER" ]]; then
  AFTER_SRC="$(checkout "$AFTER" after)"
  AFTER_LABEL="$AFTER"
else
  AFTER_SRC="$REPO/src"
  AFTER_LABEL="this checkout"
fi

run_side() {  # run_side <side> <src>
  echo "==> [$1] gxf2parquet from $2"
  "$BENCH_DIR/run.sh" "$CONFIG" --cores "$CORES" -- \
    --config "gxf2parquet_src=$2" "replicates=$REPLICATES" \
    "results_dir=$OUT_ABS/$1" "stage_dir=$OUT_ABS/stage-$1"
}

run_side before "$BEFORE_SRC"
run_side after "$AFTER_SRC"

pixi run python "$BENCH_DIR/workflow/scripts/compare_report.py" \
  --before "$OUT_ABS/before" --after "$OUT_ABS/after" \
  --before-label "$BEFORE" --after-label "$AFTER_LABEL" \
  --command "cd benchmarks && $COMMAND" \
  --out "$OUT_ABS/report.md"
echo "==> report: $OUT_ABS/report.md"
