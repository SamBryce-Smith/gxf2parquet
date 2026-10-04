#!/usr/bin/env bash
# Run the gxf2parquet benchmark pipeline.
#
#   ./run.sh                                     # default: performance config, target `all`
#   ./run.sh config/smoke.config.yaml            # smoke test on the in-repo test GTF
#   ./run.sh config/x.yaml -- build_all          # build subworkflow only
#   ./run.sh config/x.yaml -- query_all          # query subworkflow only
#   ./run.sh --dry-run                           # preview job graph
#   ./run.sh --cores 4                           # cap cores (default: all)
#
# Anything after `--` is forwarded verbatim to snakemake (targets, --config ...).
#
# Resource semantics: bench=100 is hard-coded. Every timed rule reserves the
# full pool (bench=100) so only ONE timed run happens at a time; every other
# rule reserves bench=1 so staging/aggregation can run during fill-in. Do not
# raise bench above 100: it would allow concurrent timed rules and contaminate
# the measurements.
set -euo pipefail

cd "$(dirname "$0")"

CONFIG_FILE=""
DRY_RUN=0
CORES="all"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run|-n) DRY_RUN=1; shift ;;
    --cores)      CORES="$2"; shift 2 ;;
    --cores=*)    CORES="${1#--cores=}"; shift ;;
    --) shift; EXTRA_ARGS+=("$@"); break ;;
    -h|--help)
      sed -n '2,13p' "$0"
      exit 0
      ;;
    -*)
      echo "unknown flag: $1" >&2
      exit 2
      ;;
    *)
      if [[ -z "$CONFIG_FILE" ]]; then
        CONFIG_FILE="$1"
      else
        echo "too many positional args: $1" >&2
        exit 2
      fi
      shift
      ;;
  esac
done

CONFIG_FILE="${CONFIG_FILE:-config/performance.config.yaml}"
[[ -f "$CONFIG_FILE" ]] || { echo "config not found: $CONFIG_FILE" >&2; exit 1; }

SNAKE_ARGS=(
  --configfile "$CONFIG_FILE"
  --cores "$CORES"
  --resources bench=100
  --rerun-incomplete
)
if [[ "$DRY_RUN" -eq 1 ]]; then
  SNAKE_ARGS+=(-n -p)
fi
# Bash 3.2 (macOS default) errors under `set -u` on empty-array expansion;
# the `${name[@]+...}` form expands only if name was set.
SNAKE_ARGS+=(${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"})

echo "==> pixi run snakemake ${SNAKE_ARGS[*]}"
exec pixi run snakemake "${SNAKE_ARGS[@]}"
