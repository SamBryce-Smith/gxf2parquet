# BUILD SUBWORKFLOW (target: build_all).
#
# One timed `gxf2parquet build` per (annotation, build, rep), driven by the
# named parameter sets under `builds:` in the config, plus on-disk sizes of
# the GTF, its gzip and every build.
#
# Output schema:
#   results/build/{annotation}/{build}/rep{rep}/
#     time.txt     /usr/bin/time -v output (parsed by parse_gnu_time.py)
#     cmdline.txt  the exact command line that was run
#     tool.log     stdout + stderr from gxf2parquet
#   stage/parquet/{annotation}/{build}/rep{rep}/data.parquet
#     the build itself (file, or directory when partitioned); the query
#     subworkflow reads rep1

import shlex


def _build_argv(wildcards) -> list[str]:
    cfg = BUILDS[wildcards.build]
    preset = ANNOTATION_SOURCES[wildcards.annotation].get("preset", "gencode")
    argv = [
        "gxf2parquet", "build",
        gtf_gz_path(wildcards.annotation),
        build_data(wildcards.annotation, wildcards.build, wildcards.rep),
        "--preset", preset,
        "--compression", str(cfg.get("compression", "zstd")),
    ]
    if cfg.get("partition_cols"):
        argv += ["--partition-cols", *cfg["partition_cols"]]
    return argv


rule build_parquet:
    input:
        gz=lambda w: gtf_gz_path(w.annotation),
    output:
        time=f"{RESULTS_DIR}/build/{{annotation}}/{{build}}/rep{{rep}}/time.txt",
        parquet=directory(f"{STAGE_DIR}/parquet/{{annotation}}/{{build}}/rep{{rep}}"),
    log:
        # Declared as log: rather than output: so they survive Snakemake's
        # cleanup when the tool fails.
        cmdline=f"{RESULTS_DIR}/build/{{annotation}}/{{build}}/rep{{rep}}/cmdline.txt",
        tool_log=f"{RESULTS_DIR}/build/{{annotation}}/{{build}}/rep{{rep}}/tool.log",
    params:
        argv_str=lambda w: shlex.join(_build_argv(w)),
        cache_mode=CACHE_MODE,
        gxf2parquet_src=GXF2PARQUET_SRC,
    threads: 1
    resources:
        bench=100,
    shell:
        r"""
        set -euo pipefail
        # Optional: import gxf2parquet from another source tree (compare.sh).
        # Bound to a variable first: an empty {{...:q}} renders as nothing.
        src={params.gxf2parquet_src:q}
        if [[ -n "$src" ]]; then export PYTHONPATH="$src${{PYTHONPATH:+:$PYTHONPATH}}"; fi
        mkdir -p {output.parquet:q} "$(dirname {output.time:q})"
        printf '%s\n' {params.argv_str:q} > {log.cmdline:q}
        if [[ -n "$src" ]]; then echo "# with PYTHONPATH=$src" >> {log.cmdline:q}; fi
        python {SCRIPTS}/cache_prep.py --mode {params.cache_mode} {input.gz:q}
        # `command time` bypasses the bash keyword so we get GNU time (-v / -o).
        command time -v -o {output.time:q} {params.argv_str} > {log.tool_log:q} 2>&1
        """


rule disk_usage:
    """On-disk size of the GTF, its gzip, or a build (rep1)."""
    input:
        lambda w: artifact_path(w.annotation, w.artifact),
    output:
        f"{RESULTS_DIR}/disk/{{annotation}}/{{artifact}}.json",
    wildcard_constraints:
        artifact=_alternation(["gtf", "gtf_gz", *BUILDS]),
    resources:
        bench=1,
    shell:
        r"""
        python {SCRIPTS}/disk_usage.py --path {input:q} \
            --annotation {wildcards.annotation} --artifact {wildcards.artifact} \
            --out {output:q}
        """
