# QUERY SUBWORKFLOW (target: query_all).
#
# One timed run_query.py per (annotation, engine, query, param, rep). Engines
# are the naive `pyranges1` workflow and `gxf2parquet-<build>` for each entry
# of `query_builds` (reading that build's rep1 output).
#
# Output schema:
#   results/query/{annotation}/{engine}/{query}/{param}/rep{rep}/
#     time.txt     /usr/bin/time -v output
#     result.json  run_query.py's own numbers (op_wall_s, n_rows, ...)
#     cmdline.txt  the exact command line that was run
#     tool.log     stdout + stderr

import shlex


rule sample_genes:
    """Deterministic gene_name sample for the gene_set query (untimed)."""
    input:
        lambda w: build_dir(w.annotation, GENE_SAMPLE_BUILD, 1),
    output:
        f"{STAGE_DIR}/genes/{{annotation}}/n{{n}}.txt",
    params:
        parquet=lambda w: build_data(w.annotation, GENE_SAMPLE_BUILD, 1),
        seed=config.get("gene_set_seed", 42),
    resources:
        bench=1,
    shell:
        r"""
        python {SCRIPTS}/sample_genes.py --parquet {params.parquet:q} \
            --n {wildcards.n} --seed {params.seed} --out {output:q}
        """


def _query_data(wildcards) -> str:
    """The file/dataset the engine reads (passed to run_query + cache_prep)."""
    if wildcards.engine == NAIVE_ENGINE:
        return engine_input(wildcards.annotation, wildcards.engine)
    return build_data(wildcards.annotation, wildcards.engine.removeprefix("gxf2parquet-"), 1)


def _query_inputs(wildcards) -> dict:
    inputs = {"data": engine_input(wildcards.annotation, wildcards.engine)}
    if wildcards.query == "gene_set":
        inputs["genes"] = genes_path(wildcards.annotation, wildcards.param)
    return inputs


def _query_dir(wildcards) -> str:
    w = wildcards
    return f"{RESULTS_DIR}/query/{w.annotation}/{w.engine}/{w.query}/{w.param}/rep{w.rep}"


def _query_argv(wildcards) -> list[str]:
    w = wildcards
    argv = [
        "python", str(SCRIPTS / "run_query.py"),
        "--engine", "pyranges1" if w.engine == NAIVE_ENGINE else "gxf2parquet",
        "--input", _query_data(w),
        "--query", w.query,
        "--out-json", f"{_query_dir(w)}/result.json",
    ]
    if w.query == "chrom":
        argv += ["--chrom", w.param]
    elif w.query == "columns":
        argv += ["--columns", *config["column_sets"][w.param]]
    elif w.query == "gene_set":
        argv += ["--gene-names-file", genes_path(w.annotation, w.param)]
    return argv


rule run_query:
    input:
        unpack(_query_inputs),
    output:
        time=f"{RESULTS_DIR}/query/{{annotation}}/{{engine}}/{{query}}/{{param}}/rep{{rep}}/time.txt",
        result=f"{RESULTS_DIR}/query/{{annotation}}/{{engine}}/{{query}}/{{param}}/rep{{rep}}/result.json",
    log:
        cmdline=f"{RESULTS_DIR}/query/{{annotation}}/{{engine}}/{{query}}/{{param}}/rep{{rep}}/cmdline.txt",
        tool_log=f"{RESULTS_DIR}/query/{{annotation}}/{{engine}}/{{query}}/{{param}}/rep{{rep}}/tool.log",
    params:
        argv_str=lambda w: shlex.join(_query_argv(w)),
        data=_query_data,
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
        mkdir -p "$(dirname {output.time:q})"
        printf '%s\n' {params.argv_str:q} > {log.cmdline:q}
        if [[ -n "$src" ]]; then echo "# with PYTHONPATH=$src" >> {log.cmdline:q}; fi
        python {SCRIPTS}/cache_prep.py --mode {params.cache_mode} {params.data:q}
        command time -v -o {output.time:q} {params.argv_str} > {log.tool_log:q} 2>&1
        """
