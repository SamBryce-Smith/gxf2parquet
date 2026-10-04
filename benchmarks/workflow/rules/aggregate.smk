# Aggregation: per-run files -> TSVs -> replicate summaries -> plots.
# Every rule here reserves bench=1.

import json


rule host_info:
    output:
        f"{RESULTS_DIR}/host.json",
    resources:
        bench=1,
    shell:
        "python {SCRIPTS}/host_info.py > {output:q}"


rule dump_config:
    """Resolved config (after --config overrides), for merge_results.py and
    as a record of what was run. The config is a param, so editing it re-runs
    this rule."""
    output:
        f"{RESULTS_DIR}/run_config.json",
    params:
        cfg=json.dumps(dict(config), indent=2, sort_keys=True, default=str),
    resources:
        bench=1,
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        Path(output[0]).write_text(params.cfg + "\n")


rule merge_results:
    input:
        build_times=all_build_targets(),
        query_times=all_query_targets(),
        disk=all_disk_targets(),
        host=f"{RESULTS_DIR}/host.json",
        run_config=f"{RESULTS_DIR}/run_config.json",
    output:
        build=f"{RESULTS_DIR}/bench_build.tsv",
        query=f"{RESULTS_DIR}/bench_query.tsv",
        disk=f"{RESULTS_DIR}/disk.tsv",
    resources:
        bench=1,
    shell:
        r"""
        python {SCRIPTS}/merge_results.py \
            --host {input.host:q} \
            --run-config {input.run_config:q} \
            --build-times {input.build_times:q} \
            --query-times {input.query_times:q} \
            --disk {input.disk:q} \
            --out-build {output.build:q} \
            --out-query {output.query:q} \
            --out-disk {output.disk:q}
        """


rule summarize:
    """Collapse replicates: median + min/max per cell."""
    input:
        f"{RESULTS_DIR}/bench_{{kind}}.tsv",
    output:
        f"{RESULTS_DIR}/bench_{{kind}}_summary.tsv",
    wildcard_constraints:
        kind="build|query",
    resources:
        bench=1,
    shell:
        "python {SCRIPTS}/summarize.py --kind {wildcards.kind} --in {input:q} --out {output:q}"


rule plots:
    input:
        query=f"{RESULTS_DIR}/bench_query_summary.tsv",
        build=f"{RESULTS_DIR}/bench_build_summary.tsv",
        disk=f"{RESULTS_DIR}/disk.tsv",
    output:
        [f"{RESULTS_DIR}/plots/{p}" for p in PLOTS],
    resources:
        bench=1,
    shell:
        r"""
        mkdir -p "{RESULTS_DIR}/plots"
        Rscript {SCRIPTS}/plot.R {input.query:q} {input.build:q} {input.disk:q} "{RESULTS_DIR}/plots"
        """
