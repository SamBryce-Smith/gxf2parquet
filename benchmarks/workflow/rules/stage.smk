# Staging: a verified copy of each annotation (gzip) plus a decompressed GTF.
# Untimed; every rule here reserves bench=1.


rule stage_annotation:
    """Fetch/verify the annotation via download_annotation.py (the same script
    behind `pixi run download <name>`; a verified copy in data/ is re-used)
    and place it at stage/annotations/{annotation}/annotation.gtf.gz."""
    output:
        gz=f"{STAGE_DIR}/annotations/{{annotation}}/annotation.gtf.gz",
    wildcard_constraints:
        annotation=_alternation(FETCHED_ANNOTATIONS),
    resources:
        bench=1,
    shell:
        r"""
        set -euo pipefail
        src=$(python {SCRIPTS}/download_annotation.py \
            --config {ANNOTATIONS_YAML:q} \
            --name {wildcards.annotation} \
            --outdir {DATA_DIR:q} \
            --base-dir {WORKDIR:q} \
            --print-path)
        mkdir -p "$(dirname {output.gz:q})"
        # Copy (never link) so Snakemake's output cleanup can't touch the source.
        case "$src" in
            *.gz) cp "$src" {output.gz:q} ;;
            *)    gzip -c "$src" > {output.gz:q} ;;
        esac
        gzip -t {output.gz:q}
        """


rule scale_annotation:
    """Derived annotation (`scale_from` + `copies` in annotations.yaml): the
    source GTF repeated under renamed chromosomes and IDs (make_scaled_gtf.py)."""
    input:
        gz=lambda w: gtf_gz_path(ANNOTATION_SOURCES[w.annotation]["scale_from"]),
    output:
        gz=f"{STAGE_DIR}/annotations/{{annotation}}/annotation.gtf.gz",
    wildcard_constraints:
        annotation=_alternation(DERIVED_ANNOTATIONS),
    params:
        copies=lambda w: ANNOTATION_SOURCES[w.annotation]["copies"],
    resources:
        bench=1,
    shell:
        r"""
        python {SCRIPTS}/make_scaled_gtf.py --input {input.gz:q} \
            --copies {params.copies} --out {output.gz:q}
        gzip -t {output.gz:q}
        """


rule decompress_annotation:
    """Plain-text GTF: the naive engine's input (when naive_input: gtf) and
    the `gtf` row of disk.tsv."""
    input:
        gz=f"{STAGE_DIR}/annotations/{{annotation}}/annotation.gtf.gz",
    output:
        gtf=f"{STAGE_DIR}/annotations/{{annotation}}/annotation.gtf",
    resources:
        bench=1,
    shell:
        "gzip -dc {input.gz:q} > {output.gtf:q}"
