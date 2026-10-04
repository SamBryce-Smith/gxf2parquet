# Staging: a verified copy of each annotation (gzip) plus a decompressed GTF.
# Untimed; every rule here reserves bench=1.


rule stage_annotation:
    """Fetch/verify the annotation via download_annotation.py (the same script
    behind `pixi run download <name>`; a verified copy in data/ is re-used)
    and place it at stage/annotations/{annotation}/annotation.gtf.gz."""
    output:
        gz=f"{STAGE_DIR}/annotations/{{annotation}}/annotation.gtf.gz",
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
