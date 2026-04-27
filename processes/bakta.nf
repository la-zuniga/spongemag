// processes/bakta.nf
//
// Bakta replaces Prokka for MAG / full-assembly annotation. Downstream
// KOfamScan + extract_ncycle_seqs.py consume `.faa` files only, and Bakta
// emits `${prefix}.faa` in the same layout Prokka did, so nothing
// downstream needs to change.
//
// Requires in nextflow.config:
//   params.containers.bakta    — Bakta Singularity image
//   params.bakta_db            — host path to the Bakta DB
//   params.bakta_threads       — threads for full-assembly bakta (default 16)
//   params.bakta_bin_threads   — threads per bin invocation (default 4)
//   params.bakta_bin_parallel  — how many bin invocations to run
//                                concurrently (default 4)
//
// Aggressive `--skip-*` flags are used because downstream only consumes CDS
// proteins. Disabling CRISPR / tRNA / ncRNA / plot detection is a pure
// speed win for our use case.


// ── Full-assembly Bakta (tier-agnostic, runs once per sample) ──────────────
process bakta_assembly {
    tag "BAKTA assembly: $sample_id"
    publishDir "${params.outdir}/${sample_id}/bakta/assembly", mode: 'copy'

    container "${params.containers.bakta}"
    containerOptions "--bind ${params.bakta_db}:${params.bakta_db}"

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("bakta_assembly_annotation")

    script:
    def threads = params.bakta_threads ?: 16
    """
    mkdir -p bakta_assembly_annotation

    bakta \\
        --db ${params.bakta_db} \\
        --output bakta_assembly_annotation \\
        --prefix ${sample_id} \\
        --threads ${threads} \\
        --skip-plot \\
        --skip-crispr \\
        --skip-trna \\
        --skip-tmrna \\
        --skip-ncrna \\
        --skip-ncrna-region \\
        --skip-ori \\
        --skip-gap \\
        --force \\
        ${contigs}/contigs.fasta
    """
}


// ── Per-tier bin Bakta ─────────────────────────────────────────────────────
// Bins are small, so we parallelize the per-bin loop via xargs -P.
// Each inner bakta uses fewer threads so several bins can run concurrently.
process bakta_bins {
    tag "BAKTA bins (${tier}): $sample_id"
    publishDir "${params.outdir}/${sample_id}/bakta/bins_${tier}", mode: 'copy'

    container "${params.containers.bakta}"
    containerOptions "--bind ${params.bakta_db}:${params.bakta_db}"

    input:
    tuple val(sample_id), val(tier), path(dastool_bins), path(checkm2_out)

    output:
    tuple val(sample_id), val(tier), path("bakta_bins_annotation")

    script:
    def bin_threads = params.bakta_bin_threads ?: 4
    def parallel_jobs = params.bakta_bin_parallel ?: 4
    """
    export TMPDIR=\$PWD

    awk '{ print \$1 }' ${checkm2_out}/${tier}_quality.tsv \\
        | tail -n +2 > bin_list.txt

    mkdir -p bakta_bins_annotation

    cat bin_list.txt | xargs -P ${parallel_jobs} -I {} \\
        bakta \\
            --db ${params.bakta_db} \\
            --output bakta_bins_annotation/{} \\
            --prefix {} \\
            --threads ${bin_threads} \\
            --skip-plot \\
            --skip-crispr \\
            --skip-trna \\
            --skip-tmrna \\
            --skip-ncrna \\
            --skip-ncrna-region \\
            --skip-ori \\
            --skip-gap \\
            --force \\
            ${dastool_bins}/{}.fa
    """
}
