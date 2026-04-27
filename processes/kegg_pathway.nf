process kegg_pathway {
    tag "KEGG pathway mapping (${label}) on $sample_id"
    publishDir "${params.outdir}/${sample_id}/kegg_pathways/${label}", mode: 'copy'

    container "${params.containers.python}"   // reuse — just needs python + requests

    input:
    tuple val(sample_id), val(label), path(kofamscan_filtered)

    output:
    tuple val(sample_id), val(label),
          path("${label}_kegg_pathway_completeness.tsv"),
          path("${label}_nitrogen_cycle_completeness.tsv")

    script:
    """
    map_kegg_api.py \
        -i ${kofamscan_filtered} \
        -o ${label}_kegg_pathway_completeness.tsv \
        --min_tier confirmed \
        --min_completeness 0.05 \
        --target_pathways map00910 map00920 map00630 \
        --nitrogen_output ${label}_nitrogen_cycle_completeness.tsv
    """
}
