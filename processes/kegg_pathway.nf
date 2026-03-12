process kegg_pathway {
    tag "KEGG pathway mapping on $sample_id"
    publishDir "${params.outdir}/${sample_id}/kegg_pathways", mode: 'copy'

    container "${params.containers.kofamscan}"   // reuse — just needs python + requests

    input:
    val(sample_id)
    path(kofamscan_filtered)

    output:
    path("kegg_pathway_completeness.tsv")

    script:
    """
    map_kegg_api.py \
        -i ${kofamscan_filtered} \
        -o kegg_pathway_completeness.tsv \
        --min_tier confirmed \
        --min_completeness 0.1 \
        --target_pathways map00910 map00920 map00630
    """
}
