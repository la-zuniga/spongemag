process KOfamscan {
    tag "KOfamscan on $sample_id metagenomes (translated)"
    publishDir "${params.outdir}/${sample_id}/kofamscan", mode: 'copy'

    container "${params.containers.kofamscan}"

    input:
    val(sample_id)
    path(prokka_bins_annotation) 

    output:
    path("kofamscan_bins_annotation")
    path("kofamscan_filtered.tsv")

    script:

    """
    mkdir kofamscan_bins_annotation
    mkdir -p tmp/tabular
    mkdir -p tmp/mapper
    echo "Current working directory: \$(pwd)"


    find ${prokka_bins_annotation}/ -type f -name "*.faa" \
        | parallel -j ${params.kofamscan_threads} \
            'exec_annotation \
                -o kofamscan_bins_annotation/{/.}.out \
                -p /kofamscan/db/profiles/ncycle.hal \
                --cpu 1 {}'
  filter_kofamscan.py \
    -i kofamscan_bins_annotation \
    -o kofamscan_filtered.tsv \
    --min_tier putative \
    --evalue_cutoff 1e-5 \
    --min_putative_score 50
    """
}
