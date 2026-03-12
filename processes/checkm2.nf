// processes/checkm2.nf
process checkm2 {
    tag "CheckM2 on $sample_id"
    publishDir "${params.outdir}/${sample_id}/checkm2", mode: 'copy'

    input:
    val(sample_id)
    file(concoct_out)
    file(metabat_out)

    output:
    val(sample_id)           // [0]
    file("checkm2_out")      // [1]

    script:
    """
    mkdir -p checkm2_out/concoct_bins checkm2_out/metabat_bins

    # Run CheckM2 on CONCOCT bins
    singularity exec ${params.containers.checkm2} \
        checkm2 predict \
        -x .fa \
        --threads ${params.checkm2_threads} \
        --input ${concoct_out}/fasta_bins \
        --output-directory checkm2_out/concoct_bins \
        --database_path /db/uniref100.KO.1.dmnd

    # Run CheckM2 on MetaBAT2 bins
    singularity exec ${params.containers.checkm2} \
        checkm2 predict \
        -x .fa \
        --threads ${params.checkm2_threads} \
        --input ${metabat_out} \
        --output-directory checkm2_out/metabat_bins \
        --database_path /db/uniref100.KO.1.dmnd

    # Sort quality reports — only bins passing thresholds are kept
    sort_quality.py \
        -i checkm2_out/concoct_bins/quality_report.tsv \
        -o checkm2_out/concoct_bins/concoct_sorted_quality_report.tsv

    sort_quality.py \
        -i checkm2_out/metabat_bins/quality_report.tsv \
        -o checkm2_out/metabat_bins/metabat_sorted_quality_report.tsv
    """
}
