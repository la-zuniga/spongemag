process multiqc {
    tag "MultiQC on $sample_id"
    publishDir "${params.outdir}/${sample_id}/multiqc", mode: 'copy'

    container "${params.containers.multiqc}"

    input:
    tuple val(sample_id), path(qc_dir)

    output:
    path("MultiQC_report.html")

    script:
    """
    multiqc ${qc_dir} -o . -n MultiQC_report.html
    """
}
