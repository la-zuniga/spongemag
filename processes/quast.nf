process quast {
    tag "QUAST on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/quast", mode: 'copy'

    container "${params.containers.quast}"

    input:
    val(sample_id)
    path(dastool_bins)
    
    output:
    file("quast_output")

    script:
    """
    # Running QUAST on the MetaBAT bins

    
    quast -o quast_output ${dastool_bins}/*.fa
    """
}
