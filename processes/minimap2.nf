// processes/minimap2.nf
process minimap2 {
    tag "minimap2 on $sample_id"
    publishDir "${params.outdir}/${sample_id}/minimap2", mode: 'copy'

    container "${params.containers.minimap2}"

    input:
    tuple val(sample_id), file(r1), file(r2)
    file(contigs)

    output:
    tuple val(sample_id), file("${sample_id}.sam")

    script:
    """
    minimap2 -ax sr ${contigs}/contigs.fasta ${r1} ${r2} > ${sample_id}.sam
    """
}
