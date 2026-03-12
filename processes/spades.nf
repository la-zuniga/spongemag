// processes/spades.nf
process spades {
    tag "spades on $sample_id"
    publishDir "${params.outdir}/${sample_id}/spades", mode: 'copy'

    container "${params.containers.spades}"

    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    file("${sample_id}_spades_assembly")

    script:
    """
    spades.py -1 ${r1} -2 ${r2} -o ${sample_id}_spades_assembly --meta
    """
}
