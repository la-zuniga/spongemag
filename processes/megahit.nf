// processes/megahit.nf
process megahit {
    tag "megahit on $sample_id"
    publishDir "${params.outdir}/${sample_id}/megahit", mode: 'copy'

    container "/home/luiszuniga/work/containers/megahit_v1.2.9.sif"

    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    file("${sample_id}_contigs")

    script:
    """
    megahit -1 ${r1} -2 ${r2} -o ${sample_id}_contigs
    cd ${sample_id}_contigs
    cp final.contigs.fa contigs.fasta
    """
}
