// processes/samtools.nf
process samtools {
    tag "samtools on $sample_id"
    publishDir "${params.outdir}/${sample_id}/samtools", mode: 'copy'

    container "/home/luiszuniga/work/containers/samtools_v1.21.sif"

    input:
    tuple val(sample_id), file(alignment)

    output:
    tuple val(sample_id), file("${sample_id}.bam"), file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bai")

    script:
    """
    samtools view -S -b ${alignment} > ${sample_id}.bam 
    samtools sort ${sample_id}.bam -o ${sample_id}.sorted.bam 
    samtools index ${sample_id}.sorted.bam ${sample_id}.sorted.bai
    """
}
