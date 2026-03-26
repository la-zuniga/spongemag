// processes/align.nf
process align {
    tag "SAMTOOLS & MINIMAP2: $sample_id"
    publishDir "${params.outdir}/${sample_id}/align", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)
    path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bai")

    script:
    """
    singularity exec ${params.containers.minimap2} \
        minimap2 -ax sr ${contigs}/contigs.fasta ${r1} ${r2} \
    | singularity exec ${params.containers.samtools} \
        samtools sort -o ${sample_id}.sorted.bam

    singularity exec ${params.containers.samtools} \
        samtools index ${sample_id}.sorted.bam ${sample_id}.sorted.bai
    """
}
