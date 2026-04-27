// processes/megahit.nf
process megahit {
    tag "MEGAHIT ASSEMBLYL: $sample_id"
    publishDir "${params.outdir}/${sample_id}/megahit", mode: 'copy'

    container "${params.containers.megahit}"
    // Only one MEGAHIT at a time — each run can peak at 40–80 GB on large inputs
    // and two concurrent assemblies would OOM a 124.8 GB host.
    maxForks 1

    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs")

    script:
    """
    megahit -1 ${r1} -2 ${r2} --memory 0.8 -o ${sample_id}_contigs
    cd ${sample_id}_contigs
    cp final.contigs.fa contigs.fasta
    """
}
