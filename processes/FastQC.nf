process fastqc {
    tag "FastQC on $sample_id"
    publishDir "${params.outdir}/${sample_id}/fastqc", mode: 'copy'

    container "/home/luiszuniga/work/containers/FastQC_v0.12.1.sif"

    input:
    tuple val(sample_id), file(r1), file(r2)

    output:
    tuple val(sample_id), path("fastqc_results")

    script:
    """
    mkdir fastqc_results
    fastqc ${r1} ${r2} --outdir fastqc_results
    """
}
