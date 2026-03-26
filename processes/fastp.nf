// processes/fastp.nf
process fastp {
    tag "FASTP: $sample_id"
    publishDir "${params.outdir}/${sample_id}/fastp", mode: 'copy'

    container "${params.containers.fastp}"

    input:
    tuple val(sample_id), file(x)

    output:
    tuple val(sample_id), file("${x[0].baseName}"), file("${x[1].baseName}")
    file("${sample_id}_report.html")

    script:

    """
    mkdir -p ${params.outdir}/${sample_id}/fastp
    fastp -i ${x[0]} -I ${x[1]} -o ${x[0].baseName} -O ${x[1].baseName} \
    -j ${sample_id}_report.json -h ${sample_id}_report.html 


    """
}
