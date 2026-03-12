// processes/checkm2.nf
process checkm2 {
    tag "CheckM2 on $sample_id"
    publishDir "${params.outdir}/${sample_id}/checkm2", mode: 'copy'

    input:
    val(sample_id)
    path(dastool_bins)

    output:
    val(sample_id)           // [0]
    file("checkm2_out")      // [1]

    script:
    """
    mkdir -p checkm2_out

   singularity exec ${params.containers.checkm2} \
    checkm2 predict \
    -x .fa \
    --threads ${params.checkm2_threads} \
    --input ${dastool_bins} \
    --output-directory checkm2_out \
    --database_path /db/uniref100.KO.1.dmnd

   sort_quality.py \
    -i checkm2_out/quality_report.tsv \
    -o checkm2_out/sorted_quality_report.tsv 
    """
}
