process quast {
    tag "QUAST (${tier}) on $sample_id"
    publishDir "${params.outdir}/${sample_id}/quast/${tier}", mode: 'copy'

    container "${params.containers.quast}"

    input:
    tuple val(sample_id), val(tier), path(dastool_bins), path(checkm2_out)

    output:
    tuple val(sample_id), val(tier), path("quast_output")

    script:
    """
    awk '{ print \$1 }' ${checkm2_out}/${tier}_quality.tsv \
        | tail -n +2 > bin_list.txt

    mkdir -p tier_bins
    while read bin; do
        [ -z "\$bin" ] && continue
        cp ${dastool_bins}/\${bin}.fa tier_bins/
    done < bin_list.txt

    if ls tier_bins/*.fa >/dev/null 2>&1; then
        quast -o quast_output tier_bins/*.fa
    else
        mkdir -p quast_output
        echo "No ${tier}-quality bins for ${sample_id}" \
            > quast_output/EMPTY.txt
    fi
    """
}
