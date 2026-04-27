// processes/kofamscan.nf
//
// Generic KOfamscan process. `label` drives publishDir + output filenames so
// the same module can run per bin tier (bins_hq, bins_mq) and on the full
// assembly.

process KOfamscan {
    tag "KOfamscan (${label}) on $sample_id"
    publishDir "${params.outdir}/${sample_id}/kofamscan/${label}", mode: 'copy'

    container "${params.containers.kofamscan}"

    input:
    tuple val(sample_id), val(label), path(prokka_annotation)

    output:
    tuple val(sample_id), val(label),
          path("kofamscan_${label}_annotation"),
          path("${sample_id}_${label}_kofamscan_filtered.tsv")

    script:
    """
    export TMPDIR=\$PWD
    mkdir kofamscan_${label}_annotation
    mkdir -p tmp/tabular
    mkdir -p tmp/mapper
    echo "Current working directory: \$(pwd)"


    find ${prokka_annotation}/ -type f -name "*.faa" \
        | parallel -j ${params.kofamscan_threads} \
            'mkdir -p tmp/{/.} && \
             exec_annotation \
                -o kofamscan_${label}_annotation/{/.}.out \
                -p /kofamscan/db/profiles/ncycle.hal \
                --tmp-dir tmp/{/.} \
                --cpu 1 {}'
  filter_kofamscan.py \
    -i kofamscan_${label}_annotation \
    -o ${sample_id}_${label}_kofamscan_filtered.tsv \
    --min_tier putative \
    --evalue_cutoff 1e-5 \
    --min_putative_score 50
    """
}
