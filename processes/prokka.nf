process prokka {
    tag "PROKKA on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/prokka", mode: 'copy'

    container "${params.containers.prokka}"

    input:
    val(sample_id)
    file(contigs)
    file(dastool_bins)
    file(checkm2_out)

    output:
    file("prokka_assembly_annotation")
    file("prokka_bins_annotation")
    script:

    """
    awk '{ print \$1 }' ${checkm2_out}/sorted_quality_report.tsv \
        | tail -n +2 > bin_list.txt

    mkdir -p prokka_assembly_annotation prokka_bins_annotation

    # Annotate full assembly
    prokka --force \
        --outdir prokka_assembly_annotation \
        --prefix ${sample_id} \
        ${contigs}/contigs.fasta

    # Annotate only bins passing CheckM2 quality thresholds
    while read bin; do
        prokka --outdir prokka_bins_annotation/\${bin} \
               --prefix \${bin} \
               ${dastool_bins}/\${bin}.fa
    done < bin_list.txt

    /home/luiszuniga/self/MAG/bin/./rename_bins.sh

    """
}
