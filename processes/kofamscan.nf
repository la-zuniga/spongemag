process KOfamscan {
    tag "KOfamscan on $sample_id metagenomes (translated)"
    publishDir "${params.outdir}/${sample_id}/kofamscan", mode: 'copy'

    container "/home/luiszuniga/work/containers/KOfamscan_v1.3.0.sif"

    input:
    val(sample_id)
    file(prokka_concoct_annotation)
    file(prokka_metabat_annotation)

    output:
    file("kofamscan_concoct_annotation")
    file("kofamscan_metabat_annotation")

    script:

    """
    mkdir kofamscan_concoct_annotation
    mkdir kofamscan_metabat_annotation


    find ${prokka_concoct_annotation}/ -type f -name "*.faa" | parallel -j 4 'exec_annotation -o kofamscan_concoct_annotation/{/.}.out -p /kofamscan/db/profiles/ncycle.hal {}'
    find ${prokka_metabat_annotation}/ -type f -name "*.faa" | parallel -j 4 'exec_annotation -o kofamscan_metabat_annotation/{/.}.out -p /kofamscan/db/profiles/ncycle.hal {}'


    """
}
