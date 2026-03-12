process KOfamscan {
    tag "KOfamscan on $sample_id metagenomes (translated)"
    publishDir "${params.outdir}/${sample_id}/kofamscan", mode: 'copy'

    container "${params.containers.kofamscan}"

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
    mkdir -p tmp/tabular
    mkdir -p tmp/mapper
    echo "Current working directory: \$(pwd)"
    echo "Listing contents of prokka_metabat_annotation:"
    ls -lh ${prokka_metabat_annotation}


    find ${prokka_metabat_annotation}/ -type f -name "*.faa" | parallel -j 2 'exec_annotation -o kofamscan_metabat_annotation/{/.}.out -p /kofamscan/db/profiles/ncycle.hal {}'

    find ${prokka_concoct_annotation}/ -type f -name "*.faa" | parallel -j 2 'exec_annotation -o kofamscan_concoct_annotation/{/.}.out -p /kofamscan/db/profiles/ncycle.hal {}'


    """
}
