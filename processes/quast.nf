process quast {
    tag "QUAST on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/quast", mode: 'copy'

    container "/home/luiszuniga/work/containers/QUAST_v5.3.0.sif"

    input:
    val(sample_id)
    path(metabat_path)
    path(concoct_path)


    output:
    file("quast_output/metabat_bins")
    file("quast_output/concoct_bins")

    script:
    """
    # Running QUAST on the MetaBAT bins

    echo "Concoct Path: ${concoct_path}"
    echo "Metabat Path: ${metabat_path}"

    quast -o quast_output/concoct_bins ${concoct_path}/fasta_bins/*.fa
    
    quast -o quast_output/metabat_bins ${metabat_path}/*.fa
    """
}