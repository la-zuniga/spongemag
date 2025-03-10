process KOfamscan {
    tag "KOFamScan on $sample_id metagenomes (translated)"
    publishDir "${params.outdir}/${sample_id}/kofamscan", mode: 'copy'

    container "/home/luiszuniga/work/containers/KOfamscan.sif"

    input:
    val(sample_id)
    file(prokka_metabat_annotation)
    file(prokka_concoct_annotation)

    output:
    file("prokka_bin_annotation")

    script:
    """
    mkdir prokka_annotation
    
    prokka --outdir prokka_bin_annotation --prefix ${sample_id} ${contigs}/contigs.fasta

    """
}
