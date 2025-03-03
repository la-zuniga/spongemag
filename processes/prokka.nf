process prokka {
    tag "PROKKA on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/prokka", mode: 'copy'

    container "/home/luiszuniga/work/containers/prokka.sif"

    input:
    val(sample_id)
    file(contigs)

    output:
    file("prokka_bin_annotation")

    script:
    """
    mkdir prokka_annotation
    
    prokka --outdir prokka_bin_annotation --prefix ${sample_id} ${contigs}/contigs.fasta

    """
}