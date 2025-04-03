process prokka {
    tag "PROKKA on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/prokka", mode: 'copy'

    container "/home/luiszuniga/work/containers/prokka.sif"

    input:
    val(sample_id)
    file(contigs)
    file(metabat_out)
    file(concoct_out)
    file(checkm2_out)

    output:
    file("prokka_assembly_annotation")
    file("prokka_concoct_annotation")
    file("prokka_metabat_annotation")

    script:

    """

    awk '{ print \$1 }' ${checkm2_out}/concoct_bins/concoct_sorted_quality_report.tsv >> concoct_bin_list.txt
    sed -i "1d" concoct_bin_list.txt
    awk '{ print \$1 }' ${checkm2_out}/metabat_bins/metabat_sorted_quality_report.tsv >> metabat_bin_list.txt
    sed -i "1d" metabat_bin_list.txt
    
    mkdir -p prokka_assembly_annotation
    mkdir -p prokka_metabat_annotation
    mkdir -p prokka_concoct_annotation

    echo "Checking MetaBAT input files:"
    ls -l ${metabat_out}/*.fa

    echo "Checking CONCOCT input files:"
    ls -l ${concoct_out}/fasta_bins/*.fa


    # Annotate contigs (single genome assembly)
    prokka --force --outdir prokka_assembly_annotation --prefix ${sample_id} ${contigs}/contigs.fasta
    # Annotate MetaBAT bins in parallel
    ls ${metabat_out}/*.fa | parallel -j 4 'prokka --outdir prokka_metabat_annotation/{/.} --prefix {/.} {}'
    # Annotate CONCOCT bins in parallel
    ls ${concoct_out}/fasta_bins/*.fa | parallel -j 4 'prokka --outdir prokka_concoct_annotation/{/.} --prefix {/.} {}'

    /home/luiszuniga/self/MAG/bin/./rename_bins.sh

    """
}
