process prokka {
    tag "PROKKA on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/prokka", mode: 'copy'

    container "/home/luiszuniga/work/containers/prokka.sif"

    input:
    val(sample_id)
    file(contigs)
    file(metabat_out)
    file(concoct_out)

    output:
    file("prokka_assembly_annotation")
    file("prokka_concoct_annotation")
    file("prokka_metabat_annotation")

    script:

    """

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

    for folder in prokka_concoct_annotation/* prokka_metabat_annotation/*; do
        if [ -d "$folder" ]; then
            mag_name=$$(basename "$folder")
            tool=$$(basename "$(dirname "$folder")" | cut -d'_' -f2)  # Extracts 'concoct' or 'metabat'

            # Rename each file in the folder
            for file in "$folder"/PROKKA_*; do
                ext="${file##*.}"  # Get file extension
                mv "$file" "$folder/${mag_name}_${tool}.${ext}"
            done
        fi
    done


    """
}
