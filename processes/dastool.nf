process dastool {
    tag "DASTool on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/DAStool", mode: 'copy'

    container "${params.containers.dastool}"

    input:
    tuple val(sample_id), path(concoct_contig_bin), path(metabat_contig_bin), path(metabinner_contig_bin)
    path(contigs)

    output:
    val(sample_id)
    file("DAStool_out")
    path("DAStool_out/${sample_id}_DASTool_bins")

    script:
    """
    mod_headers.py ${contigs}/contigs.fasta contigs.fasta
    
    mkdir -p DAStool_out tmp_dastool
    DAS_Tool  -i ${concoct_contig_bin},${metabat_contig_bin},${metabinner_contig_bin} \
    -l concoct,metabat,metabinner \
    -c contigs.fasta \
    -o DAStool_out/${sample_id} \
    --write_bins \
    --debug 
    """
}
