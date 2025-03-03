process dastool {
    tag "DASTool on $sample_id metagenomes"
    publishDir "${params.outdir}/${sample_id}/DAStool", mode: 'copy'

    container "/home/luiszuniga/work/containers/DAStool_v1.1.7.sif"

    input:
    val(sample_id)
    file(contigs)
    file("${sample_id}_concoct_contig_bin.tsv")
    file("${sample_id}_metabat_contig_bin.tsv")

    output:
    file("DAStool_out")

    script:
    """
    python3 /home/luiszuniga/self/MAG/bin/mod_headers.py ${contigs}/contigs.fasta contigs.fasta
    
    mkdir DAStool_out
    DAS_Tool  -i ${sample_id}_concoct_contig_bin.tsv,${sample_id}_metabat_contig_bin.tsv,\
    -l concoct,metabat \
    -c contigs.fasta \
    -o ${sample_id} \
    --debug 
    """
}