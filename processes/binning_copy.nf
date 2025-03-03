// processes/binning.nf
process MetaBinner {
    tag "MetaBinner on $sample_id"
    publishDir "${params.outdir}/${sample_id}/binning", mode: 'copy'

    input:
    tuple val(sample_id), file("${sample_id}.bam"), file("${sample_id}.sorted.bam"), file("${sample_id}.sorted.bai")
    file(contigs)

    output:
    val(sample_id)
    file("${sample_id}_coverage_table.tsv")
    file("${sample_id}_contigs_10K_kmer_4_f1000.csv")
    file("${sample_id}_coverage_table_sorted.tsv")
    file("concoct_out")
    file("metabat_out")
    file("checkm2_out")
    // file("${sample_id}_concoct_contig_bin.tsv")
    // file("${sample_id}_metabat_contig_bin.tsv")

    script:
    """

    # Run Concoct: Cut up contigs 
    singularity exec /home/luiszuniga/work/containers/concoct_v1.1.0.sif cut_up_fasta.py ${contigs}/contigs.fasta -c 10000 -o 0 --merge_last -b ${sample_id}_contigs_10K.bed > ${sample_id}_contigs_10K.fa

    # Switch to MetaBinner container to filter short contigs and generate kmer profile
    singularity exec /home/luiszuniga/work/containers/metabinner_v1.4.4.sif Filter_tooshort.py ${sample_id}_contigs_10K.fa 1000
    singularity exec /home/luiszuniga/work/containers/metabinner_v1.4.4.sif gen_kmer.py ${sample_id}_contigs_10K.fa 1000 4

    # Concoct: Generate coverage table
    singularity exec /home/luiszuniga/work/containers/concoct_v1.1.0.sif concoct_coverage_table.py ${sample_id}_contigs_10K.bed ${sample_id}.sorted.bam > ${sample_id}_coverage_table.tsv
    # sorting the concoct coverage file with a custom script
    python3 /home/luiszuniga/work/python_scripts/sortCoverage/sort_coverage.py -c ${sample_id}_coverage_table.tsv -k ${sample_id}_contigs_10K_kmer_4_f1000.csv -o ${sample_id}_coverage_table_sorted.tsv

    mkdir concoct_out
    mkdir checkm2_out

    # running concoct and extracting bins
    singularity exec /home/luiszuniga/work/containers/concoct_v1.1.0.sif concoct --composition_file ${sample_id}_contigs_10K_1000.fa --coverage_file ${sample_id}_coverage_table_sorted.tsv -b concoct_out --threads 26 --length_threshold 1000 --read_length 250 --clusters 5

    singularity exec /home/luiszuniga/work/containers/concoct_v1.1.0.sif merge_cutup_clustering.py concoct_out/clustering_gt1000.csv > concoct_out/clustering_merged.csv
    # making the contig bin file
    concoct_contig2bin.sh concoct_out/clustering_gt1000.csv ${sample_id}_concoct_contig_bin.tsv
    singularity exec /home/luiszuniga/work/containers/concoct_v1.1.0.sif mkdir concoct_out/fasta_bins

    singularity exec /home/luiszuniga/work/containers/concoct_v1.1.0.sif extract_fasta_bins.py ${contigs}/contigs.fasta concoct_out/clustering_merged.csv --output_path concoct_out/fasta_bins

    # running metabat

    singularity exec /home/luiszuniga/work/containers/metabat_v2.17.sif jgi_summarize_bam_contig_depths --outputDepth ${sample_id}_depth.txt ${sample_id}.sorted.bam
    singularity exec /home/luiszuniga/work/containers/metabat_v2.17.sif metabat2 -i ${contigs}/contigs.fasta -a ${sample_id}_depth.txt -o metabat_out/
    
    # renaming bin files so they are not hidden files
    cd metabat_out
    for file in .*.fa; do
        cp "\$file" "\$(basename "\$file" | sed 's/^\\.//')"
    done

    # removing hidden bins
    rm .*.fa
    cd ..
    # generating the contig bin file for metabat bins
    python3 /home/luiszuniga/self/MAG/bin/metabat_contig2bin.py metabat_out/ metabat${sample_id}_metabat_contig_bin.tsv

    # running checkm2 on the bins from concoct and metabat
    cd ..
    singularity exec /home/luiszuniga/work/containers/checkm2_v1.0.2.sif checkm2 predict -x .fa --threads 26 --input concoct_out/fasta_bins --output-directory checkm2_out/concoct_bins --database_path /db/uniref100.KO.1.dmnd
    singularity exec /home/luiszuniga/work/containers/checkm2_v1.0.2.sif checkm2 predict -x .fa --threads 26 --input metabat_out/ --output-directory checkm2_out/metabat_bins --database_path /db/uniref100.KO.1.dmnd
    """
}
