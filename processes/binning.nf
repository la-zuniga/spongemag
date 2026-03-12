process Binning {
    tag "Binning $sample_id"
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
    file("${sample_id}_concoct_contig_bin.tsv")
    file("${sample_id}_metabat_contig_bin.tsv")
    file("${sample_id}_contigs_10K.fa")
    file("metabinner_out")

    script:
    """
    # Cut up contigs for CONCOCT
    singularity exec ${params.containers.concoct} \
        cut_up_fasta.py ${contigs}/contigs.fasta -c 10000 -o 0 --merge_last \
        -b ${sample_id}_contigs_10K.bed > ${sample_id}_contigs_10K.fa

    # Filter short contigs and generate kmer profile (MetaBinner)
    singularity exec ${params.containers.metabinner} \
        Filter_tooshort.py ${sample_id}_contigs_10K.fa 1000
    singularity exec ${params.containers.metabinner} \
        gen_kmer.py ${sample_id}_contigs_10K.fa 1000 4

    # Generate coverage table (CONCOCT)
    singularity exec ${params.containers.concoct} \
        concoct_coverage_table.py ${sample_id}_contigs_10K.bed ${sample_id}.sorted.bam \
        > ${sample_id}_coverage_table.tsv

    # Sort coverage table using helper script (from bin/)
    sort_coverage.py \
        -c ${sample_id}_coverage_table.tsv \
        -k ${sample_id}_contigs_10K_kmer_4_f1000.csv \
        -o ${sample_id}_coverage_table_sorted.tsv

    mkdir -p concoct_out concoct_out/fasta_bins checkm2_out

    # Run CONCOCT
    singularity exec ${params.containers.concoct} \
        concoct \
        --composition_file ${sample_id}_contigs_10K_1000.fa \
        --coverage_file ${sample_id}_coverage_table_sorted.tsv \
        -b concoct_out \
        --threads ${params.concoct_threads} \
        --length_threshold 1000 \
        --read_length ${params.concoct_read_len} \
        --clusters ${params.concoct_clusters}

    singularity exec ${params.containers.concoct} \
        merge_cutup_clustering.py concoct_out/clustering_gt1000.csv \
        > concoct_out/clustering_merged.csv

    singularity exec ${params.containers.concoct} \
        extract_fasta_bins.py ${contigs}/contigs.fasta \
        concoct_out/clustering_merged.csv \
        --output_path concoct_out/fasta_bins

    # Run MetaBAT2
    singularity exec ${params.containers.metabat} \
        jgi_summarize_bam_contig_depths \
        --outputDepth ${sample_id}_depth.txt ${sample_id}.sorted.bam
    singularity exec ${params.containers.metabat} \
        metabat2 -i ${contigs}/contigs.fasta -a ${sample_id}_depth.txt -o metabat_out/

    # Rename hidden MetaBAT output files
    cd metabat_out
    for file in .*.fa; do
        cp "\$file" "\$(basename "\$file" | sed 's/^\\.//')"
    done
    rm -f .*.fa
    cd ..

    # Run MetaBinner
    singularity exec ${params.containers.metabinner} \
        run_metabinner.sh \
        -a ${sample_id}_contigs_10K.fa \
        -o metabinner_out \
        -d ${sample_id}_coverage_table_sorted.tsv \
        -k ${sample_id}_contigs_10K_kmer_4_f1000.csv \
        -p /opt/MetaBinner \
        -t ${params.metabinner_threads}

    # Generate contig-to-bin TSV for MetaBAT (from bin/)
    metabat_contig2bin.py metabat_out/ ${sample_id}_metabat_contig_bin.tsv


    # Generate CONCOCT contig-to-bin TSV
    singularity exec ${params.containers.concoct} \
        sed 's/,/\\t/g' concoct_out/clustering_merged.csv \
        | tail -n +2 \
        | awk -F'\\t' '{print \$1 "\\t" \$2 "_concoct"}' \
        > ${sample_id}_concoct_contig_bin.tsv
    """
}
