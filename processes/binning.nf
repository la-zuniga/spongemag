// processes/binning.nf
//
// Split into 4 processes so the three binners can run in parallel:
//   binning_prep      — shared preprocessing (cut contigs, coverage table, kmer profile)
//   binning_concoct   — CONCOCT
//   binning_metabat   — MetaBAT2
//   binning_metabinner — MetaBinner


// ── Shared preprocessing ──────────────────────────────────────────────────────
// Produces the files that CONCOCT and MetaBinner both need.
// Requires two containers (concoct + metabinner), so singularity exec is called
// manually — no container directive here.
process binning_prep {
    tag "Binning prep: $sample_id"
    publishDir "${params.outdir}/${sample_id}/binning/prep", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai), path(contigs)

    output:
    tuple val(sample_id),
          path("${sample_id}_contigs_10K.fa"),
          path("${sample_id}_contigs_10K.bed"),
          path("${sample_id}_contigs_10K_1000.fa"),
          path("${sample_id}_contigs_10K_kmer_4_f1000.csv"),
          path("${sample_id}_coverage_table.tsv"),
          path("${sample_id}_coverage_table_sorted.tsv")

    script:
    """
    # Cut up contigs for CONCOCT (produces .fa + .bed)
    singularity exec ${params.containers.concoct} \
        cut_up_fasta.py ${contigs}/contigs.fasta -c 10000 -o 0 --merge_last \
        -b ${sample_id}_contigs_10K.bed > ${sample_id}_contigs_10K.fa

    # Filter short contigs (produces ${sample_id}_contigs_10K_1000.fa)
    singularity exec ${params.containers.metabinner} \
        Filter_tooshort.py ${sample_id}_contigs_10K.fa 1000

    # Generate kmer profile for MetaBinner
    singularity exec ${params.containers.metabinner} \
        gen_kmer.py ${sample_id}_contigs_10K.fa 1000 4

    # Generate coverage table from sorted BAM
    singularity exec ${params.containers.concoct} \
        concoct_coverage_table.py ${sample_id}_contigs_10K.bed ${sorted_bam} \
        > ${sample_id}_coverage_table.tsv

    # Filter coverage table to only contigs present in the kmer profile
    sort_coverage.py \
        -c ${sample_id}_coverage_table.tsv \
        -k ${sample_id}_contigs_10K_kmer_4_f1000.csv \
        -o ${sample_id}_coverage_table_sorted.tsv
    """
}


// ── CONCOCT ───────────────────────────────────────────────────────────────────
process binning_concoct {
    tag "CONCOCT: $sample_id"
    publishDir "${params.outdir}/${sample_id}/binning/concoct", mode: 'copy'

    container "${params.containers.concoct}"

    input:
    tuple val(sample_id), path(contigs_10K_1000), path(coverage_sorted), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_concoct_contig_bin.tsv")
    path("concoct_out")

    script:
    """
    mkdir -p concoct_out/fasta_bins

    concoct \
        --composition_file ${contigs_10K_1000} \
        --coverage_file ${coverage_sorted} \
        -b concoct_out \
        --threads ${params.concoct_threads} \
        --length_threshold 1000 \
        --read_length ${params.concoct_read_len} \
        --clusters ${params.concoct_clusters}

    merge_cutup_clustering.py concoct_out/clustering_gt1000.csv \
        > concoct_out/clustering_merged.csv

    extract_fasta_bins.py ${contigs}/contigs.fasta \
        concoct_out/clustering_merged.csv \
        --output_path concoct_out/fasta_bins

    # Generate contig-to-bin TSV
    sed 's/,/\\t/g' concoct_out/clustering_merged.csv \
        | tail -n +2 \
        | awk -F'\\t' '{print \$1 "\\t" \$2 "_concoct"}' \
        > ${sample_id}_concoct_contig_bin.tsv
    """
}


// ── MetaBAT2 ──────────────────────────────────────────────────────────────────
process binning_metabat {
    tag "MetaBAT2: $sample_id"
    publishDir "${params.outdir}/${sample_id}/binning/metabat", mode: 'copy'

    container "${params.containers.metabat}"

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_metabat_contig_bin.tsv")
    path("metabat_out")

    script:
    """
    mkdir -p metabat_out

    jgi_summarize_bam_contig_depths \
        --outputDepth ${sample_id}_depth.txt ${sorted_bam}

    metabat2 -i ${contigs}/contigs.fasta \
             -a ${sample_id}_depth.txt \
             -o metabat_out/

    # MetaBAT2 sometimes hides output files with a leading dot — unhide them
    for file in metabat_out/.*.fa; do
        [ -f "\$file" ] || continue
        newname=\$(basename "\$file")
        newname=\${newname#.}
        cp "\$file" "metabat_out/\${newname}"
    done
    rm -f metabat_out/.*.fa

    metabat_contig2bin.py metabat_out/ ${sample_id}_metabat_contig_bin.tsv
    """
}


// ── MetaBinner prep ──────────────────────────────────────────────────────────
// Generates kmer profile and coverage table from the ORIGINAL assembly contigs
// (not the 10kb cut-ups that CONCOCT uses). MetaBinner was designed for
// original contigs — feeding it pre-fragmented contigs hurts its clustering.
// Requires two containers (metabinner + metabat), so no container directive.
process metabinner_prep {
    tag "MetaBinner prep: $sample_id"
    publishDir "${params.outdir}/${sample_id}/binning/metabinner_prep", mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam), path(sorted_bai), path(contigs)

    output:
    tuple val(sample_id),
          path("${sample_id}_contigs.fasta"),
          path("${sample_id}_contigs_kmer_4_f1000.csv"),
          path("${sample_id}_coverage_sorted.tsv")

    script:
    """
    # Copy original contigs and normalize headers to ID-only
    # (strip anything after first whitespace — MetaBinner's split_hhbins.py 
    # expects the full header to match the contig ID in bin files)
    awk '/^>/ {print \$1; next} {print}' ${contigs}/contigs.fasta \\
        > ${sample_id}_contigs.fasta

    # Generate kmer profile from normalized contigs
    singularity exec ${params.containers.metabinner} \\
        gen_kmer.py ${sample_id}_contigs.fasta 1000 4

    # Generate per-contig coverage from BAM
    singularity exec ${params.containers.metabat} \\
        jgi_summarize_bam_contig_depths \\
            --outputDepth ${sample_id}_jgi_depth.txt ${sorted_bam}

    # Extract contig name + mean depth (columns 1 and 3, skip header)
    awk -F'\\t' 'NR>1 {print \$1 "\\t" \$3}' ${sample_id}_jgi_depth.txt \\
        > ${sample_id}_coverage.tsv

    # Filter coverage to contigs present in kmer profile
    sort_coverage.py \\
        -c ${sample_id}_coverage.tsv \\
        -k ${sample_id}_contigs_kmer_4_f1000.csv \\
        -o ${sample_id}_coverage_sorted.tsv
    """
     
    
}


// ── MetaBinner ────────────────────────────────────────────────────────────────
process binning_metabinner {
    tag "MetaBinner: $sample_id"
    publishDir "${params.outdir}/${sample_id}/binning/metabinner", mode: 'copy'

    container "${params.containers.metabinner}"
    // Bind-mount patched split_hhbins.py over the container's copy to fix a
    // KeyError on the 2nd post-process pass (see bin/split_hhbins.py header).
    containerOptions "--bind ${projectDir}/bin/split_hhbins.py:/opt/MetaBinner/scripts/split_hhbins.py"

    input:
    tuple val(sample_id), path(contigs), path(coverage_sorted), path(kmer_csv)

    output:
    tuple val(sample_id), path("${sample_id}_metabinner_contig_bin.tsv")
    path("metabinner_out")

    script:
    """
    run_metabinner.sh \
        -a ${contigs} \
        -o metabinner_out \
        -d ${coverage_sorted} \
        -k ${kmer_csv} \
        -p /opt/MetaBinner \
        -t ${params.metabinner_threads}

    METABINNER_BINS=\$(ls -d metabinner_out/metabinner_res/unitem_profile/kmeans_length_weight_*_t_logtrans_result.tsv_bins)

    metabinner_contig2bin.py \
        \${METABINNER_BINS} \
        ${sample_id}_metabinner_contig_bin.tsv
    """
}
