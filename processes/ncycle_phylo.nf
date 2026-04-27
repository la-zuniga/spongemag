// processes/ncycle_phylo.nf
//
// Split per-sample phylogenies so the bin path can run per quality tier
// while the full-assembly path runs once (tier-agnostic — all contigs are
// the same regardless of which tier of bins we're looking at).
// Cross-sample phylogenies are split the same way.


// ── Per-sample bin phylogenies (per tier) ──────────────────────────────────
process ncycle_phylo_bins {
    tag "N-cycle bins phylo (${tier}): $sample_id"
    publishDir "${params.outdir}/${sample_id}/ncycle_phylo/bins_${tier}", mode: 'copy'

    input:
    tuple val(sample_id), val(tier),
          path(bins_kofamscan_filtered),
          path(prokka_bins_annotation)

    output:
    tuple val(sample_id), val(tier),
          path("${sample_id}_${tier}_sequences"),
          path("alignments"),
          path("trees")

    script:
    """
    extract_ncycle_seqs.py \
        --sample_id ${sample_id} \
        --bins_tsv ${bins_kofamscan_filtered} \
        --bins_faa ${prokka_bins_annotation} \
        -o ${sample_id}_${tier}_sequences \
        --key_only

    mkdir -p alignments trees

    for faa in ${sample_id}_${tier}_sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)

        singularity exec ${params.containers.mafft} \
            mafft --auto "\$faa" > "alignments/\${gene}.aln"

        singularity exec ${params.containers.fasttree} \
            FastTree "alignments/\${gene}.aln" > "trees/\${gene}.nwk"

        echo "  \${gene}: alignment + tree done"
    done

    if [ -z "\$(ls -A ${sample_id}_${tier}_sequences/*.faa 2>/dev/null)" ]; then
        echo "No genes with >= 2 sequences found for ${tier}"
        touch alignments/.empty trees/.empty
    fi
    """
}


// ── Per-sample assembly phylogeny (runs once, tier-agnostic) ──────────────
process ncycle_phylo_assembly {
    tag "N-cycle assembly phylo: $sample_id"
    publishDir "${params.outdir}/${sample_id}/ncycle_phylo/assembly", mode: 'copy'

    input:
    tuple val(sample_id),
          path(assembly_kofamscan_filtered, stageAs: 'assembly_kofamscan_filtered.tsv'),
          path(prokka_assembly_annotation)

    output:
    tuple val(sample_id),
          path("${sample_id}_assembly_sequences"),
          path("alignments"),
          path("trees")

    script:
    """
    extract_ncycle_seqs.py \
        --sample_id ${sample_id} \
        --assembly_tsv ${assembly_kofamscan_filtered} \
        --assembly_faa ${prokka_assembly_annotation} \
        -o ${sample_id}_assembly_sequences \
        --key_only

    mkdir -p alignments trees

    for faa in ${sample_id}_assembly_sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)

        singularity exec ${params.containers.mafft} \
            mafft --auto "\$faa" > "alignments/\${gene}.aln"

        singularity exec ${params.containers.fasttree} \
            FastTree "alignments/\${gene}.aln" > "trees/\${gene}.nwk"

        echo "  \${gene}: alignment + tree done"
    done

    if [ -z "\$(ls -A ${sample_id}_assembly_sequences/*.faa 2>/dev/null)" ]; then
        echo "No genes with >= 2 sequences found for assembly"
        touch alignments/.empty trees/.empty
    fi
    """
}


// ── Cross-sample bin phylogenies (per tier) ────────────────────────────────
process ncycle_phylo_cross_bins {
    tag "Cross-sample N-cycle bins phylo (${tier})"
    publishDir "${params.outdir}/cross_sample/ncycle_phylo/bins_${tier}", mode: 'copy'

    input:
    tuple val(tier), path(seq_dirs)

    output:
    tuple val(tier), path("sequences"), path("alignments"), path("trees")

    script:
    """
    mkdir -p sequences alignments trees

    for faa in *_sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)
        cat "\$faa" >> "sequences/\${gene}.faa"
    done

    for faa in sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)

        singularity exec ${params.containers.mafft} \
            mafft --auto "\$faa" > "alignments/\${gene}.aln"

        singularity exec ${params.containers.fasttree} \
            FastTree "alignments/\${gene}.aln" > "trees/\${gene}.nwk"

        echo "  \${gene}: cross-sample ${tier} alignment + tree done"
    done

    if [ -z "\$(ls -A sequences/*.faa 2>/dev/null)" ]; then
        echo "No genes with >= 2 sequences across ${tier} samples"
        touch alignments/.empty trees/.empty
    fi
    """
}


// ── Cross-sample assembly phylogeny (runs once, tier-agnostic) ────────────
process ncycle_phylo_cross_assembly {
    tag "Cross-sample N-cycle assembly phylo"
    publishDir "${params.outdir}/cross_sample/ncycle_phylo/assembly", mode: 'copy'

    input:
    path(seq_dirs)

    output:
    path("sequences")
    path("alignments")
    path("trees")

    script:
    """
    mkdir -p sequences alignments trees

    for faa in *_sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)
        cat "\$faa" >> "sequences/\${gene}.faa"
    done

    for faa in sequences/*.faa; do
        [ -f "\$faa" ] || continue
        gene=\$(basename "\$faa" .faa)

        singularity exec ${params.containers.mafft} \
            mafft --auto "\$faa" > "alignments/\${gene}.aln"

        singularity exec ${params.containers.fasttree} \
            FastTree "alignments/\${gene}.aln" > "trees/\${gene}.nwk"

        echo "  \${gene}: cross-sample assembly alignment + tree done"
    done

    if [ -z "\$(ls -A sequences/*.faa 2>/dev/null)" ]; then
        echo "No genes with >= 2 sequences across assembly samples"
        touch alignments/.empty trees/.empty
    fi
    """
}
