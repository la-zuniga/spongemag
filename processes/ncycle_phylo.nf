// processes/ncycle_phylo.nf
process ncycle_phylo {
    tag "N-cycle phylogenies: $sample_id"
    publishDir "${params.outdir}/${sample_id}/ncycle_phylo", mode: 'copy'

    input:
    val(sample_id)
    path(bins_kofamscan_filtered)
    path(assembly_kofamscan_filtered, stageAs: 'assembly_kofamscan_filtered.tsv')
    path(prokka_bins_annotation)
    path(prokka_assembly_annotation)

    output:
    path("sequences")
    path("alignments")
    path("trees")

    script:
    """
    # 1. Extract per-gene amino acid FASTAs from both bins and full assembly
    singularity exec ${params.containers.python} \
        extract_ncycle_seqs.py \
        --bins_tsv ${bins_kofamscan_filtered} \
        --bins_faa ${prokka_bins_annotation} \
        --assembly_tsv ${assembly_kofamscan_filtered} \
        --assembly_faa ${prokka_assembly_annotation} \
        -o sequences \
        --key_only

    mkdir -p alignments trees

    # 2. Align each gene with MAFFT, then build tree with FastTree
    for faa in sequences/*.faa; do
        gene=\$(basename "\$faa" .faa)

        singularity exec ${params.containers.mafft} \
            mafft --auto "\$faa" > "alignments/\${gene}.aln"

        singularity exec ${params.containers.fasttree} \
            FastTree "alignments/\${gene}.aln" > "trees/\${gene}.nwk"

        echo "  \${gene}: alignment + tree done"
    done

    # Handle case where no genes had >= 2 sequences
    if [ -z "\$(ls -A sequences/*.faa 2>/dev/null)" ]; then
        echo "No genes with >= 2 sequences found — skipping alignment/tree steps"
        touch alignments/.empty trees/.empty
    fi
    """
}
