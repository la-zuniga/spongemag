nextflow.enable.dsl=2

// =============================
// Project Directories & Params
// =============================
projectDir = "/home/luis/spongemag"
params.input_dir = "/home/luis/data"
params.reads = "$params.input_dir/*_{R1,R2}.fastq.gz"
params.outdir = "$projectDir/out"
containers = "/home/luis/containers"
bin = "$projectDir/bin"

// =============================
// Include Process Modules
// =============================
include { fastp }       from './processes/fastp.nf'
include { fastqc }      from './processes/FastQC.nf'
include { multiqc }     from './processes/multiqc.nf'
include { megahit }     from './processes/megahit.nf'
include { spades }      from './processes/spades.nf'
include { align }       from './processes/align.nf'
include { binning_prep; binning_concoct; binning_metabat; metabinner_prep; binning_metabinner } from './processes/binning.nf'
include { checkm2 }     from './processes/checkm2.nf'
include { bakta_assembly; bakta_bins } from './processes/bakta.nf'
include { KOfamscan as KOfamscan_bins }     from './processes/kofamscan.nf'
include { KOfamscan as KOfamscan_assembly } from './processes/kofamscan.nf'
include { quast }       from './processes/quast.nf'
include { dastool }     from './processes/dastool.nf'
include { kegg_pathway as kegg_pathway_bins }     from './processes/kegg_pathway.nf'
include { kegg_pathway as kegg_pathway_assembly } from './processes/kegg_pathway.nf'
include { ncycle_phylo_bins;
          ncycle_phylo_assembly;
          ncycle_phylo_cross_bins;
          ncycle_phylo_cross_assembly } from './processes/ncycle_phylo.nf'

// =============================
// Input Channels
// =============================
Channel.fromFilePairs(params.reads).set { read_pairs_ch }

// =============================
// Workflow Definition
// =============================
workflow {

    // Step 1: Preprocessing
    fastpOutput   = fastp(read_pairs_ch)
    fastQCoutput  = fastqc(fastpOutput[0])
    multiqcReport = multiqc(fastQCoutput[0])

    // Step 2: Assembly — emits [sample_id, contigs_dir]
    if (params.assembly == 'megahit') {
        assembly = megahit(fastpOutput[0])
    } else if (params.assembly == 'spades') {
        assembly = spades(fastpOutput[0])
    } else {
        error "Invalid assembler choice: ${params.assembly}"
    }

    // Step 3a: Align reads to assembly (minimap2 | samtools — no SAM on disk)
    alignOutput = align(fastpOutput[0].join(assembly))

    // Step 3b: Shared binning prep (cuts contigs, builds coverage + kmer profile)
    prepOut = binning_prep(alignOutput.join(assembly))

    // Step 3c: Run three binners in parallel
    concoctOut    = binning_concoct(
        prepOut.map { sid, fa10K, bed, fa1000, kmer, cov, cov_sorted -> tuple(sid, fa1000, cov_sorted) }
            .join(assembly)
    )
    metabatOut    = binning_metabat(alignOutput.join(assembly))

    mbPrepOut = metabinner_prep(alignOutput.join(assembly))
    metabinnerOut = binning_metabinner(
        mbPrepOut.map { sid, contigs, kmer, cov -> tuple(sid, contigs, cov, kmer) }
    )

    // Step 3d: DAS Tool consensus
    joined_bins = concoctOut[0].join(metabatOut[0]).join(metabinnerOut[0])
    dastoolOutput = dastool(joined_bins.join(assembly))
    // dastoolOutput[0] = sample_id val channel
    // dastoolOutput[1] = DAStool_out dir
    // dastoolOutput[2] = DASTool_bins dir

    // Step 4a: CheckM2 (emits tier TSVs inside checkm2_out/)
    checkm2Output = checkm2(dastoolOutput[0], dastoolOutput[2])
    // checkm2Output[0] = sample_id val channel
    // checkm2Output[1] = checkm2_out dir

    // ── Build keyed [sid, artifact] channels for joining ──────────────────
    dastool_bins_keyed = dastoolOutput[0].merge(dastoolOutput[2])
    checkm2_keyed      = checkm2Output[0].merge(checkm2Output[1])

    // ── Tier fan-out for bin-path processes ───────────────────────────────
    tier_ch = Channel.of('high', 'medium')

    // [sid, dastool_bins, checkm2_out, tier] → [sid, tier, dastool_bins, checkm2_out]
    bins_base = dastool_bins_keyed
        .join(checkm2_keyed)
        .combine(tier_ch)
        .map { sid, db, co, tier -> tuple(sid, tier, db, co) }

    // Step 4b: Full-assembly Bakta — runs once per sample, gated on dastool
    assembly_for_bakta = dastoolOutput[0]
        .map { sid -> [sid, sid] }
        .join(assembly)
        .map { sid, dummy, contigs -> tuple(sid, contigs) }
    baktaAssemblyOut = bakta_assembly(assembly_for_bakta)
    // baktaAssemblyOut: [sid, bakta_assembly_annotation]

    // Step 4c: Per-tier bin Bakta
    baktaBinsOut = bakta_bins(bins_base)
    // baktaBinsOut: [sid, tier, bakta_bins_annotation]

    // Step 4d: KOfamscan — bins (per tier) + assembly (once)
    kofamBinsOut = KOfamscan_bins(
        baktaBinsOut.map { sid, tier, faa -> tuple(sid, "bins_${tier}", faa) }
    )
    // kofamBinsOut: [sid, "bins_${tier}", annotation_dir, filtered_tsv]

    kofamAssemblyOut = KOfamscan_assembly(
        baktaAssemblyOut.map { sid, faa -> tuple(sid, 'assembly', faa) }
    )
    // kofamAssemblyOut: [sid, 'assembly', annotation_dir, filtered_tsv]

    // Step 5: QUAST per tier
    quast(bins_base)

    // Step 6a: KEGG pathway — bins (per tier) + assembly (once)
    kegg_pathway_bins(
        kofamBinsOut.map { sid, label, ann, tsv -> tuple(sid, label, tsv) }
    )
    kegg_pathway_assembly(
        kofamAssemblyOut.map { sid, label, ann, tsv -> tuple(sid, label, tsv) }
    )

    // Step 6b: Per-sample N-cycle phylogenies — bins (per tier) + assembly (once)
    // Join kofam filtered TSV with prokka annotation dir on (sid, tier) for bins
    phylo_bins_in = kofamBinsOut
        .map { sid, label, ann, tsv ->
            def tier = label - 'bins_'
            tuple(sid, tier, tsv)
        }
        .join(
            baktaBinsOut.map { sid, tier, faa -> tuple(sid, tier, faa) },
            by: [0, 1]
        )
        // [sid, tier, bins_tsv, bakta_bins_dir]
    ncycleBinsOut = ncycle_phylo_bins(phylo_bins_in)
    // ncycleBinsOut: [sid, tier, sequences_dir, alignments, trees]

    // Assembly phylo: join on sid
    phylo_assembly_in = baktaAssemblyOut
        .join(kofamAssemblyOut.map { sid, label, ann, tsv -> tuple(sid, tsv) })
        .map { sid, faa, tsv -> tuple(sid, tsv, faa) }
    ncycleAssemblyOut = ncycle_phylo_assembly(phylo_assembly_in)
    // ncycleAssemblyOut: [sid, sequences_dir, alignments, trees]

    // Step 6c: Cross-sample N-cycle phylogenies
    // Bins: group sequences dirs by tier across all samples
    cross_bins_in = ncycleBinsOut
        .map { sid, tier, seqs, aln, trees -> tuple(tier, seqs) }
        .groupTuple()
    ncycle_phylo_cross_bins(cross_bins_in)

    // Assembly: collect sequences dirs across all samples
    ncycle_phylo_cross_assembly(
        ncycleAssemblyOut
            .map { sid, seqs, aln, trees -> seqs }
            .collect()
    )

    // Optional
    // gtdbtk_classification = gtdbtk(dastoolOutput[2])
}
