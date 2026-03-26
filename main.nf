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
include { binning_prep; binning_concoct; binning_metabat; binning_metabinner } from './processes/binning.nf'
include { checkm2 }     from './processes/checkm2.nf'
include { prokka }      from './processes/prokka.nf'
include { KOfamscan }                        from './processes/kofamscan.nf'
include { KOfamscan as KOfamscan_assembly }  from './processes/kofamscan.nf'
include { quast }       from './processes/quast.nf'
include { dastool }  from './processes/dastool.nf'
include { kegg_pathway } from './processes/kegg_pathway.nf'
include { ncycle_phylo } from './processes/ncycle_phylo.nf'
// include { gtdbtk }   from './processes/gtdbtk.nf'
// I'm working on trying to immplement DAStool to get more robust consensus bins, but haven't worked it out yet, so it's commented out.
// the gtdbtk-container I've built is quite large, but I know it works, so I comment it out when testing main.nf
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

    // Step 2: Assembly
    if (params.assembly == 'megahit') {
        assembly = megahit(fastpOutput[0])
    } else if (params.assembly == 'spades') {
        assembly = spades(fastpOutput[0])
    } else {
        error "Invalid assembler choice: ${params.assembly}"
    }

    // Step 3a: Align reads to assembly (minimap2 | samtools — no SAM on disk)
    alignOutput = align(fastpOutput[0], assembly)

    // Step 3b: Shared binning prep (cuts contigs, builds coverage + kmer profile)
    prepOut = binning_prep(alignOutput, assembly)

    // Step 3c: Run three binners in parallel
    //   prepOut tuple: [sample_id, contigs_10K_fa, contigs_10K_bed,
    //                   contigs_10K_1000_fa, kmer_csv, coverage_table, coverage_sorted]
    concoctOut    = binning_concoct(
        prepOut.map { sid, fa10K, bed, fa1000, kmer, cov, cov_sorted -> tuple(sid, fa1000, cov_sorted) },
        assembly
    )
    metabatOut    = binning_metabat(alignOutput, assembly)
    metabinnerOut = binning_metabinner(
        prepOut.map { sid, fa10K, bed, fa1000, kmer, cov, cov_sorted -> tuple(sid, fa10K, cov_sorted, kmer) }
    )

    // Step 3d: Join binner outputs by sample_id, then run DASTool consensus binning
    // joined_bins: [sample_id, concoct_tsv, metabat_tsv, metabinner_tsv]
    joined_bins = concoctOut[0].join(metabatOut[0]).join(metabinnerOut[0])

    // dastoolOutput[0] = sample_id
    // dastoolOutput[1] = DAStool_out (full output dir)
    // dastoolOutput[2] = DASTool_bins (refined bins dir — fed downstream)
    dastoolOutput  = dastool(joined_bins, assembly)

    // Step 4a: Quality Assessment on refined bins only
    // checkm2Output[0] = sample_id
    // checkm2Output[1] = checkm2_out
    checkm2Output  = checkm2(dastoolOutput[0], dastoolOutput[2])

    // Step 4: Annotation of quality-filtered refined bins
    // protein_annotation[0] = prokka_assembly_annotation
    // protein_annotation[1] = prokka_bins_annotation
    protein_annotation = prokka(
        dastoolOutput[0],
        assembly,
        dastoolOutput[2],
        checkm2Output[1]
    )
    kofam_annotation = KOfamscan(
        dastoolOutput[0],
        protein_annotation[1]
    )

    // Step 4b: KOfamscan on full assembly (captures unbinned contigs too)
    kofam_assembly = KOfamscan_assembly(
        dastoolOutput[0],
        protein_annotation[0]
    )

    // Step 5: Assembly QC on refined bins
    quastReport = quast(dastoolOutput[0], dastoolOutput[2])

    kegg_results = kegg_pathway(
    dastoolOutput[0],
    kofam_annotation[1]    // kofamscan_filtered.tsv
    )

    // Step 6: Nitrogen cycle gene phylogenies (bins + assembly)
    ncycle_phylo(
        dastoolOutput[0],
        kofam_annotation[1],       // bins filtered TSV
        kofam_assembly[1],         // assembly filtered TSV
        protein_annotation[1],     // prokka_bins_annotation
        protein_annotation[0]      // prokka_assembly_annotation
    )

    // Optional
    // gtdbtk_classification = gtdbtk(dastoolOutput[2])}
  }















































