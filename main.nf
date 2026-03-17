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
include { minimap2 }    from './processes/minimap2.nf'
include { samtools }    from './processes/samtools.nf'
include { Binning }  from './processes/binning.nf'
include { checkm2 }     from './processes/checkm2.nf'
include { prokka }      from './processes/prokka.nf'
include { KOfamscan }   from './processes/kofamscan.nf'
include { quast }       from './processes/quast.nf'
include { dastool }  from './processes/dastool.nf'
include { kegg_pathway } from './processes/kegg_pathway.nf'
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

    // Step 3a: Mapping & Binning
    // binningOutput[0]  = sample_id
    // binningOutput[6]  = concoct_contig_bin.tsv
    // binningOutput[7]  = metabat_contig_bin.tsv
    // binningOutput[10] = metabinner_contig_bin.tsv
    minimap2Output = minimap2(fastpOutput[0], assembly)
    samtoolsOutput = samtools(minimap2Output)
    binningOutput  = Binning(samtoolsOutput, assembly)

    // Step 3b: Consensus Binning
    // dastoolOutput[0] = sample_id
    // dastoolOutput[1] = DAStool_out (full output dir)
    // dastoolOutput[2] = DASTool_bins (refined bins dir — fed downstream)
    dastoolOutput  = dastool(
        binningOutput[0],
        assembly,
        binningOutput[6],
        binningOutput[7],
        binningOutput[10]
    )

    // Step 3c: Quality Assessment on refined bins only
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

    // Step 5: Assembly QC on refined bins
    quastReport = quast(dastoolOutput[0], dastoolOutput[2])

    kegg_results = kegg_pathway(
    dastoolOutput[0],
    kofam_annotation[1]    // kofamscan_filtered.tsv
    )

    // Optional
    // gtdbtk_classification = gtdbtk(dastoolOutput[2])}
  }















































