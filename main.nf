nextflow.enable.dsl=2

// =============================
// Project Directories & Params
// =============================
projectDir = "/home/luiszuniga/self/MAG"
params.input_dir = "/home/luiszuniga/self/data"
params.reads = "$params.input_dir/*_{R1,R2}.fastq.gz"
params.outdir = "$projectDir/out"

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
include { MetaBinner }  from './processes/binning.nf'
include { prokka }      from './processes/prokka.nf'
include { KOfamscan }   from './processes/kofamscan.nf'
include { quast }       from './processes/quast.nf'
// include { dastool }  from './processes/dastool.nf'
// include { gtdbtk }   from './processes/gtdbtk.nf'

// =============================
// Input Channels
// =============================
Channel.fromFilePairs(params.reads).set { read_pairs_ch }

// =============================
// Workflow Definition
// =============================
workflow {

    // Step 1: Preprocessing
    fastpOutput    = fastp(read_pairs_ch)
    fastQCoutput   = fastqc(fastpOutput[0])
    multiqcReport  = multiqc(fastQCoutput[0])

    // Step 2: Assembly
    if (params.assembly == 'megahit') {
        assembly = megahit(fastpOutput[0])
    } else if (params.assembly == 'spades') {
        assembly = spades(fastpOutput)
    } else {
        error "Invalid assembler choice: ${params.assembly}"
    } 

    // Step 3: Mapping & Binning
    minimap2Output    = minimap2(fastpOutput[0], assembly)
    samtoolsOutput    = samtools(minimap2Output)
    metabinnerOutput  = MetaBinner(samtoolsOutput, assembly)

    // Step 4: Annotation
    protein_annotation = prokka(metabinnerOutput[0], assembly, metabinnerOutput[5], metabinnerOutput[4], metabinnerOutput[6])
    kofam_annotation   = KOfamscan(metabinnerOutput[0], protein_annotation[1], protein_annotation[2])

    // Step 5: Quality Control & Reporting
    quastReport = quast(metabinnerOutput[0], metabinnerOutput[5], metabinnerOutput[4])

    // Optional steps (commented out)
    // dastoolOutput = dastool(metabinnerOutput[0], assembly, metabinnerOutput[4], metabinnerOutput[3])
    // gtdbtk_classification = gtdbtk(metabinnerOutput[4], metabinnerOutput[3])
}















































