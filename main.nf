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
include { Binning }  from './processes/binning.nf'
include { prokka }      from './processes/prokka.nf'
include { KOfamscan }   from './processes/kofamscan.nf'
include { quast }       from './processes/quast.nf'
// include { dastool }  from './processes/dastool.nf'
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
    binningOutput  = Binning(samtoolsOutput, assembly)

    // Step 4: Annotation
    protein_annotation = prokka(binningOutput[0], assembly, binningOutput[5], binningOutput[4], binningOutput[6])
    kofam_annotation   = KOfamscan(binningOutput[0], protein_annotation[1], protein_annotation[2])

    // Step 5: Quality Control & Reporting
    quastReport = quast(binningOutput[0], binningOutput[5], binningOutput[4])

    // Optional steps (commented out)
    // dastoolOutput = dastool(binningOutput[0], assembly, binningOutput[4], binningOutput[3])
    // gtdbtk_classification = gtdbtk(binningOutput[4], binningOutput[3])
}















































