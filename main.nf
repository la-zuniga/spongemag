// main.nf
nextflow.enable.dsl=2

// Define project directory and output directory
projectDir = "/home/luiszuniga/self/MAG"
params.input_dir = "/home/luiszuniga/self/data"  // Default location
params.reads = "$params.input_dir/*_{R1,R2}.fastq.gz"
params.outdir = "$projectDir/out"

include { fastp } from './processes/fastp.nf'
include { fastqc } from './processes/FastQC.nf'
include { multiqc } from './processes/multiqc.nf'
include { megahit } from './processes/megahit.nf'
include { spades } from './processes/spades.nf'
include { minimap2 } from './processes/minimap2.nf'
include { samtools } from './processes/samtools.nf'
include { MetaBinner } from './processes/binning.nf'
include { prokka } from './processes/prokka.nf'
// include { dastool } from './processes/dastool.nf'
include { KOfamscan } from './processes/kofamscan.nf'
include { quast } from './processes/quast.nf'

// Channel for read pairs
Channel.fromFilePairs(params.reads).set { read_pairs_ch }

workflow {
    fastpOutput = fastp(read_pairs_ch)
    fastQCoutput = fastqc(fastpOutput[0])
    multiqc = multiqc(fastQCoutput[0])

    if (params.assembly == 'megahit') {
        assembly = megahit(fastpOutput[0])
    } else if (params.assembly == 'spades') {
        assembly = spades(fastpOutput)
    } else {
        error "Invalid assembler choice: ${params.assembly}"
    }
    minimap2Output = minimap2(fastpOutput[0], assembly)
    samtoolsOutput = samtools(minimap2Output)
    
    // Correctly combine the BAM file and contigs file as a tuple for MetaBinner
    metabinner = MetaBinner(samtoolsOutput, assembly)
    protein_annotation = prokka(metabinner[0], assembly, metabinner[5], metabinner[4])
    kofam_annotation = KOfamscan(metabinner[0], protein_annotation[1], protein_annotation[2])
    // dastool = dastool(metabinner[0], assembly, metabinner[4], metabinner[3])
    quast = quast(metabinner[0], metabinner[5], metabinner[4])
    // gtdbtk_classification = gtdbtk(metabinner[4], metabinner[3])
}
