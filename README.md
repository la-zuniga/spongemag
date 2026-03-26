# spongeMAG

A [Nextflow](https://www.nextflow.io/) pipeline for recovering metagenome-assembled genomes (MAGs) from shotgun sequencing data, with a focus on functional annotation of nitrogen cycle genes.

## Overview

spongeMAG takes paired-end shotgun metagenomic reads and carries them through quality control, assembly, binning, and functional annotation. It was inspired by the pipeline I developed during my master's thesis, which used metagenomics to profile nitrogen cycling potential in microbial symbionts associated with marine sponges.

The part I find most interesting is the annotation of genes using **KOfamScan**, which assigns KEGG orthologs (KOs) to proteins predicted from the recovered MAGs. These annotations can be used to infer the functional potential of a genome. The current version specifically targets KOs involved in nitrogen cycle reactions but it can be adapted to search for any KOs of interest.

## Pipeline steps

1. **Preprocessing** — Adapter trimming and quality filtering with fastp, followed by FastQC and MultiQC reporting.
2. **Assembly** — De novo metagenome assembly using either MEGAHIT or metaSPAdes (configurable).
3. **Binning** — Reads are mapped back to the assembly with minimap2. Three binners (CONCOCT, MetaBAT2, MetaBinner) run in parallel, and their results are reconciled into consensus bins using DAS Tool.
4. **Quality assessment** — CheckM2 evaluates bin completeness and contamination. QUAST reports assembly statistics on the refined bins.
5. **Annotation** — Prokka predicts genes on both the refined bins and the full assembly. KOfamScan then assigns KEGG orthologs to the predicted proteins, and results are mapped to KEGG pathways.
6. **Nitrogen cycle phylogenies** — Phylogenetic trees are built for nitrogen cycle genes recovered from both the bins and the full assembly.

GTDB-Tk taxonomic classification is implemented but currently commented out due to container size.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- [Singularity](https://sylabs.io/singularity/) — each module runs inside a Singularity container (containers are not currently hosted publicly)

## Usage

```bash
nextflow run main.nf --assembly megahit
# or
nextflow run main.nf --assembly spades
```

Input reads are expected as paired-end FASTQ files matching the pattern `*_{R1,R2}.fastq.gz`. Edit `nextflow.config` and the paths at the top of `main.nf` to point to your data and container directories.

## Project status

This is an active work in progress — still rough around the edges, but I'm iterating on it steadily. I'm making it public because I think it's fun, even in its unpolished state. Contributions, feedback, and suggestions are welcome.
