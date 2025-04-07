### Overview

April 7, 2025
This is a pipeline I've been working on to recover MAGs from shot-gun sequencing data.
It has been a personal project to help me learn nextflow and how to build and run containerized tools for reproducible results.
This project is still rough around the edges, but I am actively working to improve it. I also think it's fun, so I'm making it public despite its unpolished state. Each module relies on one or more singulartiy containers that I've built. Currently, I don't have these containers hosted anywhere. 

It is inspired by the pipeline I worked on during my master's thesis, and so the part of this that I think is most interesting is the annotation of genes using KOFamScan. This assigns proteins translated from the nuclecid acid sequences from the recovered MAGs a kegg ortholog (KO). These types of annotations can be used to try and infer the functional potential of a genome. 

In this version that I have up right now, the KOFamScan module specifically looks for KO's that are from nitrogen cycle reactions (i.e. the proteins that mediate chemical reactions that move nitrogen between biotic and abiotic compartments). My thesis was a metagenomic profiling of the nitrogen cycling potential in marine sponges, so that's the reason for the specificity. One could search for any KO's the are interested in, however. 

TO DO:
Proabably need to split checkm2 away from the binning module. Currently, checkm2 runs at the end of the binning, but I think splitting it would allow for a more organized output.

