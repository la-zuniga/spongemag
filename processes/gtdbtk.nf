process gtdbtk {
    tag "GTDB-Tk classification on $sample_id"
    publishDir "${params.outdir}/${sample_id}/gtdbtk_classification", mode: 'copy'

    input:
    file("metabat_out")
    file("concoct_out")

    output:
    file("gtdbtk_classify_out/metabat_bins")
    file("gtdbtk_classify_out/concoct_bins")

    script:
    """
    # Running GTDB-Tk on the MetaBAT bins
    singularity exec /home/luiszuniga/work/containers/gtdbtk/gtdbtk_2.4.0.sif gtdbtk classify_wf --genome_dir metabat_out/ --out_dir gtdbtk_classify_out/metabat_bins --cpus 2 --extension fa --skip_ani_screen

    singularity exec /home/luiszuniga/work/containers/gtdbtk/gtdbtk_2.4.0.sif gtdbtk classify_wf --genome_dir concoct_out/fasta_bins/ --out_dir gtdbtk_classify_out/concoct_bins --cpus 2 --extension fa --skip_ani_screen
    """
}
