process gtdbtk {
    tag "GTDB-Tk classification on $sample_id"
    publishDir "${params.outdir}/${sample_id}/gtdbtk_classification", mode: 'copy'

    input:
    path(dastool_bins)
    output:
    file("gtdbtk_classify_out/dastool_bins")

    script:
    """
    # Running GTDB-Tk on the MetaBAT bins
    singularity exec /home/luiszuniga/work/containers/gtdbtk/gtdbtk_2.4.0.sif gtdbtk classify_wf --genome_dir dastool_bins --out_dir gtdbtk_classify_out --cpus 2 --extension fa --skip_ani_screen

    """
}
