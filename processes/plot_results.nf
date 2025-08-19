process plot_results {
        
    label 'plot_results'
       
    publishDir "${params.outdir}/output/final_plots"
    container "${params.apptainer}/plot_admixture.sif"

    input:
    tuple path(admixture_out),path(CV_file)
    tuple path(eigenval), path(eigenvec)

    output:
    path("*.png"), emit: plots

    script:

    """
    plot_admixture.R -Q ${admixture_out} -V ./merged_filtered -E ${CV_file}
      
    """

}