process plot_results {
        
    label 'plot_results'
       
    publishDir "${params.outdir}/output/final_plots"
    container "${params.apptainer}/plot_admixture.sif"

    input:
    tuple path(P_files), path(Q_files), path(CV_file), path(samples)
    tuple path(eigenval), path(eigenvec)
    val kmin
    val kmax

    output:
    path("*.png"), emit: plots

    script:

    """
    plot_admixture.R -S ${samples} -A ./merged_filtered -E ${CV_file} -L ${kmin} -H ${kmax}
      
    """

}