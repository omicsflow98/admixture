process plot_results {
        
    label 'plot_results'
       
    publishDir "${params.outdir}/output/final_plots"
    container "${params.apptainer}/plot_admixture.sif"

    input:
    tuple path(admixture_out),path(CV_file)
    tuple path(eigenval), path(eigenvec)
    path newick_tree
    path plot_admixture
    val mink 
    val maxk
    
    output:
    path("*.png"), emit: plots

    script:

    """
    Rscript ${plot_admixture} \
    -Q ${admixture_out} \
    -V ./merged.filtered \
    -E ${CV_file} \
    -T ${newick_tree} \
    -S ${mink} \
    -B ${maxk}
      
    """

}
