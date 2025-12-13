process plot_results {
        
    label 'plot_results'
       
    publishDir "${params.outdir}/output/final_plots"

    input:
    tuple path(admixture_out),path(CV_file)
    tuple path(eigenval), path(eigenvec)
    path phylip_file
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
    -P ${phylip_file} \
    -T ${newick_tree} \
    -S ${mink} \
    -B ${maxk}
      
    """

}
