process phylo_tree {
        
    label 'iqtree'
       
    publishDir "${params.outdir}/output/phylip_tree"
    container "${params.apptainer}/iqtree.sif"

    input:
    path phylip_file

    output:
    path("*.treefile"), emit: newick_tree

    script:

    """
    iqtree -s ${phylip_file} -nt AUTO
      
    """

}