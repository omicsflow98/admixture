process run_admixture {
        
    label 'plink'
       
    publishDir "${params.outdir}/output/plink"
    container "${params.apptainer}/plink.sif"

    input:
    path(pruned_vcf)

    output:
    tuple path("*.bed"), path("*bim"), path("*fam"), path("*.log"), path("*.nosex"), emit: admixture_bed
    tuple path("*.eigenval"), path("*.eigenvec"), emit: eigenvalues
   
    script:

    """
        
    """

}