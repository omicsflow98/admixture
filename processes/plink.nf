process plink {
        
    label 'plink'
       
    publishDir "${params.outdir}/output/plink"
    container "${params.apptainer}/plink.sif"

    input:
    path(pruned_vcf)

    output:
    tuple path("*.bed"), path("*bim"), path("*fam"), emit: admixture_bed
    tuple path("*.eigenval"), path("*.eigenvec"), emit: eigenvalues
   
    script:

    """
    plink \
    --vcf ${pruned_vcf} \
    --make-bed \
    --allow-extra-chr \
    --pca \
    --out \$(basename ${pruned_vcf} .vcf.gz)
        
    """

}