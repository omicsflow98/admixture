process prune_snp {
        
    label 'prune_snp'
       
    publishDir "${params.outdir}/output/prune_snp"
    container "${params.apptainer}/prune_snp.sif"

    input:
    path merged_vcf

    output:
    path("*.vcf.gz"), emit: pruned_vcf

    script:

    """
    prune_snp.R -f ${merged_vcf}
      
    """

}