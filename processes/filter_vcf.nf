process filter_vcf {
        
    label 'filter_vcf'
       
    publishDir "${params.outdir}/output/gvcf_merged", mode : 'copy'
    container "${params.apptainer}/gatk.sif"

    input:
    path merged_vcf

    output:
    path("merged.filtered.vcf.gz"), emit: filtered_vcf

    script:

    """
    bcftools +fill-tags ${merged_vcf} \
    -Oz -o temp.vcf.gz \
    -- -t MAF

    bcftools index temp.vcf.gz

    bcftools view temp.vcf.gz \
    -v snps \
    -m2 -M2 \
    -i 'QUAL>30 && F_MISSING=0 && MAF>=0.1 && INFO/DP>10 && INFO/DP<6000' \
    -Oz -o merged.filtered.vcf.gz
      
    """

}