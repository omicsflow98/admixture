process vcf2phylip {
        
    label 'vcf2phylip'
       
    publishDir "${params.outdir}/output/phylip_tree"
    container "${params.apptainer}/python.sif"

    input:
    path filtered_vcf
    path vcf2phylip_file

    output:
    path("sample.min4.phy"), emit: phylip_file

    script:

    """
    python ${vcf2phylip_file} \
    -i ${filtered_vcf} \
    --output-prefix sample
      
    """

}