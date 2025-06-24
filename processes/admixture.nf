process run_admixture {
        
    label 'admixture'
       
    publishDir "${params.outdir}/output/admixture"
    container "${params.apptainer}/admixture.sif"

    input:
    tuple path(bed_file), path(bim_file), path(fam_file)
    val min_k
    val max_k

    output:
    tuple path("*.P"), path("*.Q"), path("CV_values.log"), path("samples.txt"), emit: admixture_out
   
    script:

    """
    chromosomes=\$(cut -f1 ${bim_file} | uniq)
    i=0

    for chrom in \$chromosomes; 
    do
    ((++i))
    sed "s/\${chrom}/\${i}/g" ${bim_file} > test.out
    mv test.out ${bim_file}
    done

    for (( k=$min_k; k<=$max_k; k++ ));
    do
    admixture --cv ${bed_file} \$k | tee \$(basename ${bed_file} .bed).\$k.log
    done
       
    grep -h CV *.log | cut -f3,4 -d' ' | sed -E 's/.*K=([0-9]+)\\): ([0-9.eE+-]+)/\\1 \\2/' > CV_values.log

    cut -f1 -d' ' ${fam_file} > samples.txt

    """

}