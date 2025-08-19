process run_admixture {
        
    label 'admixture'
       
    publishDir "${params.outdir}/output/admixture"
    container "${params.apptainer}/admixture.sif"

    input:
    tuple path(bed_file), path(bim_file), path(fam_file)
    val min_k
    val max_k

    output:
    tuple path("admixture.Q"), path("CV_values.log"), emit: admixture_out
   
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
    header=""
    for ((i=1; i<=k; i++));
    do
        header+="\$k:\$i "
    done
    header=\${header::-1}
    admixture --cv ${bed_file} \$k | tee \$(basename ${bed_file} .bed).\$k.log
    sed -i "1i\${header}" merged_filtered.\${k}.Q
    done
       
    grep -h CV *.log | cut -f3,4 -d' ' | sed -E 's/.*K=([0-9]+)\\): ([0-9.eE+-]+)/\\1 \\2/' > CV_values.log

    awk '{if (\$1 == \$2) print \$1; else print \$1"_"\$2}' ${fam_file} > samples.txt
    sed -i "1iSamples" samples.txt

    paste -d ' ' samples.txt *.Q > admixture.Q

    """

}