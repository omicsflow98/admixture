#!/usr/bin/env nextflow

include { plink } from '../../processes/plink.nf'
include { run_admixture } from '../../processes/admixture.nf'
include { prune_snp } from '../../processes/prune_snp.nf'

workflow admixture {
	take:
	merged_vcf

	main:

	prune_snp(merged_vcf)

	plink(prune_snp.out.pruned_vcf)


	emit:
	vcf_file = plink.out.admixture_bed

}