#!/usr/bin/env nextflow

include { plink } from '../../processes/plink.nf'
include { run_admixture } from '../../processes/admixture.nf'
include { prune_snp } from '../../processes/prune_snp.nf'
include { plot_results } from '../../processes/plot_results.nf'

workflow admixture {
	take:
	merged_vcf

	main:

	prune_snp(merged_vcf)

	plink(prune_snp.out.pruned_vcf)

	run_admixture(plink.out.admixture_bed, params.kmin, params.kmax)

	plot_results(run_admixture.out.admixture_out, plink.out.eigenvalues)


	emit:
	plots = plot_results.out.plots

}