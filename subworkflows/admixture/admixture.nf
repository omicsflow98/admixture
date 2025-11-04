#!/usr/bin/env nextflow

include { plink } from '../../processes/plink.nf'
include { run_admixture } from '../../processes/admixture.nf'
include { filter_vcf } from '../../processes/filter_vcf.nf'
include { plot_results } from '../../processes/plot_results.nf'
include { vcf2phylip } from '../../processes/vcf2phylip.nf'
include { phylo_tree } from '../../processes/iqtree.nf'

workflow admixture {
	take:
	merged_vcf

	main:

	plot_admix = file("${projectDir}/scripts/plot_results.R")

	filter_vcf(merged_vcf)

	plink(filter_vcf.out.filtered_vcf)

	run_admixture(plink.out.admixture_bed, params.kmin, params.kmax)

	vcf2phylip(filter_vcf.out.filtered_vcf)

	phylo_tree(vcf2phylip.out.phylip_file)

	plot_results(run_admixture.out.admixture_out, plink.out.eigenvalues, phylo_tree.out.newick_tree, plot_admix, params.mink, params.maxk)


	emit:
	plots = plot_results.out.plots

}