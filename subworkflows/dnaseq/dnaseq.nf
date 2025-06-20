#!/usr/bin/env nextflow

include { alignment } from './alignment.nf'
include { variant_calling } from './variant_calling.nf'

workflow dnaseq {

	if (!params.start_vcf) {
		alignment()
		bam_files = alignment.out
	} else {
		bam_files = Channel.fromPath(params.data_csv, checkIfExists: true)
		| splitCsv(header: true, sep: '\t')
		| map { row -> tuple(row.SampName,
				file(row.File1)) 
			}
	}

	variant_calling(bam_files)
}