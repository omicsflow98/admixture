#!/usr/bin/env nextflow

include {runfastqc} from './processes/fastqc.nf'
include {trim_galore} from './processes/trim_galore.nf'
include {bwa} from './processes/bwa.nf'
include {markduplicates} from './processes/markdup.nf'
include {BQSR} from './processes/bqsr.nf'
include {gvcf} from './processes/VariantCalling.nf'
include {mergedgvcf} from './processes/VariantCalling.nf'

workflow {

    first_ch = Channel.fromPath(params.data_csv, checkIfExists: true)
		| splitCsv(header: true, sep: '\t')
		| map { row -> tuple(row.SampName,
				file(row.File1),
				file(row.File2),
				row.Adapter1,
				row.Adapter2,
				row.LibName,
				row.Barcode,
				row.Platform) }

	bwa_index = file( "${params.reference}/index/bwa/genome.fa.{,amb,ann,bwt,pac,sa}" )

	runfastqc(first_ch)

	trim_galore(first_ch)

	bwa(bwa_index, trim_galore.out.readgroup, trim_galore.out.trimmed_fastq)

	markduplicates(params.temp_dir, bwa.out.bam_files)

	BQSR(params.reference, params.temp_dir, markduplicates.out.mark_dup)

	//gvcf(params.reference, params.temp_dir, BQSR.out.bam_bqsr)

	gvcf(params.reference, params.temp_dir, markduplicates.out.mark_dup)

	gvcf.out.vcf_files
	| collect
	| set {collected_vcf}

	mergedgvcf(params.reference, params.temp_dir, collected_vcf)
}
