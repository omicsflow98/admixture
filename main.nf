#!/usr/bin/env nextflow

include { dnaseq } from './subworkflows/dnaseq/dnaseq.nf'

workflow {
	
    dnaseq()
}
