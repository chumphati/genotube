#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process removeFastq {
	label 'script'

	input:
	val outalign
	val(results)

	when:
	!params.help && mf.checkFORCE('FASTQ', params.REMOVE)

	script:
	"""
	rm ${results}/FASTQ/RAW/* ${results}/FASTQ/TRIMMED/*
	"""
}

process removeBam {
	label 'script'

	input:
	val outvcf
	val(results)

	when:
	!params.help && mf.checkFORCE('BAM', params.REMOVE)

	script:
	"""
	rm ${results}/BAM/RAW/*
	if test -f '$results/BAM/FILTERED/*.bam'; then
		rm ${results}/BAM/FILTERED/*
	fi
	"""
}

process removeVCF {
	label 'script'

	input:
	val outallprocessvcf
	val(results)

	when:
	!params.help && mf.checkFORCE('VCF', params.REMOVE)

	script:
	"""
	rm ${results}/VCF/RAW/* ${results}/VCF/FILTERED/*
	"""
}

mf = new myFunctions()

workflow cleaner {
	take: all_mapping
	take: all_ann_vcf
	take: strain_info

	main:

	results = file(params.results)

	removeFastq(all_mapping, results)

	removeBam(all_ann_vcf, results)

	removeVCF(strain_info, results)

}