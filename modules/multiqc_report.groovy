#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process multiqc {
	label 'multiqc'

	input:
	val(results)
	file(multiqc_config)
	val(flag)

	script:
	"""
	{
		multiqc $results/ALL_REPORTS -c $multiqc_config
		mv multiqc_report.html $results/report.html
	} || {
		:
	} #find something to do
	"""
}

workflow multiqc_report {
	take: wait_end_process
	take: filterVariants

	main:
	results = file(params.results)
	multiqc_config_file = file(params.multiqc_config)
	multiqc(results, multiqc_config_file, wait_end_process.collect().ifEmpty(true))    
}
