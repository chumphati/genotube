#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//delete fastq files when finished to be used
process deleteFastq {
	label 'script'

	input:
	tuple val(sample), file(f1), file(f2)
	val(results)
  
	when:
	!params.keep_fastq && !params.keep_raw_fastq

	script:  
	"""
	rm -rf ${results}/FQ_PE_RAW/${sample}*.fastq.gz ${results}/FQ_SE_RAW/${sample}*.fastq.gz
	"""
}

//delete fastq trimmed files when finished to be used
process deleteTrimmedFastq {
	label 'script'

	input:
	tuple val(sample), file(bam)
  	val(results)

	when:
	!params.keep_fastq && !params.keep_trimmed_fastq

	script:
 	"""
	rm -rf ${results}/FQ_PE_TRIMMED/${sample}*.fastq.gz ${results}/FQ_SE_TRIMMED/${sample}*.fastq.gz
  	"""
}

//delete raw bam if flag
process deleteBam {
	label 'script'

	input:
	tuple val(sample), file("${sample}.vcf.gz")
	val(results)
  
	when:
	!params.keep_bam && !params.keep_raw_bam

	script:
	"""
	rm -rf ${results}/BAM_RAW/${sample}*.bam
	"""
}

//delete processed bam if flag
process deleteProcessedBam {
	label 'script'

	input:
	tuple val(sample), file(vcf)
	val(results)
  
	when:
	!params.keep_bam && !params.keep_processed_bam

	script:
	"""
	rm -rf ${results}/BAM/${sample}*.bam
	"""
}

process deleteVCF {
	label 'script'

	input:
	tuple val(sample), file(vcf)
	val(results)
  
	when:
	!params.keep_vcf

	script:
	"""
	rm -rf ${results}/VCF_RAW/${sample}*.vcf.gz
	"""
}

process deleteCoverage {
	label 'script'

	input:
	tuple val(sample), file(vcf)
	val(results)

	when:
	!params.keep_reports && !params.keep_coverage

	script:
	"""
	rm -rf ${results}/BAM_COVERAGE/${sample}*_QC.txt
	"""
}

process deleteMetrics {
	label 'script'

	input:
	tuple val(sample), file(vcf)
	val(results)
  
	when:
	!params.keep_reports && !params.keep_metricsreports
  	
	script:
	"""
	rm -rf ${results}/ALL_REPORTS/${sample}*_dedup_metrics.txt
	"""
}

workflow cleaner {
	take: trim_info
	take: bam
	take: process_bam
	take: process_vcf
	take: strain_info

	main:
		results = file(params.results)

		deleteFastq(trim_info, results)
		deleteTrimmedFastq(bam, results)

		deleteBam(process_bam, results)
		deleteProcessedBam(process_vcf, results)

		deleteVCF(process_vcf, results)

		deleteCoverage(process_vcf, results)
		deleteMetrics(process_vcf, results)
}
