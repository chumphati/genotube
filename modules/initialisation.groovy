#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process printHelp {
	input:
	val errorMessage

	output:
	val false

	when:
	params.help

	exec:
	"cat $errorMessage".execute().text.readLines().each{println it}
	// a bit ugly but would not work otherwise
}

process createFolders {
	label 'script'

	input:
	val results

	output:
	val false

	when:
	params.dry && !params.help

	exec:
	//Fasta folders for assemble genome input
	file("$results/FASTA").mkdirs()

	//Fastq folders
	file("$results/FASTQ/RAW").mkdirs()
	file("$results/FASTQ/FILTERED").mkdirs()

	file("$results/ALL_REPORTS/FASTQ/ABUNDANCE").mkdirs()
	file("$results/ALL_REPORTS/FASTQ/COMPET_MAPPING").mkdirs()
	file("$results/ALL_REPORTS/FASTQ/KMERFILTER").mkdirs()
	file("$results/ALL_REPORTS/FASTQ/TAXOCLASS").mkdirs()
	file("$results/ALL_REPORTS/FASTQ/TRIMMING").mkdirs()

	//BAM folders
	file("$results/BAM/RAW").mkdirs()
	file("$results/BAM/FILTERED").mkdirs()
	file("$results/ALL_REPORTS/BAM/RAW").mkdirs()
	file("$results/ALL_REPORTS/BAM/QUAL").mkdirs()
	file("$results/ALL_REPORTS/BAM/DEDUP").mkdirs()

	//Variant results
	file("$results/VCF_RAW").mkdirs()
	file("$results/FILTERED_VCF").mkdirs()
	file("$results/OUT_VCF").mkdirs()

	//Quality control folders
	file("$results/QC_RESUME_LINEAGE_SNP").mkdirs()
	file("$results/QC_FULL_LINEAGE_SNP").mkdirs()

	//Taxonomical classification
	file("$results/OUT_TAXONOMY").mkdirs()

	//Phylogenetic placement & taxonomical classification
	file("$results/OUT_PLACEMENT").mkdirs()
	file("$results/OUT_TAXONOMY").mkdirs()

	//M. tuberculosis specific outputs
	file("$results/OUT_LINEAGE").mkdirs()
	file("$results/ANTIBIO_PROFILE").mkdirs()

	file("$results/REPORTS").mkdirs()

	println "\n\nThe result folder environnement has been created !\n"
}

workflow initialisation {
	main:
	results = file(params.results)

	errorMessage = file("./data/help.txt")
	printHelp(errorMessage)

	createFolders(results)
}
