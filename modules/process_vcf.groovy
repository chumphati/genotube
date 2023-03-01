#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process variantAnnotation {
	tag "$sample"
	label 'snpeff'

	input:
	tuple val(sample), file(vcf), file(csi)
	tuple file(fasta), file(fai), file(dict)
	file(snpeffConfig_file)
	val(start)
	val(resultsAnn)
	val(results)

	output:
	tuple val(sample), file("${sample}.vcf")
	file("${sample}_SNP_summary.txt")

	when:
	!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") || mf.checkFORCE('ANN', params.FORCE)

	script:
	"""
	snpEff ann -c $snpeffConfig_file ${fasta.simpleName} $vcf \
		-no-upstream -no-downstream -no-utr -no SPLICE_SITE_REGION \
		-noStats -csvStats ${sample}_SNP_summary.txt > ${sample}.vcf

	cp ${sample}_SNP_summary.txt $results/ALL_REPORTS/VCF/SNP_SUM/${sample}_SNP_summary.txt
	"""
}

process compressVCF {
	tag "$sample"
	label 'samtools'

	input:
	tuple val(sample), file(vcf)
	file("SNP_summary.txt")
	val(results)

	output:
	tuple val(sample), file("${vcf.simpleName}.vcf.gz")

	script:
	"""
	rm -rf $results/VCF/RAW/${sample}.vcf.gz
	bgzip $vcf
	cp ${vcf}.gz $results/VCF/RAW/${sample}.vcf.gz
	"""
}

process filterVariants {
	tag "$sample"
	label 'bcftools'

	input:
	tuple val(sample), file(vcf)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.filtered.vcf.gz"), emit: all_variants

	when:
	!mf.checkFile("$results/VCF/FILTERED", sample, "vcf.gz") || mf.checkFORCE('ANN', params.FORCE)

	script:
	"""
	bcftools filter $vcf \
		-i "(AF[0]>0.95 || AF[1] > 0.95) && (FORMAT/AD[:0]>=5 || FORMAT/AD[:1]>=5) && QUAL>30" \
		--SnpGap 2 --IndelGap 9 --soft-filter FAIL --mode x \
		-Oz -o ${sample}.filtered.vcf.gz

	rm -rf $results/VCF/FILTERED/${sample}.filtered.vcf.gz
	mv ${sample}.filtered.vcf.gz $results/VCF/FILTERED/${sample}.filtered.vcf.gz
	ln -s $results/VCF/FILTERED/${sample}.filtered.vcf.gz ${sample}.filtered.vcf.gz
	"""
}

mf = new myFunctions()

workflow process_vcf {
	take: raw_vcf
	take: coverage_info
	take: index
	take: snpeff_emit_signal

	main:
		snpeffConfig = file("data/snpeff/snpEff.config")
		results = file(params.results)
		file("$results/ALL_REPORTS/VCF/SNP_SUM").mkdirs()
		file("$results/VCF/FILTERED").mkdirs()

		variantAnnotation(raw_vcf, index, snpeffConfig, snpeff_emit_signal, results, results)
		compressVCF(variantAnnotation.out, results)

		filterVariants(compressVCF.out, index, results)

		all_ann_vcf = compressVCF.out
		end_signal = all_ann_vcf.ifEmpty(true)
	emit:
		all_ann_vcf
		all_raw_vcf = raw_vcf
		end_signal
		garbage = filterVariants.out
}
