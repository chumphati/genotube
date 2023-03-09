#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process paired_alignement {
	tag "$sample"
	label 'bwa'

	input:
	tuple val(sample), file(f1), file(f2)
	tuple val(indexName), file(indexFiles)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.bam")

	when:
	((! mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	! mf.checkFile("$results/BAM", sample, "bam") && \
	! mf.checkFile("$results/VCF_FILTERED", sample, "vcf.gz") && \
	! mf.checkFile("$results/VCF_RAW", sample, "vcf.gz")) && !params.dry && !params.help) || mf.checkFORCE('MAP', params.FORCE)

	script:
	//if flag param_region is set, convert to file
	if(params.target_region){region = file(params.target_region)}
	if(!params.target_region)
		"""
		bwa-mem2 mem -t ${task.cpus} $indexName $f1 $f2 -R "@RG\tID:$sample\tLB:LIB\tPL:Illumina\tPU:UNIT\tSM:$sample" 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samclip --ref ${fasta} | \
			samtools view -@ ${task.cpus} -bSh -q 30 -F0x4 --no-PG 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samtools sort -@ ${task.cpus} --no-PG -o ${sample}.bam

		rm -rf $results/BAM/RAW/${sample}.bam
		mv ${sample}.bam $results/BAM/RAW/${sample}.bam
		ln -s $results/BAM/RAW/${sample}.bam ${sample}.bam
		"""
	else
	//filter out reads non overlapping with target regions
		"""
		bwa-mem2 mem -t ${task.cpus} $indexName $f1 $f2 -R "@RG\tID:$sample\tLB:LIB\tPL:Illumina\tPU:UNIT\tSM:$sample" 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samclip --ref ${fasta} | \
			samtools view -@ ${task.cpus} -bSh -q 30 -F0x4 -L $region --no-PG 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samtools sort -@ ${task.cpus} --no-PG -o ${sample}.bam

		rm -rf $results/BAM/RAW/${sample}.bam
		mv ${sample}.bam $results/BAM/RAW/${sample}.bam
		ln -s $results/BAM/RAW/${sample}.bam ${sample}.bam
		"""
}

process single_alignement {
	tag "$sample"
	label 'bwa'

	input:
	tuple val(sample), file(f1)
	tuple val(indexName), file(indexFiles)
	tuple file(fasta), file(fai), file(dict)
	val(results)

	output:
	tuple val(sample), file("${sample}.bam")

	when:
	((! mf.checkFile("$results/BAM/RAW", sample, "bam") && \
	! mf.checkFile("$results/BAM", sample, "bam") && \
	! mf.checkFile("$results/VCF_FILTERED", sample, "vcf.gz") && \
	! mf.checkFile("$results/VCF_RAW", sample, "vcf.gz")) && !params.help && !params.dry) || mf.checkFORCE('MAP', params.FORCE)

	script:
	//if flag param_region is set, convert to file
	if(params.target_region){region = file(params.target_region)}
	if(!params.target_region)
		"""
		bwa-mem2 mem -t ${task.cpus} $indexName $f1 -R "@RG\tID:$sample\tLB:LIB\tPL:Illumina\tPU:UNIT\tSM:$sample" 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samclip --ref ${fasta} | \
			samtools view -@ ${task.cpus} -bSh -q 30 -F0x4 --no-PG 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samtools sort -@ ${task.cpus} --no-PG -o ${sample}.bam

		rm -rf $results/${sample}.bam
		mv ${sample}.bam $results/BAM/RAW/${sample}.bam
		ln -s $results/BAM/RAW/${sample}.bam ${sample}.bam
		"""
	else
	//filter out reads non overlapping with target regions
		"""
		bwa-mem2 mem -t ${task.cpus} $indexName $f1 -R "@RG\tID:$sample\tLB:LIB\tPL:Illumina\tPU:UNIT\tSM:$sample" 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samclip --ref ${fasta} | \
			samtools view -@ ${task.cpus} -bSh -q 30 -F0x4 --no-PG -L $region 2>> $results/ALL_REPORTS/BAM/RAW/${sample}.align.log | \
			samtools sort -@ ${task.cpus} --no-PG -o ${sample}.bam

		rm -rf $results/${sample}.bam
		mv ${sample}.bam $results/BAM/RAW/${sample}.bam
		ln -s $results/BAM/RAW/${sample}.bam ${sample}.bam
		"""
}

mf = new myFunctions()

workflow align {
	take: single_trimmed_fastq
	take: paired_trimmed_fastq
	take: bwa_index
	take: samtools_index

	main:
		results = file(params.results)
		file("$results/BAM/RAW").mkdirs()
		file("$results/ALL_REPORTS/BAM/RAW").mkdirs()

		paired_alignement(paired_trimmed_fastq, bwa_index, samtools_index, results)
		single_alignement(single_trimmed_fastq, bwa_index, samtools_index, results)

		//merge newly mapped bam channel with bam alread stored in the BAM_RAW folder into a new flux
		if( mf.checkFORCE('MAP', params.FORCE) ){
			old_bam = Channel.empty()
		} else if( params.bam == false ) {
			old_bam = Channel.fromPath(results+"/BAM/RAW/*.bam").map{it -> [it.simpleName, it]}
		} else {
			old_bam = Channel.fromPath([results+"/BAM/RAW/*.bam", params.bam+"/*.bam"]).map{it -> [it.simpleName, it]}
			old_bam.map{it -> it[1]}.flatten()
				.subscribe{ it -> mf.createSymLink(it.toString(), results.toString()+"/BAM/RAW") }
		}

		all_mapping = old_bam.mix(single_alignement.out.mix(paired_alignement.out))
	emit:
		all_mapping
		garbage = paired_alignement.out
}
