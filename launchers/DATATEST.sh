#!/bin/bash

nextflow main.nf -resume \
	--referenceSequence "../datasets/Leprae/reference_genomes/mycobrowser_leprae.fasta.gz" \
	--referenceAnnotation "../datasets/Leprae/reference_genomes/mycobrowser_leprae.gff.gz" \
	--results "../datasets/Leprae/DATATEST" \
	--species "other" $@
