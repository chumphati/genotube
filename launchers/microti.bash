#!/usr/bin/env bash

nextflow main.nf -resume \
	--sra "/home/alm/Desktop/extending_the_barcodes/microti/sra.txt" \
	--results "/home/alm/Desktop/extending_the_barcodes/microti/genotube_results" \
	--tree_build "classic" $@
