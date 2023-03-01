#!/bin/bash

nextflow main.nf -resume \
	--sra "/home/alm/Desktop/datasets/MTBC/DATATEST/SRA.txt" \
	--results "/home/alm/Desktop/datasets/MTBC/DATATEST/25_CLASSIC_STRAINS" \
	$@
