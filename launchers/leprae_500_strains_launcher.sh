

nextflow main.nf -resume \
	--referenceSequence "../datasets/Leprae/reference_genomes/mycobrowser_leprae.fasta.gz" \
	--referenceAnnotation "../datasets/Leprae/reference_genomes/mycobrowser_leprae.gff.gz" \
	--results "../datasets/Leprae/500_strains" --vcf "../datasets/Leprae/outgroup" --outgroup "lepromatosis" \
	--species "other" --tree_build "fastest" --SKIP 'COVERAGE_CHECK'
