
# Genotube v3.0
Version : Tokyo Drift

<https://github.com/adrienlemeur/genotube-td/>

Genotube is a genomic workflow for the study of Mycobacterium tuberculosis and affiliated species (only MTBC supported for now). It includes download, trimming, alignement, variant calling, lineage caracterrisation and antibioresistance predition for M. tuberculosis.

## Installation
Genotube requires :
- Nextflow (> 21.0)
- Singularity (or Docker) (any recent version should do the trick)

Both Nextflow and Singularity can be installed with Conda or Mamba. Docker must be installed with root rights.


## Documentation

[Genotube wiki](https://github.com/adrienlemeur/genotube/wiki) (coming later)


## Quick start

**Running Genotube** :

Command :
```
nextflow main.nf \\
  --sra data/SRA.txt \\
  --referenceSequence data/mtuberculosis/genomes/reference/H37Rv.fa.gz \\
  --referenceAnnotation data/mtuberculosis/genomes/reference/H37Rv.fa.gz
```
**MAIN ARGUMENTS**
* --sra : path to a line separated file with NCBI run ID to be download. Default : a small dataset of 25 run ID of strains from different lineage
* --referenceSequence : reference sequence for mapping
* --referenceAnnotation : GFF3 annotation file for the provided reference sequence
* --results : Folder where intermediary and output files are stored. Genotube will use these files to resume the analysis if interrupted.
> NOTE : SRA WILL NOT be redownloaded if either raw fastq, trimmed fastq, raw bam, bam, raw vcf or annotated vcf are available. To force re-download, use --FORCE_DOWNLOAD flag

**IMPORTING FILES**
* --fasta : path to a folder with fasta samples. Raw reads will be simulated from the fasta and input in Genotube.
* --single_fastq : path to a folder containing single fastq to import to the analysis
* --paired_fastq : path to a folder containing paired fastq to import to the analysis

**GENE ORIENTED**
* --target_region : bed file specifying specific regions to analyse
> NOTE : Make sure bed file CHROM column must match reference annotation.

**OTHER**
* --no_profiling : don't perform lineage classification and in silico drug resistance
* --no_quality_check : by default, Genotube discards samno_profilingple with mean cov < 80% and sequencing depth < 10. Disable this filter.

**PHYLOGENY**
* --tree : build a SNP tree with RAxML. **If --target_region was specified, only SNP overlapping with target regions will be used**
* --bigtree : build a SNP tree with fasttree. Faster for large datasets but less consistent.

**PHYLOGENETICAL PLACEMENT (UNDER DEVELOPPMENT) **
Do not change unless you want to use you own reference samples / taxonomical ranks
* --referenceSNP : Bed file with the position of reference SNP
* --referenceMSA : MSA file with full genome sequence of reference samples
> Note : the MSA width must match the size of the reference genome ! (by default : H37Rv strain)
* --referenceTree : reference tree for phylogenetical placement. Must match the MSA file.
* --referenceModel : RAxML evolution model of the reference tree
* --referenceTaxo : reference strain taxonomy
