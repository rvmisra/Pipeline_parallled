# Pipeline_paralled
USE AT OWN RISK.
Pipeline geared towards my server setup, which pulls together various 3rd party tools such as gatk, bwa and snpeff.  Each has optimised to make use of paralled processing, to make it as fast as posssible.  Its based on GATK recommendations.

NOTE: Its best to use a reference genome already in snpeff as adding a new snpeff reference genome is a pain to do.
So, check snpeff first, make sure the fasta reference name is the same as the snpeff genome i.e. >Chromosome.

NOTE: This currently setup for a single interlaced fastq, change accordingly for PE fastq data.
