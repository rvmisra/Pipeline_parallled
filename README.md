# Pipeline_parallled
Bash script

Use at own risk.

Bsed on GATK recommendations.

Cange file inputs in the script and run pipeline on HPC servers

NOTE:

1.For snpEFF to work effectively, need to select the appropriate snpEFF database and set that as the reference for in the pipeline.
If new/unusual species, comment out snpEFF and run it seperately.  Creating a snpEFF database is non-trivial

2.Currently setup for a single interlaced fastq file, change pipeline where needed to process individual fastqs from a PE run
