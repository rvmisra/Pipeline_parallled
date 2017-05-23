#!/bin/sh
#pipe1
#Define variables
#Reference_sequence_fasta/fna - using full path.
REFpath="/srv/data0/raju/PB/pips/Ref/Peptoclostridium_difficile_r20291.fna"
#Fastq sequence if sinle end use this only
FQ1="/srv/data0/raju/PB/cons_gens/fa2fq/497_Q2_BB_PE.fq"
#Fastq sequence, for PE.  Change pipeline for PE.
#FQ2=
#Output path for all outputs from each step of the pipeline exc. reference sequence
OUTPATH="/srv/data0/raju/PB/pips/"
#Sample prefix name for all file outputs
SamplePrefix="497_Q2_BB"

#trimiing
#Trimming not tested yet as part of this pipeline, as it was trialled using insilico fastq made from PB fastas.
#bbduk.sh -Xmx10g in1=read1.fq in2=read2.fq out1=clean1.fq out2=clean2.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 hdist=2

#1)See snpEFF reference sequence, use that name for the reference sequence file name.  As a first try, rename the fasta header to 'Chromosome' only.
#Run through and see what sequence name its looking for and change accordingly if it fails.
#2) Fasta reference has to be be '.fna'
#Reference seq formating

#Picard Create .dict
java -jar /home/raju.misra/programs/picard/picard-tools-2.6.0/picard.jar CreateSequenceDictionary R=$REFpath O=$REFpath.dict
#/usr/local/src/samtools-1.4/Samtools index
samtools faidx $REFpath
#BWA index
bwa index $REFpath

#BWA mapping and readgrouping
bwa mem -t 10 -M -R '@RG\tID:'$SamplePrefix'\tLB:'$SamplePrefix'\tPL:ILLUMINA\tPM:HISEQ\tSM:'$SamplePrefix'' $REFpath $FQ1 > $OUTPATH$SamplePrefix.sam

#Picard SAM to sorted BAM
java -jar /home/raju.misra/programs/picard/picard-tools-2.6.0/picard.jar SortSam INPUT=$OUTPATH$SamplePrefix.sam OUTPUT=$OUTPATH$SamplePrefix.sorted.bam SORT_ORDER=coordinate

#QC STATS
java -jar /home/raju.misra/programs/picard/picard-tools-2.6.0/picard.jar CollectAlignmentSummaryMetrics R=$REFpath I=$OUTPATH$SamplePrefix.sorted.bam O=$OUTPATH$SamplePrefix\_alignment_metrics.txt
java -jar /home/raju.misra/programs/picard/picard-tools-2.6.0/picard.jar CollectInsertSizeMetrics INPUT=$OUTPATH$SamplePrefix.sorted.bam OUTPUT=$OUTPATH$SamplePrefix\_insert_metrics.txt HISTOGRAM_FILE=$OUTPATH$SamplePrefix\_insertsize_histogram.pdf
/usr/local/src/samtools-1.4/samtools depth -a $OUTPATH$SamplePrefix.sorted.bam > $OUTPATH$SamplePrefix\_depthout.txt


#Mark duplicates
java -jar /home/raju.misra/programs/picard/picard-tools-2.6.0/picard.jar MarkDuplicates INPUT=$OUTPATH$SamplePrefix.sorted.bam OUTPUT=$OUTPATH$SamplePrefix\_sorted_dededup_reads.bam METRICS_FILE=$OUTPATH$SamplePrefix\_metrics.txt

#Build BAM index
java -jar /home/raju.misra/programs/picard/picard-tools-2.6.0/picard.jar BuildBamIndex INPUT=$OUTPATH$SamplePrefix\_sorted_dededup_reads.bam

#GATK create realignment targets
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 10 -R $REFpath -I $OUTPATH$SamplePrefix\_sorted_dededup_reads.bam -o $OUTPATH$SamplePrefix\_realignment_targets.list


#GATK realign indels
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T IndelRealigner -R $REFpath -I $OUTPATH$SamplePrefix\_sorted_dededup_reads.bam -targetIntervals $OUTPATH$SamplePrefix\_realignment_targets.list -o $OUTPATH$SamplePrefix\_realigned_reads.bam


#GATK call variants
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 10 -R $REFpath -I $OUTPATH$SamplePrefix\_realigned_reads.bam -o $OUTPATH$SamplePrefix\_RAW_variants.vcf


#GATK extract SNPS and InDels
#SNPS
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T SelectVariants -nt 10 -R $REFpath -V $OUTPATH$SamplePrefix\_RAW_variants.vcf -selectType SNP -o $OUTPATH$SamplePrefix\_raw_snps.vcf

#InDels
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T SelectVariants -nt 10 -R $REFpath -V $OUTPATH$SamplePrefix\_RAW_variants.vcf -selectType INDEL -o $OUTPATH$SamplePrefix\_raw_indels.vcf


#GATK Filter SNPs
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T VariantFiltration -R $REFpath -V $OUTPATH$SamplePrefix\_raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $OUTPATH$SamplePrefix\_filtered_snps.vcf


#GATK Filter InDels
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T VariantFiltration -R $REFpath -V $OUTPATH$SamplePrefix\_raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $OUTPATH$SamplePrefix\_filtered_indels.vcf


#GATK Base Quality Score Recalibration (BQSR)
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 10 -R $REFpath -I $OUTPATH$SamplePrefix\_realigned_reads.bam -knownSites $OUTPATH$SamplePrefix\_filtered_snps.vcf -knownSites $OUTPATH$SamplePrefix\_filtered_indels.vcf -o $OUTPATH$SamplePrefix\_recal_data.table


#GATK BQSR cont. 2
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T BaseRecalibrator -nct 10 -R $REFpath -I $OUTPATH$SamplePrefix\_realigned_reads.bam -knownSites $OUTPATH$SamplePrefix\_filtered_snps.vcf -knownSites $OUTPATH$SamplePrefix\_filtered_indels.vcf -BQSR $OUTPATH$SamplePrefix\_recal_data.table -o $OUTPATH$SamplePrefix\_post_recal_data.table


#GATK Analyse covariates - recalibration_plot
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $REFpath -before $OUTPATH$SamplePrefix\_recal_data.table -after $OUTPATH$SamplePrefix\_post_recal_data.table -csv $OUTPATH$SamplePrefix\_tmp_BQSR.csv -plots $OUTPATH$SamplePrefix\_recalibration_plots.pdf


#GATK Apply BQSR
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T PrintReads -nct 10 -R $REFpath -I $OUTPATH$SamplePrefix\_realigned_reads.bam -BQSR $OUTPATH$SamplePrefix\_recal_data.table -o $OUTPATH$SamplePrefix\_recal_reads.bam

#GATK call variants - based on BQSR tidies bam
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 10 -R $REFpath -I $OUTPATH$SamplePrefix\_recal_reads.bam -o $OUTPATH$SamplePrefix\_raw_variants_recal.vcf

#SNPS post BQSR
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T SelectVariants -nt 10 -R $REFpath -V $OUTPATH$SamplePrefix\_raw_variants_recal.vcf -selectType SNP -o $OUTPATH$SamplePrefix\_recal_raw_snps.vcf
#InDels
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T SelectVariants -nt 10 -R $REFpath -V $OUTPATH$SamplePrefix\_raw_variants_recal.vcf -selectType INDEL -o $OUTPATH$SamplePrefix\_recal_raw_indels.vcf

#GATK Filter SNPs
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T VariantFiltration -R $REFpath -V $OUTPATH$SamplePrefix\_recal_raw_snps.vcf  --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o $OUTPATH$SamplePrefix\_recal_filtered_snps_FINAL.vcf
#GATK Filter InDels
java -jar /home/raju.misra/programs/gatk3_7/GenomeAnalysisTK.jar -T VariantFiltration -R $REFpath -V $OUTPATH$SamplePrefix\_recal_raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o $OUTPATH$SamplePrefix\_filtered_indels_recal_FINAL.vcf
#COM1

#SNPeff annotate SNPs and predict effects
#See snpEFF reference sequence, use that name for the reference sequence file name.  As a first try, rename the fasta header to 'Chromosome' only.
#Run through and see what sequence name its looking for and change accordingly if it fails.
java -jar /home/raju.misra/programs/snpEFF/snpEff/snpEff.jar -no-downstream -no-upstream -no-utr -no-intergenic -v -o gatk Peptoclostridium_difficile_r20291 $OUTPATH$SamplePrefix\_recal_filtered_snps_FINAL.vcf -stats $OUTPATH$SamplePrefix\_recal_filtered_snps_FINAL_SNPEFF.html > $OUTPATH$SamplePrefix\_recal_filtered_snps_FINAL_SNPEFF.ann.vcf

