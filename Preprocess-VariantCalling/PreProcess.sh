#! /bin/bash

GENOME='../genome/canFam3.fa.gz'
VCF='../vcf1/Filtred_Published1.vcf.gz'

SAMPLE=$1

prefetch --max-size 400000000000 -p $SAMPLE 
fasterq-dump -e 20 -p './'$SAMPLE 
bgzip -@ 20 $SAMPLE'_1.fastq' 
bgzip -@ 20 $SAMPLE'_2.fastq'

bwa mem -t 20 $GENOME $SAMPLE'_1.fastq.gz' $SAMPLE'_2.fastq.gz' > Align/$SAMPLE'.aln.sam'
rm $SAMPLE'_1.fastq.gz' && rm $SAMPLE'_2.fastq.gz'

gatk SamFormatConverter -I Align/$SAMPLE'.aln.sam' -O Convert/$SAMPLE'.aln.bam' 
rm $SAMPLE'.aln.sam'

gatk SortSamSpark -I Convert/$SAMPLE'.aln.bam' -O Sorted/$SAMPLE'.sort.bam' -SO coordinate
rm $SAMPLE'.aln.bam'

gatk AddOrReplaceReadGroups -I Sorted/$SAMPLE'.sort.bam' -O Groups/$SAMPLE'.group.bam' -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM $SAMPLE
rm $SAMPLE'.sort.bam'

gatk MarkDuplicatesSpark -I Groups/$SAMPLE'.group.bam' -O MDups/$SAMPLE'.mark.bam' -M MDups/$SAMPLE'.fasta'
rm $SAMPLE'.group.bam'

gatk BaseRecalibratorSpark -I MDups/$SAMPLE'.mark.bam' -O BQRS/$SAMPLE'.table' -R $GENOME --known-sites $VCF
gatk ApplyBQSRSpark -I MDups/$SAMPLE'.mark.bam' -O BQRS/$SAMPLE'.ready.bam' -R $GENOME --bqsr-recal-file BQRS/$SAMPLE'.table'
rm $SAMPLE'.mark.bam'

