#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

NORMAL=$1

gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '1_'$NORMAL'.vcf.gz'  -L chr1 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '2_'$NORMAL'.vcf.gz'  -L chr2 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '3_'$NORMAL'.vcf.gz'  -L chr3 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '4_'$NORMAL'.vcf.gz'  -L chr4 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '5_'$NORMAL'.vcf.gz'  -L chr5 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '6_'$NORMAL'.vcf.gz'  -L chr6 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '7_'$NORMAL'.vcf.gz'  -L chr7 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '8_'$NORMAL'.vcf.gz'  -L chr8 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '9_'$NORMAL'.vcf.gz'  -L chr9 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '10_'$NORMAL'.vcf.gz' -L chr10 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '11_'$NORMAL'.vcf.gz' -L chr11 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '12_'$NORMAL'.vcf.gz' -L chr12 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '13_'$NORMAL'.vcf.gz' -L chr13 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '14_'$NORMAL'.vcf.gz' -L chr14 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '15_'$NORMAL'.vcf.gz' -L chr15 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '16_'$NORMAL'.vcf.gz' -L chr16 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '17_'$NORMAL'.vcf.gz' -L chr17 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '18_'$NORMAL'.vcf.gz' -L chr18 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '19_'$NORMAL'.vcf.gz' -L chr19 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '20_'$NORMAL'.vcf.gz' -L chr20 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '21_'$NORMAL'.vcf.gz' -L chr21 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '22_'$NORMAL'.vcf.gz' -L chr22 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '23_'$NORMAL'.vcf.gz' -L chr23 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '24_'$NORMAL'.vcf.gz' -L chr24 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '25_'$NORMAL'.vcf.gz' -L chr25 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '26_'$NORMAL'.vcf.gz' -L chr26 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '27_'$NORMAL'.vcf.gz' -L chr27 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '28_'$NORMAL'.vcf.gz' -L chr28 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '29_'$NORMAL'.vcf.gz' -L chr29 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '30_'$NORMAL'.vcf.gz' -L chr30 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '31_'$NORMAL'.vcf.gz' -L chr31 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '32_'$NORMAL'.vcf.gz' -L chr32 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '33_'$NORMAL'.vcf.gz' -L chr33 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '34_'$NORMAL'.vcf.gz' -L chr34 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '35_'$NORMAL'.vcf.gz' -L chr35 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '36_'$NORMAL'.vcf.gz' -L chr36 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '37_'$NORMAL'.vcf.gz' -L chr37 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O '38_'$NORMAL'.vcf.gz' -L chr38 & \ 
gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O 'X_'$NORMAL'.vcf.gz'  -L chrX 
#gatk Mutect2 -R $GENOME -I 'BQRS/'$NORMAL'.ready.bam'  -max-mnp-distance 0 -O 'UN_'$NORMAL'.vcf.gz' -XL chr.list
