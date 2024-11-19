#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

NORMAL=$1
TUMOR=$2

samtools mpileup -q 1 -f $GENOME 'BQRS/'$NORMAL'.ready.bam' -o $NORMAL'.mpileup' & \ 
samtools mpileup -q 1 -f $GENOME 'BQRS/'$TUMOR'.ready.bam' -o $TUMOR'.mpileup'
wait
rm -r $NORMAL'-'$TUMOR
mkdir 'Varscan/'$NORMAL'-'$TUMOR

varscan somatic $NORMAL'.mpileup' \
$TUMOR'.mpileup' \
'Varscan/'$NORMAL'-'$TUMOR'/'$NORMAL'-'$TUMOR'.vcf' \
--output-vcf 1 

varscan processSomatic 'Varscan/'$NORMAL'-'$TUMOR'/'$NORMAL'-'$TUMOR'.vcf.indel'
varscan processSomatic 'Varscan/'$NORMAL'-'$TUMOR'/'$NORMAL'-'$TUMOR'.vcf.snp'

rm $NORMAL'.mpileup' 
rm $TUMOR'.mpileup'
