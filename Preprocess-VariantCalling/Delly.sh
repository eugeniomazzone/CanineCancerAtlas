#!/bin/bash

GENOME='../genome/canFam3.fa.gz'
VCF='../vcf1/Filtred_Published1.vcf.gz'

mkdir Delly

cat Couples.txt | while read line 
do
        NORMAL=$( echo $line | cut -d, -f1)
        TUMOR=$( echo $line | cut -d, -f2)
        echo $NORMAL 
        echo $TUMOR
	delly call -o 'Delly/'$TUMOR.bcf -g $GENOME 'BQRS/'$TUMOR 'BQRS/'$NORMAL 
done
