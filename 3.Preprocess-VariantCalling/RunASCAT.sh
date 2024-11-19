#! /bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

mkdir CNA

cat Couples.txt | while read line 
do
        NORMAL=$( echo $line | cut -d, -f1)
        TUMOR=$( echo $line | cut -d, -f2)
        ls BQRS | grep $NORMAL 
        ls BQRS | grep $TUMOR
	Rscript ASCAT.r $NORMAL $TUMOR 
done
