#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

mkdir Mutect

cat Couples.txt | while read COUPLE
do
	NORMAL=$(echo $COUPLE | cut -f1 -d, )
	TUMOR=$(echo $COUPLE | cut -f2 -d, )
	sh Mutect.sh $NORMAL $TUMOR
done

./endMutect.sh
