#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

mkdir Pon

cat AllSample.txt | awk -F, ' $2=="Normal"{print $1}' | while read NORMAL
do
	echo $NORMAL
	sh PON.sh $NORMAL
done

./mergeChPON.sh
./endPON.sh

rm -f $(ls *_SRR*)
