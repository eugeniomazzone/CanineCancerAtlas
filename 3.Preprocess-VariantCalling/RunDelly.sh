#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

mkdir Delly

cat Couples.txt | while read line 
do
        NORMAL=$( echo $line | cut -d, -f1)
        TUMOR=$( echo $line | cut -d, -f2)
        echo $NORMAL 
        echo $TUMOR
	#delly call -o 'Delly/'$TUMOR.bcf -g $GENOME 'BQRS/'$TUMOR'.ready.bam' 'BQRS/'$NORMAL'.ready.bam' 
done

mkdir 'PostDelly/'

cat Couples.txt | while read line
do 
	echo -e $(echo $line | cut -f1 -d,)'\t''control''\n'$(echo $line | cut -f2 -d,)'\t''tumor' > sample.tsv; 
	delly filter -f somatic -o 'PostDelly/'$(echo $line | cut -f2 -d,).bcf -s sample.tsv 'Delly/'$(echo $line | cut -f2 -d,).bcf ;
	bcftools view 'PostDelly/'$(echo $line | cut -f2 -d,).bcf > 'PostDelly/'$(echo $line | cut -f2 -d,).vcf
	rm 'PostDelly/'$(echo $line | cut -f2 -d,).bcf
	rm 'PostDelly/'$(echo $line | cut -f2 -d,).bcf.csi
done


mkdir 'Dellyfiltered/'
ls 'PostDelly' | grep -v csi | while read SAMPLE
do
	bcftools view -e INFO/IMPRECISE=1 -f PASS 'PostDelly'/$SAMPLE > 'Dellyfiltered/'$(echo $SAMPLE | cut -f1 -d.).vcf
done

