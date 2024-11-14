#!/bin/bash

GENOME='/genome/canFam3.fa.gz'
VCF='/vcf1/Filtred_Published1.vcf.gz'

cat Couples.txt | cut -f2 -d, | while read TUMOR 
do

	echo gatk MergeVcfs '\' > to_do.sh
	ls | grep 'fi.*'$TUMOR | grep -vE 'tbi|tsv' | while read file
	do
		echo -I $file '\' >> to_do.sh
	done
	echo -O 'Mutect/'$TUMOR'_fi_.vcf.gz'  >> to_do.sh
	sh to_do.sh
	rm to_do.sh
done 

echo '------------------------------------------'
