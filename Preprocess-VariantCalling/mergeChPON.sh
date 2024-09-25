#!/bin/bash

GENOME='../genome/canFam3.fa.gz'
VCF='../vcf1/Filtred_Published1.vcf.gz'

cat AllSample.txt | awk -F, ' $2=="Normal"{print $1}' | while read NORMAL
do
	echo gatk MergeVcfs '\' > to_do.sh
	ls | grep $NORMAL | grep -vE 'stats|tbi|tsv' | while read file
	do
		echo -I $file '\' >> to_do.sh
	done
	echo -O 'Pon/'$NORMAL'.vcf.gz'  >> to_do.sh
	sh to_do.sh
	rm to_do.sh
done 
