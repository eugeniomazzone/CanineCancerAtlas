#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

ls Pon/ | grep -v tbi | while read line 
do
	echo -e $( echo $line | cut -f1 -d'.' )'\t''Pon/'$line 
done > gen_pon.list

mkdir tmp1
rm -fr pon_db

gatk GenomicsDBImport -R $GENOME -L 'chr.list'  --genomicsdb-workspace-path 'pon_db' --sample-name-map 'gen_pon.list' --tmp-dir 'tmp1'
gatk CreateSomaticPanelOfNormals -R $GENOME -V gendb://pon_db -O 'pon.vcf.gz'

