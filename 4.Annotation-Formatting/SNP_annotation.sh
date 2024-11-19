#!/bin/bash

mkdir temp
mkdir temp_filt
mkdir temp_annovar
mkdir SNP

ls /WorkDir/SNP_raw/ | while read line
do 
	bcftools norm -d none /WorkDir/SNP_raw/$line > temp/$line
done

ls temp | while read FILE
do
	COV=$(cat NRPCC.tsv | grep $(echo $FILE | cut -f2 -d'-' | cut -f1 -d.) | awk '{print $NF}')
	eval bcftools view -H temp/$FILE -i "'"FORMAT/DP'>'$COV"'" '>' temp_filt/$FILE
done

ls temp_filt/ | cut -d'.' -f1 | while read line
do 
	/refFiles/annovar/convert2annovar.pl --format vcf4 -allsample -withfreq --withfilter --keepindelref --includeinfo temp_filt/$line'.vcf' > temp_annovar/$line.vcf  
done
ls temp_annovar/ | cut -d'.' -f1 | uniq | while read line 
do 
	/refFiles/annovar/table_annovar.pl -buildver canFam3 -protocol refGene -operation g  temp_annovar/$line.vcf /refFiles/annovar/dogdb/  
	mv temp_annovar/$line.vcf.canFam3_multianno.txt SNP/ 
done

rm -rf temp
rm -rf temp_annovar
