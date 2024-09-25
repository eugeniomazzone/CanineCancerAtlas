#!/bin/bash

mkdir temp
mkdir temp_annovar
mkdir SNP

ls SNP_raw/ | while read line
do 
	bcftools norm -d none SNP_raw/$line > temp/$line
done

ls temp/ | cut -d'.' -f1 | while read line
do 
	/home/$($whoami)/annovar/annovar/convert2annovar.pl --format vcf4 -allsample -withfreq --withfilter --keepindelref --includeinfo temp/$line'.vcf' > temp_annovar/$line.vcf  
done
ls temp_annovar/ | cut -d'.' -f1 | uniq | while read line 
do 
	/home/$($whoami)/annovar/annovar/table_annovar.pl -buildver canFam3 -protocol refGene -operation g  temp_annovar/$line.vcf /home/$($whoami)/annovar/annovar/dogdb/dogdb_canFam3/  
	mv temp_annovar/$line.vcf.canFam3_multianno.txt SNP/ 
done

rm -rf temp
rm -rf temp_annovar
