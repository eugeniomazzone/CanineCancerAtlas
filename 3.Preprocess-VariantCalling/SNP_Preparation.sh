#!/bin/bash

mkdir Mutect_OP/ 
ls Mutect/ | grep -v tbi | while read line; do bcftools view -f PASS Mutect/$line > Mutect_OP/$line; done
ls Mutect_OP/ | cut -f1 -d'_' | while read line; do  mv Mutect_OP/$line'_fi_.vcf.gz' Mutect_OP/$line.vcf; done
#ls Mutect_OP/ | cut -f1 -d'.' | uniq |while read line  ; do bgzip Mutect_OP/$line.vcf; tabix Mutect_OP/$line.vcf.gz ; done

mkdir Mutect_Reheader/
cat Couples.txt | while read line
do
	bcftools view -s $( echo $line | cut -f1 -d',' ) Mutect_OP/$( echo $line | cut -f2 -d',' )'.vcf' > samp1.vcf 
	bcftools view -s $( echo $line | cut -f2 -d',' ) Mutect_OP/$( echo $line | cut -f2 -d',' )'.vcf' > samp2.vcf 
	bgzip samp1.vcf && bgzip samp2.vcf 
	tabix samp1.vcf.gz && tabix samp2.vcf.gz 
	bcftools merge samp1.vcf.gz samp2.vcf.gz >  Mutect_Reheader/$( echo $line | cut -f2 -d',' ).vcf
	bgzip Mutect_Reheader/$( echo $line | cut -f2 -d',' ).vcf
	tabix Mutect_Reheader/$( echo $line | cut -f2 -d',' ).vcf.gz
	rm samp1.vcf.gz && rm samp1.vcf.gz.tbi
	rm samp2.vcf.gz && rm samp2.vcf.gz.tbi
done

mkdir Strelka_SNPInDel/
ls Strelka/ | while read line
do
	bcftools concat -a Strelka/$line/results/variants/somatic.indels.vcf.gz Strelka/$line/results/variants/somatic.snvs.vcf.gz > Strelka_SNPInDel/$line.vcf
done

mkdir Strelka_OP/
ls Strelka_SNPInDel/ | while read line
do
	bcftools view -f PASS Strelka_SNPInDel/$line > Strelka_OP/$line
done

mkdir Strelka_Reheader/
ls Strelka_OP/ | cut -f1 -d'.' | uniq | while read line
do
	echo -e $(echo $line | cut -f1 -d'-' )'\n'$(echo $line | cut -f2 -d'-') > strelka.head.txt
	bcftools reheader -s strelka.head.txt Strelka_OP/$line.vcf > Strelka_Reheader/$line.vcf
	mv Strelka_Reheader/$line.vcf Strelka_Reheader/$(echo $line | cut -f2 -d'-').vcf
done

ls Strelka_Reheader/ | cut -f1 -d'.' | uniq | while read line
do
	bgzip Strelka_Reheader/$line.vcf
	tabix Strelka_Reheader/$line.vcf.gz
done

rm $( ls Varscan/ | while read line ; do ls Varscan/$line'/' | while read line2 ; do echo Varscan/$line'/'$line2; done; done | grep -v Somatic )
rm $( ls Varscan/ | while read line ; do ls Varscan/$line'/' | while read line2 ; do echo Varscan/$line'/'$line2; done; done | grep hc )
ls Varscan/ | while read line ; do ls Varscan/$line'/' | while read line2 ; do bgzip Varscan/$line'/'$line2; tabix Varscan/$line'/'$line2'.gz'; done; done

mkdir Varscan_SNPInDel/
ls Varscan/ | while read line
do
	bcftools concat -a Varscan/$line/$line.vcf.indel.Somatic.gz Varscan/$line/$line.vcf.snp.Somatic.gz > Varscan_SNPInDel/$line.vcf
done

mkdir Varscan_OP/
ls Varscan_SNPInDel/ | while read line ; do bcftools view -f PASS Varscan_SNPInDel/$line > Varscan_OP/$line ; done

mkdir Varscan_Reheader/
ls Varscan_OP/ | cut -f1 -d'.' | while read line; do echo -e $(echo $line | cut -f1 -d'-' )'\n'$(echo $line | cut -f2 -d'-') > strelka.head.txt; bcftools reheader -s strelka.head.txt Varscan_OP/$line.vcf > Varscan_Reheader/$line.vcf; mv Varscan_Reheader/$line.vcf Varscan_Reheader/$(echo $line | cut -f2 -d'-').vcf ; done
ls Varscan_Reheader/| cut -f1 -d'.' | uniq  | while read line ; do bgzip Varscan_Reheader/$line.vcf; tabix Varscan_Reheader/$line.vcf.gz; done


