#! /bin/bash

biop=$1

mkdir SNP_raw
ls Mutect_Reheader | cut -f1 -d'.' | while read line
do 
bcftools isec -n+2 Mutect_Reheader/$line.vcf.gz Strelka_Reheader/$line.vcf.gz Varscan_Reheader/$line.vcf.gz -p SNP_raw

rm SNP_raw/0000.vcf.gz
rm SNP_raw/0001.vcf.gz
rm SNP_raw/0002.vcf.gz

rm SNP_raw/0000.vcf.gz.tbi
rm SNP_raw/0001.vcf.gz.tbi
rm SNP_raw/0002.vcf.gz.tbi

bgzip SNP_raw/0000.vcf
bgzip SNP_raw/0001.vcf
bgzip SNP_raw/0002.vcf

tabix SNP_raw/0000.vcf.gz
tabix SNP_raw/0001.vcf.gz
tabix SNP_raw/0002.vcf.gz

bcftools concat -a SNP_raw/0000.vcf.gz SNP_raw/0002.vcf.gz > SNP_raw/$line.vcf

done

rm SNP_raw/0000.vcf.gz
rm SNP_raw/0001.vcf.gz
rm SNP_raw/0002.vcf.gz

rm SNP_raw/0000.vcf.gz.tbi
rm SNP_raw/0001.vcf.gz.tbi
rm SNP_raw/0002.vcf.gz.tbi

rm SNP_raw/README.txt
rm SNP_raw/sites.txt
