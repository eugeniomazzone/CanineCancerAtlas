#! /bin/bash

biop=$1

mkdir $biop'_mut_OP' 
ls $biop'_mut/' | grep -v tbi | while read line; do bcftools view -f PASS $biop'_mut/'$line > $biop'_mut_OP/'$line; done
ls $biop'_mut_OP' | cut -f1 -d'_' | while read line; do  mv $biop'_mut_OP/'$line'_fi_.vcf.gz' $biop'_mut_OP/'$line.vcf; done
#ls $biop'_mut_OP/' | cut -f1 -d'.' | uniq |while read line  ; do bgzip $biop'_mut_OP/'$line.vcf; tabix $biop'_mut_OP/'$line.vcf.gz ; done

mkdir $biop'_mut_reheader'
ls $biop'_mut_OP' | cut -f1 -d'.' | uniq | while read line; do  bcftools view -s $( echo $line | cut -f1 -d'-' ) $biop'_mut_OP/'$line'.vcf' > samp1.vcf ; \
bcftools view -s $( echo $line | cut -f2 -d'-' ) $biop'_mut_OP/'$line'.vcf' > samp2.vcf ;  \
bgzip samp1.vcf; bgzip samp2.vcf ; \
tabix samp1.vcf.gz; tabix samp2.vcf.gz ; \
bcftools merge samp1.vcf.gz samp2.vcf.gz >  $biop'_mut_reheader/'$line.vcf; \
bgzip $biop'_mut_reheader/'$line.vcf ; \
tabix $biop'_mut_reheader/'$line.vcf.gz; \
rm samp1.vcf.gz; rm samp1.vcf.gz.tbi; \
rm samp2.vcf.gz && rm samp2.vcf.gz.tbi; done

#mv to_download $biop'_strelka'
mkdir $biop'_strelka_snpindel'
ls $biop'_strelka/' | while read line ; do bcftools concat -a $biop'_strelka/'$line/varinats/somatic.indels.vcf.gz $biop'_strelka/'$line/varinats/somatic.snvs.vcf.gz > $biop'_strelka_snpindel/'$line.vcf ; done
mkdir $biop'_strelka_OP'
ls $biop'_strelka_snpindel' | while read line ; do bcftools view -f PASS $biop'_strelka_snpindel/'$line > $biop'_strelka_OP/'$line ; done
mkdir $biop'_strelka_reheader'
ls $biop'_strelka_OP/' | cut -f1 -d'.' | uniq | while read line; do echo $(echo $line | cut -f1 -d'-' )'\n'$(echo $line | cut -f2 -d'-') > strelka.head.txt; bcftools reheader -s strelka.head.txt $biop'_strelka_OP/'$line.vcf > $biop'_strelka_reheader/'$line.vcf; done
ls $biop'_strelka_reheader/' | cut -f1 -d'.' | uniq | while read line  ; do bgzip $biop'_strelka_reheader/'$line.vcf; tabix $biop'_strelka_reheader/'$line.vcf.gz ; done

rm $( ls $biop'_varscan/' | while read line ; do ls $biop'_varscan/'$line'/' | while read line2 ; do echo $biop'_varscan/'$line'/'$line2; done; done | grep -v Somatic )
rm $( ls $biop'_varscan/' | while read line ; do ls $biop'_varscan/'$line'/' | while read line2 ; do echo $biop'_varscan/'$line'/'$line2; done; done | grep hc )
ls $biop'_varscan/' | while read line ; do ls $biop'_varscan/'$line'/' | while read line2 ; do echo $biop'_varscan/'$line'/'$line2; done; done

mkdir $biop'_varscan_snpindel'
ls $biop'_varscan/' | while read line ; do bcftools concat -a $biop'_varscan/'$line/$line.vcf.indel.Somatic.gz $biop'_varscan/'$line/$line.vcf.snp.Somatic.gz > $biop'_varscan_snpindel/'$line.vcf ; done
mkdir $biop'_varscan_OP'
ls $biop'_varscan_snpindel' | while read line ; do bcftools view -f PASS $biop'_varscan_snpindel/'$line > $biop'_varscan_OP/'$line ; done
mkdir $biop'_varscan_reheader/'
ls $biop'_varscan_OP/' | cut -f1 -d'.' | while read line; do echo $(echo $line | cut -f1 -d'-' )'\n'$(echo $line | cut -f2 -d'-') > strelka.head.txt; bcftools reheader -s strelka.head.txt $biop'_varscan_OP/'$line.vcf > $biop'_varscan_reheader/'$line.vcf; done
ls $biop'_varscan_reheader/'| cut -f1 -d'.' | uniq  | while read line ; do bgzip $biop'_varscan_reheader/'$line.vcf; tabix $biop'_varscan_reheader/'$line.vcf.gz; done
