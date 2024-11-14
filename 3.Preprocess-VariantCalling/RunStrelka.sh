#!/bin/bash

GENOME='/genome/canFam3.fa.gz'
VCF='/vcf1/Filtred_Published1.vcf.gz'

mkdir Manta
mkdir Strelka

cat Couples.txt | while read line 
do
        NORMAL=$( echo $line | cut -d, -f1)
        TUMOR=$( echo $line | cut -d, -f2)
        echo '========' $NORMAL '========' $TUMOR '========'
        rm -fr 'Manta/'$NORMAL'-'$TUMOR
        sh Manta.sh $NORMAL $TUMOR
	rm -fr 'Strelka/'$NORMAL'-'$TUMOR
	sh Strelka.sh $NORMAL $TUMOR 
done 
