#!/bin/bash

genome='/genome/canFam3.fa.gz'
vcf='/vcf1/Filtred_Published1.vcf.gz'

mkdir Align
mkdir Convert
mkdir Sorted
mkdir Groups
mkdir MDups
mkdir BQRS

IDENTIFIER=28
BQRS_CONT=$(ls BQRS/ | cut -f1 -d.)
cat SraRunTable_test.csv | sed -n '2,$ p' | cut -f1,$IDENTIFIER -d, | grep DLBCL | sed 's/cf_DLBCL/Tumor,/g' > AllSample.txt
IDENTIFIER=28
cat SraRunTable_test.csv | sed -n '2,$ p' | cut -f1,$IDENTIFIER -d, | grep Punch | sed 's/cf_Punch/Normal,/g' >> AllSample.txt

cat AllSample.txt | awk -F, '{print $3}' | sort | uniq | while read line 
do 
	Normal=$(cat AllSample.txt | grep Normal | awk -F, ' $3=="'$line'" {print $1}' ) 
	Tumor=$(cat AllSample.txt | grep Tumor | awk -F, ' $3=="'$line'" {print $1}') 
	echo $Normal,$Tumor   
done > Couples.txt

cat AllSample.txt | cut -f1 -d, | while read line
do
	SAMPLE=$line
	
	if [[ "$BQRS_CONT" =~ $SAMPLE ]] ;then
		echo ==================== $SAMPLE is already processed ====================  
	else
		echo ==================== $SAMPLE will NOW be processed ==================== 
		sh PreProcess.sh $SAMPLE
	fi
done
