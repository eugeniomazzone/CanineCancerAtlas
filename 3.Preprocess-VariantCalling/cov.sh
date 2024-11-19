#!/bin/bash

genome='/refFiles/genome/canFam3.fa.gz'
vcf='/refFiles/vcf1/Filtred_Published1.vcf.gz'
bed='/refFiles/refined_bed.bed'

mkdir StatsMDups
mkdir Stats
mkdir Depth
ls BQRS/ | cut -f1 -d. | uniq |  while read line 
do
        echo -e '##############' $line '#############'
        samtools stats BQRS/$line'.ready.bam' -t $bed -c 0,1000,1 > Stats/$line'.stat.txt'
        samtools view -L $bed -b BQRS/$line'.ready.bam' | samtools flagstat -  > StatsMDups/$line'.flagstat.txt'
	samtools depth -b $bed BQRS/$line'.ready.bam' > Depth/$line'.depth.txt'
done

