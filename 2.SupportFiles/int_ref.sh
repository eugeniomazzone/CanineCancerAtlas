#!/bin/bash

wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/genes/canFam3.ensGene.gtf.gz

wget http://ftp.ensembl.org/pub/release-104/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.CanFam3.1.104.gtf.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/database/cytoBandIdeo.txt.gz && gzip -d cytoBandIdeo.txt.gz

mkdir vcf1 
cd vcf1 && wget ftp://download.big.ac.cn/idog/dogsd/vcf/Filtred_Published.vcf.bz2 && bzip2 -d Filtred_Published.vcf.bz2 && gzip Filtred_Published.vcf && tabix Filtred_Published.vcf.gz

mkdir genome
cd genome/ && wget https://hgdownload.soe.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz && gzip -d canFam3.fa.gz && bgzip canFam3.fa
bwa index genome/canFam3.fa.gz
samtools faidx genome/canFam3.fa.gz 
gatk CreateSequenceDictionary -R genome/canFam3.fa.gz

Rscript make_bed.r

wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz && \
tar -xvzf annovar.latest.tar.gz && rm -f annovar.latest.tar.gz
cd annovar && \
./annotate_variation.pl --downdb -buildver canFam3 refGene dogdb && \
./annotate_variation.pl --buildver canFam3 --downdb seq dogdb/canFam3_seq && \
./retrieve_seq_from_fasta.pl --format refGene --seqfile dogdb/canFam3_seq/canFam3.fa dogdb/canFam3_refGene.txt  -outfile dogdb/canFam3_refGeneMrna.fa 

