#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

NORMAL=$1
TUMOR=$2

python2.7 /manta/bin/configManta.py \
	--exome \
	--normalBam 'BQRS/'$NORMAL'.ready.bam' \
	--tumorBam 'BQRS/'$TUMOR'.ready.bam' \
	--referenceFasta $GENOME \
	--runDir 'Manta/'$NORMAL'-'$TUMOR
# execution on a single local machine with 22 parallel jobs
'Manta/'$NORMAL'-'$TUMOR/runWorkflow.py -m local -j 2

