#!/bin/bash

GENOME='/refFiles/genome/canFam3.fa.gz'
VCF='/refFiles/vcf1/Filtred_Published1.vcf.gz'

NORMAL=$1
TUMOR=$2

python2.7 /strelka/bin/configureStrelkaSomaticWorkflow.py \
	--exome \
	--normalBam 'BQRS/'$NORMAL'.ready.bam' \
	--tumorBam 'BQRS/'$TUMOR'.ready.bam' \
	--referenceFasta $GENOME \
	--runDir 'Strelka/'$NORMAL'-'$TUMOR \
	--indelCandidates 'Manta/'$NORMAL'-'$TUMOR'/results/variants/candidateSmallIndels.vcf.gz'
# execution on a single local machine with 20 parallel jobs
'Strelka/'$NORMAL'-'$TUMOR/runWorkflow.py -m local -j 2

