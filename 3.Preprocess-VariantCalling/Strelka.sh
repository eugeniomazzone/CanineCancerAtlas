#!/bin/bash

GENOME='/genome/canFam3.fa.gz'
VCF='/vcf1/Filtred_Published1.vcf.gz'

NORMAL=$1
TUMOR=$2

python2.7 strelka/bin/configureStrelkaSomaticWorkflow.py \
	--exome \
	--NORMALBam 'BQRS/'$NORMAL \
	--TUMORBam 'BQRS/'$TUMOR \
	--referenceFasta $GENOME \
	--runDir 'Strelka/'$NORMAL'-'$TUMOR \ 
	--indelCandidates 'Manta/'$NORMAL'-'$TUMOR'/results/variants/candidateSmallIndels.vcf.gz'
# execution on a single local machine with 20 parallel jobs
'Strelka/'$NORMAL'-'$TUMOR/runWorkflow.py -m local -j 20

