#! /bin/bash

gtf='../canFam3.ensGene.gtf.gz'


mkdir Delly_sansa_bcf
mkdir Delly_sansa_tsv
ls Dellyfiltered | grep -v csi | while read SAMPLE
do
	sansa annotate -t 0 -r 0 -f CDS -i gene_id -c -g $gtf \
		-a Delly_sansa_bcf/$(echo $SAMPLE | cut -f1 -d.).vcf \
		-o Delly_sansa_tsv/$(echo $SAMPLE | cut -f1 -d.).tsv.gz \
       		Dellyfiltered/$(echo $SAMPLE | cut -f1 -d.).vcf 

done
