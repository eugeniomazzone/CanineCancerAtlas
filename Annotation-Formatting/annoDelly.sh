#! /bin/bash

gtf='canFam3.ensGene.gtf.gz'

biop=$1

mkdir $biop'_annoDelly/'

ls $biop'_postDelly/' | while read line ;
do sansa annotate -i gene_id -f CDS --gtf $gtf -r 0 -t 10  -c $biop'_postDelly/'$line -o $biop'_annoDelly/'$( echo $line | cut -f1 -d'.').tsv.gz ;
done

