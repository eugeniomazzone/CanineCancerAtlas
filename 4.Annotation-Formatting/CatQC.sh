#!/bin/bash

cat NRPCC_manualAnno.tsv | awk -F'\t' '$10!~"-" {print $1"\t"$NF}' | while read line
do
	QCFAIL=$( echo $line | cut -f2 -d' ' )
	NAME=$( echo $line | cut -f1 -d' ' )
	FILE=$(ls SNP/* | grep $NAME )
	
	if [[ "$QCFAIL" =~ "Tumor" ]]; then
		mv $FILE $FILE.outQC
	fi
	if [[ "$QCFAIL" =~ "Normal" ]]; then
                mv $FILE $FILE.outQC
        fi
	if [[ "$QCFAIL" =~ "Metastasis" ]]; then
                mv $FILE $FILE.meta
        fi
	if [[ "$QCFAIL" =~ "Duplicate" ]]; then
                mv $FILE $FILE.dups
        fi
done
