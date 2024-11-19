#!/bin/bash

echo -e Sample"\t"Tumor Depth"\t"Normal Depth > Depth.txt
cat /WorkDir/Couples.txt | while read line 
do
	NORMAL=$( echo $line | cut -f1 -d,)
	TUMOR=$(echo $line | cut -f2 -d,)
	echo -e $TUMOR"\t"$(cat /WorkDir/Depth/$NORMAL'.depth.txt' | awk '{sum=sum+$3} END{print sum/NR}')"\t"$(cat /WorkDir/Depth/$TUMOR'.depth.txt' | awk '{sum=sum+$3} END{print sum/NR}')
done >>Depth.txt
