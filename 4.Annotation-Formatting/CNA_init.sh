#!/bin/bash

mkdir TCN
cat GoFresults.txt | grep -v ^gamma | cut -f9 -d'	' | while read line
do
	cp $line TCN/
done
