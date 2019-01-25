#!/bin/bash

threads=22
inputdir="20180817_deepbinner"
outputdir="20180817_porechop"
barcodes=`echo barcode{01..06}`

if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; fi

for i in $barcodes ; do
	porechop -t $threads \
		 --format fastq.gz \
		 -i ${inputdir}/${i}.fastq.gz \
		 -o ${outputdir}/${i}.fastq.gz > ${outputdir}/${i}.log 2>&1	
done
