#!/bin/bash

# alorenzetti 20180616
# second round of nanopore analysis
# im doing everything again

threads=15
inputdir="20180209_2242_20180209_hsalinarum_genomes_run"
outputdir="20180816_basecalled"

# starting basecalling
read_fast5_basecaller.py -r -i $inputdir \
                         -t $threads \
	                 -f FLO-MIN106 \
			 -k SQK-LSK108 \
			 --barcoding \
			 -o fastq -q 0 -n 0 -s $outputdir \
			 --disable_filtering  
