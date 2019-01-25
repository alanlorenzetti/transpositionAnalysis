#!/bin/bash

barcodes=`echo barcode{01..06}`
deepbinnerdir=20180817_deepbinner
porechopdir=20180817_porechop
outputdir=20181002_filtlong_stats

if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; fi

for i in $barcodes ; do
  filtlong --min_mean_q 1 \
           --verbose \
           ${porechopdir}/${i}.fastq.gz 2>&1 > /dev/null |\
           grep -P "length = \d+"  |\
           sed 's/\s\+length = //' |\
           sed 's/\s\+\w\+ quality = /\t/g' |\
           awk -v barcode=$i -v OFS="\t" '{print barcode,$0}' > ${outputdir}/${i}_stats.txt
done

filtlong --min_mean_q 1 \
         --verbose \
         ${deepbinnerdir}/unclassified.fastq.gz 2>&1 > /dev/null |\
         grep -P "length = \d+"  |\
         sed 's/\s\+length = //' |\
         sed 's/\s\+\w\+ quality = /\t/g' |\
         awk -v barcode=unclassified -v OFS="\t" '{print barcode,$0}' > ${outputdir}/unclassified_stats.txt

echo -e "barcode\tlength\tmeanQuality\twindowQuality" > ${outputdir}/stats.txt
cat ${outputdir}/*_stats.txt >> ${outputdir}/stats.txt

Rscript --slave --args $outputdir run_filtlong_stats_ggplot2_20181002.R
