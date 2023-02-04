#!/bin/bash

for cell in b mono nk t

do

awk 'BEGIN{FS=OFS="\t"} {gsub("chr","",$1); print "chr"$1,$2,$3,$1":"$2"-"$3}' pbmc.sorted.merged.bare.bed \
        | bedtools intersect -a - -b tagalign/${cell}.tagAlign -wa -wb \
        | cut -f 4,8 | sort -S 96G -T `pwd` \
        | uniq -c \
        | awk 'BEGIN{OFS="\t"; print "peak","barcode","value"} {print $2,$3,$1}' \
        | gzip -c \
        > longf_mats/${cell}.merged_peaks.long_fmt.mtx.gz

done


