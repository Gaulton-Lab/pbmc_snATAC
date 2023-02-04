#!/bin/bash

for cell in B-cell Megakaryocyte Monocyte NK-cell pDC CD4-T-cell CD8-T-cell

do

awk -v OFS="\t" -v n="$cell" '{print $1,$2,$3, n }'  pbmc1.${cell}_peaks.narrowPeak > pbmc1.${cell}_peaks.bed

done

cat pbmc1.B-cell_peaks.bed pbmc1.Megakaryocyte_peaks.bed pbmc1.Monocyte_peaks.bed pbmc1.NK-cell_peaks.bed pbmc1.pDC_peaks.bed pbmc1.CD4-T-cell_peaks.bed pbmc1.CD8-T-cell_peaks.bed > pbmc1.bed
sort -k1,1 -k2,2n pbmc1.bed > pbmc1.sorted.bed

bedtools intersect -a pbmc1.sorted.bed -b /home/pbenaglio/general_files/ENCODE.hg19.blacklist.bed -v > pbmc1.sorted.filtered.bed
bedtools merge -i pbmc1.sorted.filtered.bed -c 4 -o collapse > pbmc1.sorted.merged.bed
awk -v OFS="\t" '{print $1,$2,$3 }' pbmc1.sorted.merged.bed > pbmc1.sorted.merged.bare.bed
