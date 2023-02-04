awk 'BEGIN{FS=OFS="\t"} {gsub("chr","",$1); print "chr"$1,$2,$3,$1":"$2"-"$3}' pbmc1.sorted.merged.bed \
	| bedtools intersect -a - -b ../pbmc_merged.filt.rmdup.tagAlign.gz -wa -wb \
	| cut -f 4,8 | sort -S 96G -T `pwd` \
	| uniq -c \
	| awk 'BEGIN{OFS="\t"; print "peak","barcode","value"} {print $2,$3,$1}' \
	| gzip -c \
	> pbmc1.merged_peaks.long_fmt.mtx.gz

