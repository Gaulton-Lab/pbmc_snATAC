#!/usr/bin/env python3

import os
import gzip
import argparse
import subprocess
import logging
import numpy as np
import pandas as pd
from multiprocessing import Pool


def load_clusters(args):
	df = pd.read_table(args.cluster_file, sep='\t', header=0, index_col=0)
	cluster_map = df[args.cluster_col].to_dict()
	unique_clusters = sorted(set(df[args.cluster_col]))
	pass_barcodes = df.index
	return cluster_map, unique_clusters, pass_barcodes

def split_reads_and_call_peaks(args, c):
	cluster_prefix  = '{}.{}'.format(args.output_prefix, c)
	cluster_tagalign = '{}.{}'.format(cluster_prefix, 'bed')
	cluster_barcodes = [bc for bc in args.cluster_map if args.cluster_map[bc]==c]
	if not os.path.isfile(cluster_tagalign):
		with gzip.open(args.tagalign_file, 'rt') as f, open(cluster_tagalign, 'w') as f_out:
			for line in f:
				fields = line.rstrip('\n').split('\t')
				barcode = fields[3]
				if barcode in cluster_barcodes:
					print(line.rstrip('\n'), file=f_out)
	macs2_cmd = ['macs2', 'callpeak', '-t', cluster_tagalign, '--outdir', os.getcwd(), '-n', cluster_prefix, '-q', '0.05', '--nomodel', '--keep-dup', 'all', '-g', 'hs']
	with open(cluster_prefix + '.macs2_callpeak.log', 'w') as f:
		subprocess.call(macs2_cmd, stdout=f, stderr=f)

	total_reads = subprocess.check_output('cat {} | wc -l'.format(cluster_tagalign), shell=True)	
	genomecov_cmd = ['bedtools', 'genomecov', '-i', cluster_tagalign, '-bg', '-scale', '{}'.format(1e6/int(total_reads.decode())), '-g', args.genome_file]
	with open(output_prefix + '.scale_1e6.bdg', 'w') as f:
		subprocess.call(genomecov_cmd, stdout=f)
	subprocess.call(['bedGraphToBigWig', output_prefix + '.scale_1e6.bdg', args.genome_file, output_prefix + '.scale_1e6.ATAC.bw'])
	os.rename(cluster_tagalign, '{}.tagalign'.format(cluster_prefix))
	subprocess.call(['gzip', '{}.tagalign'.format(cluster_prefix)])
	os.remove(output_prefix + '_summits.bed')
	os.remove(output_prefix + '_peaks.xls')
	return

def merge_peaks(args):
	peak_files = ['{}.{}_peaks.narrowPeak'.format(args.output_prefix, c) for c in args.unique_clusters]
	allpeaks_file = '{}.all_peaks.anno.bed'.format(args.output_prefix)
	mergedpeaks_file = '{}.merged_peaks.anno.bed'.format(args.output_prefix)
	
	with open(allpeaks_file, 'w') as apf:
		for peak_file in peak_files:
			cluster_name = peak_file.split('.')[1].split('_peaks')[0]
			with open(peak_file) as pf:
				for line in pf:
					fields = line.rstrip().split('\t')
					print(fields[0], fields[1], fields[2], cluster_name, sep='\t', file=apf)
	
	subprocess.call(['sort', '-k1,1', '-k2,2n', '-o', allpeaks_file, allpeaks_file])
	
	with open(mergedpeaks_file, 'w') as mpf:
		merge1 = subprocess.Popen(['bedtools', 'merge', '-i', allpeaks_file], stdout=subprocess.PIPE)
		intersect = subprocess.Popen(['bedtools', 'intersect', '-a', '-', '-b', allpeaks_file, '-wa', '-wb'], stdin=merge1.stdout, stdout=subprocess.PIPE)
		merge2= subprocess.call(['bedtools', 'merge', '-i' ,'-', '-c', '7', '-o', 'distinct'], stdin=intersect.stdout, stdout=mpf)
	return mergedpeaks_file

def create_peak_matrix(args):
	lfmatrix_file = '{}.lf_mtx.gz'.format(args.output_prefix)

	intersect = subprocess.Popen(['bedtools', 'intersect', '-a', args.mergedpeaks_file, '-b', args.tagalign_file, '-wa', '-wb', '-sorted'], stdout=subprocess.PIPE)
	with gzip.open(lfmatrix_file, 'wt') as lf:
		for line in intersect.stdout:
			fields = line.decode().rstrip().split('\t')
			print('{}:{}-{}\t{}'.format(fields[0].replace('chr',''), fields[1], fields[2], fields[7]), file=lf)
	lfm = pd.read_table(lfmatrix_file, sep='\t', header=None, names=['peak','cell'])
	lfm = lfm.loc[lfm['cell'].isin(args.pass_barcodes)]
	lfm = lfm.groupby(lfm.columns.tolist()).size().reset_index().rename(columns={0:'value'})
	lfm.to_csv(lfmatrix_file, sep='\t', index=False, compression='gzip')
	
	mtx_file = '{}.mtx'.format(args.output_prefix)
	barcodes_file = '{}.barcodes'.format(args.output_prefix)
	peaks_file = '{}.peaks'.format(args.output_prefix)

	tmp_R = '{}.tmp.R'.format(args.output_prefix)
	with open(tmp_R, 'w') as tR:
		print('''library(Matrix)''', file=tR)
		print('''data <- read.table('{}', sep='\\t', header=TRUE)'''.format(lfmatrix_file), file=tR)
		print('''sparse.data <- with(data, sparseMatrix(i=as.numeric(peak), j=as.numeric(cell), x=value, dimnames=list(levels(peak), levels(cell))))''', file=tR)
		print('''t <- writeMM(sparse.data, '{}')'''.format(mtx_file), file=tR)
		print('''write.table(data.frame(colnames(sparse.data)), file='{}', col.names=FALSE, row.names=FALSE, quote=FALSE)'''.format(barcodes_file), file=tR)
		print('''write.table(data.frame(rownames(sparse.data)), file='{}', col.names=FALSE, row.names=FALSE, quote=FALSE)'''.format(peaks_file), file=tR)
	subprocess.call(['Rscript', tmp_R])
	subprocess.call(['gzip', mtx_file])
	os.remove(tmp_R)
	return

def find_marker_peaks(args):
	import scipy.io
	import statsmodels.api as sm
	
	clusters = pd.read_table(args.cluster_file, sep='\t', header=0, index_col=0)
	mtx = scipy.io.mmread('{}.mtx.gz'.format(args.output_prefix)).tocsr()
	barcodes = open('{}.barcodes'.format(args.output_prefix)).read().splitlines()
	peaks = open('{}.peaks'.format(args.output_prefix)).read().splitlines()
	zscores = pd.DataFrame(columns=sorted(set(clusters[args.cluster_col])))
	chunks = np.arange(0, mtx.shape[0], 10000)
	chunks = np.append(chunks, mtx.shape[0])
	for i in range(len(chunks)-1):
		z = pd.DataFrame(columns=sorted(set(clusters[args.cluster_col])))
		dmtx = pd.DataFrame(mtx[chunks[i]:chunks[i+1],:].todense(), columns=barcodes, index=peaks[chunks[i]:chunks[i+1]])
		dmtx = (dmtx > 0).astype(int)
		for cluster in sorted(set(clusters[args.cluster_col])):
			label_cov = clusters.copy().loc[dmtx.columns]
			label_cov['OvR'] = [1 if c==cluster else -1 for c in label_cov[args.cluster_col]]
			label_cov = label_cov[['OvR'] + ['log_usable_counts']]
			label_cov = sm.add_constant(label_cov)
			z[cluster] = dmtx.apply(lambda x: sm.Logit(x, label_cov.values).fit().tvalues[1] if sum(x[dmtx.columns.isin(label_cov.loc[label_cov['OvR']==1].index)])>0 else np.nan, axis=1)
		zscores = zscores.append(z)
	zscores.to_csv('{}.logit_zscores.txt'.format(args.output_prefix), sep='\t', header=True, index=True, float_format='%.4f')
	return

def main(args):
	args.cluster_map, args.unique_clusters, args.pass_barcodes = load_clusters(args)
	arglist = [(args, c) for c in args.unique_clusters]
	if not args.skip_split:
		with Pool(processes=min(36, len(args.unique_clusters))) as pool:
			pool.starmap(split_reads_and_call_peaks, arglist)
	if not args.skip_merge:
		args.mergedpeaks_file = merge_peaks(args)
	if not args.skip_matrix:
		if args.skip_merge:
			args.mergedpeaks_file = '{}.merged_peaks.anno.bed'.format(args.output_prefix)
		create_peak_matrix(args)
	if args.find_marker_peaks:
		find_marker_peaks(args)
	return

def process_args():
	parser = argparse.ArgumentParser(description='Call peaks for each cluster from snATAC-seq data')
	io_group = parser.add_argument_group('I/O arguments')
	io_group.add_argument('-c', '--cluster-file', required=True, type=str, help='Path to cluster file')
	io_group.add_argument('-t', '--tagalign-file', required=True, type=str, help='Path to tagAlign file')
	io_group.add_argument('-n', '--cluster-col', required=False, default='cluster_name', type=str, help='Name of the cluster column')
	io_group.add_argument('-o', '--output-prefix', required=True, type=str, help='Output prefix to prepend')
	io_group.add_argument('-g', '--genome-file', required=False, default='/home/joshchiou/references/hg19.chrom.sizes', type=str, help='Output prefix to prepend')
	
	skip_group = parser.add_argument_group('Skip steps')
	skip_group.add_argument('--skip-split', required=False, action='store_true', default=False, help='Skip read split and peak calling step')
	skip_group.add_argument('--skip-merge', required=False, action='store_true', default=False, help='Skip peak merging step')
	skip_group.add_argument('--skip-matrix', required=False, action='store_true', default=False, help='Skip peak matrix generation step')
	skip_group.add_argument('--find-marker-peaks', required=False, action='store_true', default=False, help='Find marker peaks for each cluster')
	return parser.parse_args()

if __name__ == '__main__':
	logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.DEBUG)
	args = process_args()
	main(args)
