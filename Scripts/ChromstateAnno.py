#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''
Description:
Annotating OGs with chromatin state

Chromatin state combined segmentations produced by chromHMM and Segway of HeLa cells were retrieved from UCSC at 
https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=wgEncodeAwgSegmentation&db=hg19

'''


import pysam
import sys, os
import argparse
import gzip
import py2bit
import numpy as np
import pandas as pd
import matplotlib
from collections import defaultdict
from matplotlib import pyplot as plt




def get_arguments():

	parser = argparse.ArgumentParser()

	parser.add_argument('-a','--group_a', nargs = '+', type=str, help=
		'Group A bamfiles, usually treat sample bam.')

	parser.add_argument('-b','--group_b', nargs = '+', type=str, help=
		'Group A bamfiles, usually control sample bam.')

	parser.add_argument('--chromstate', help = 'Chromatin State file. '
		'e.g. wgEncodeAwgSegmentationCombinedHelas3.bed.gz')

	parser.add_argument('-g','--genome',help='Genome in two bit format. [.2bit]')

	parser.add_argument('-o','--outfigure', help='File name to save figure.')

	parser.add_argument('--outmatrix', help='File name to save results.')

	parser.add_argument('-t','--title', help='Set title for the figure. default is None')

	return parser



def get_base_content(chrom,start,end,tbit,fraction=True, base = 'G'):

	# check if base valid
	base = base.upper().strip()
	for i in range(len(base)):
		try:
			assert base[i] in ['A','C','G','T']
		except AssertionError:
			sys.stderr.write('Invalid base type.\n')
			return None

	bases = tbit.bases(chrom,start,end,fraction=False)
	if end > tbit.chroms(chrom):
		end = tbit.chroms(chrom)
	if sum(bases.values()) < 0.95 * (start - end):
		raise Exception("WARNING: too many NNNs present in {}:{}-{}".format(chrom, start, end))
		return None

	if len(base) == 1:
		if fraction:
			return (bases[base]) / float(end - start)
		return bases[base]

	else:
		if fraction:
			return sum([bases[base[i]] for i in range(len(base))]) / float(end - start)
		return sum([bases[base[i]] for i in range(len(base))])



def _cal_read_counts(a_bam_list, b_bam_list, annofile, outfile, tb):

	a_bamfiles = [ pysam.AlignmentFile(bam) for bam in a_bam_list ]
	b_bamfiles = [ pysam.AlignmentFile(bam) for bam in b_bam_list ]

	fo = open(outfile,'w')
	fo.write('\t'.join(map(str, ['Chrom','Start','End','Name','gA_mean_counts','gB_mean_counts', 'GC_ratio'])))
	fo.write('\n')


	with gzip.open(annofile) as fi:

		for line in fi:
			line = line.decode().rstrip().split('\t')

			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			anno_type = line[3]

			a_counts, b_counts = [], []
			for bf in a_bamfiles:
				a_counts.append(
					bamCountinRegion(bf, chrom, start, end, offset=5, strand=None))

			for bf in b_bamfiles:
				b_counts.append(
					bamCountinRegion(bf, chrom, start, end, offset=5, strand=None))

			a_mean = np.mean(a_counts)
			b_mean = np.mean(b_counts)

			GCratio = get_base_content(chrom, start, end, tb, fraction=True, base='GC')

			out_line = '\t'.join(map(str, [chrom, start, end, anno_type, a_mean, b_mean, GCratio]))
			print(out_line, file=fo)


	# close bamfiles
	for bf in a_bamfiles:
		bf.close()

	for bf in b_bamfiles:
		bf.close()


	return



def bamCountinRegion(bf,chrom,start,end,offset=0,strand=None):

	if start >= end:
		sys.exit("bamCountinRegion Error:\n\tstart >= end.")
	if start <0:
		sys.stderr.write('{}:{}-{}\nInput start site less than zero.'.format(chrom,start,end))
		return 0

	chrmax = bf.lengths[bf.references.index(chrom)]

	if end > chrmax:
		end = chrmax
	if start >= chrmax:
		return 0

	reads = bf.fetch(chrom,start,end)
	p, m = 0, 0
	for read in reads:
		if read.is_unmapped:
			continue
		if read.is_reverse:
			if read.reference_end - offset > end or read.reference_end - offset < start :
				continue
			else:
				m += 1
		else:
			if read.reference_start + offset < start or read.reference_start + offset >end:
				continue
			else:
				p += 1

	if strand == '+':
		return p
	elif strand == '-':
		return m
	else:
		a = p + m
		return a



	
def get_color(name):

	color_pairs = {
		"TSS":"#FF0000",
		"PF":"#FA8072",
		"E":"#FFA500",
		"WE":"#FFFF00",
		"CTCF":"#0000FF",
		"T":"#006400",
		"R":"#BEBEBE"
	}
	select = color_pairs.get(name)

	return select




def boxplot(data, title, outfig):

	sort_key = ['TSS','PF','E','WE','CTCF','T','R']

	# sort input data
	resort_dat = []
	colors = []

	for idx in range(len(sort_key)):
		resort_dat.append(data.get(sort_key[idx]))
		colors.append(get_color(sort_key[idx]))

	fig = plt.figure(figsize=(9,6),dpi=300)
	ax = fig.add_subplot(111)

	bp = ax.boxplot(resort_dat,patch_artist=True,labels=None,
		flierprops={'markerfacecolor':'black','marker':','},notch=0,whis=0.3,showcaps=True,showfliers=False)
	plt.setp(bp['whiskers'], color='black', linestyle='dashed',linewidth=0.8)
	plt.setp(bp['medians'], color='black',linewidth=0.6)

	[bp['boxes'][i].set(facecolor=colors[i], alpha=0.7) for i in range(len(resort_dat))]

	ax.set_xticks(range(1,len(resort_dat)+1))
	ax.set_xticklabels(sort_key,rotation=0,fontsize=12,fontweight='bold')
	ax.set_ylabel(r"Normalized OG signal",fontsize=10,fontstyle='normal',fontweight='bold')
	ax.set_xlabel("")
	ax.set_title(title,fontsize=15)
	ax.set_ylim(-2.5,2.5)
	ax.set_yticks([-2.5,0,2.5])
	ax.set_yticklabels([-2.5,0,2.5])
	plt.axhline(y=0,c='red',lw=1.5,ls='--')

	plt.savefig(outfig, bbox_inches='tight', dpi=300)
	plt.close()

	return



def main():

	args = get_arguments().parse_args()

	
	group_a = args.group_a
	group_b = args.group_b

	annofile = args.chromstate
	tb = py2bit.open(args.genome)
	outfile = args.outmatrix

	# output results
	_cal_read_counts(group_a, group_b, annofile, outfile, tb)

	tb.close()

	# load previous results and plot
	df = pd.read_table(outfile, sep='\t', header=0)

	# filter region with count=0 in both groups.
	filter_df = df[(df.gA_mean_counts>0) & (df.gB_mean_counts>0)]

	# prepare data for boxplot

	ratio = filter_df.gB_mean_counts.sum() / filter_df.gA_mean_counts.sum()
	norm_v = np.log2(filter_df.gA_mean_counts / filter_df.gB_mean_counts * ratio)
	filter_df['Norm_V'] = norm_v
	pdata = defaultdict(list)

	for i in range(filter_df.shape[0]):
		sub_line = filter_df.iloc[i]
		TYPE = str(sub_line.Name)
		NV = sub_line.Norm_V
		pdata[TYPE].append(NV)

	# into plotting
	outfig = args.outfigure
	title = " "
	if args.title:
		title = args.title

	boxplot(pdata, title, outfig)


	# Done.


if __name__ == '__main__':
	main()
