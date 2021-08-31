#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''
Description:

Used to generate background G/C content profiles
Calculating background G/C content around reference point

You can set --gc paramater depending on the type of region to decide which Base content you want to plot. 
Also, you can specify a file to save raw values and replot using the results.

'''

import sys,argparse
import py2bit
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


def get_arguments():

	parser = argparse.ArgumentParser(description=
		'Calculating background G/C content around reference point.')

	parser.add_argument('-r','--region', required = True, help=
		'Genome region file with strand-information in Bed-format.')

	parser.add_argument('-g','--genome', required = True, help=
		'Genome in two bit format. [.2bit]')

	parser.add_argument('-refpoint', type=str, default='tss',help='Options: tss, tes, center. '
		'if not consider strand, this parameter corresponding to the 2nd column, 3rd column in the region file or their average.')

	parser.add_argument('-a', '--upstream', type=int, default=5000, help='Distance upstream of the reference-point [bp].')
	
	parser.add_argument('-d', '--downstream', type=int, default=5000, help='Distance downstream of the reference-point [bp].')

	parser.add_argument('-bs', '--binsize', type=int, default=50, help='Length in bases of the binsize')

	parser.add_argument('-o', '--outfigure', help='File name to save the image. '
		'The image format is depend on the suffix of input name.')

	parser.add_argument('--outmatrix', required = False, help=
		'If true, raw g/c content score matrix would be saved into this file.')

	parser.add_argument('--gc', type = bool, required = False, help=
		'By default, strand-split is not considered, that is to say, GC content over the genome will be calculated. '
		'Set this parameter to false to compute G content only. Note that the 6th column in REGION file requires strand information.')


	return parser





def run_g(args):

	# reload parameters
	RegionFile = args.region
	GenomeFile = args.genome
	binsize = args.binsize
	upstream = args.upstream
	downstream = args.downstream
	refpoint = args.refpoint.lower()

	# create DataFrame to save values
	nbins = (upstream + downstream) // binsize
	ts_df = pd.DataFrame(0,index=[],columns=range(nbins),dtype='float64')
	nts_df = pd.DataFrame(0,index=[],columns=range(nbins),dtype='float64')


	global global_vars
	global_vars = {}
	global_vars['2bit'] = GenomeFile

	#
	counter = 0
	with open(RegionFile) as fi:
		for line in fi:
			line = line.rstrip().split('\t')
			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			strand = line[5]
			up,down = g_process_refsite(start,end,upstream,downstream,strand,refpoint)

			if up < 0 or down < 0:
				continue

			ts, nts = count_G_by_bin(chrom,up,down,binsize,strand)
			if ts == 0 or nts == 0:
				continue

			sub_ts_df = pd.DataFrame(np.array(ts).reshape(1,nbins),index=[counter])
			ts_df = ts_df.append(sub_ts_df)
			sub_nts_df = pd.DataFrame(np.array(nts).reshape(1,nbins),index=[counter])
			nts_df = nts_df.append(sub_nts_df)

			counter += 1


	# calculating score completed

	# output results
	if args.outmatrix:
		combin_df = pd.concat([ts_df, nts_df], axis=1)
		combin_df.to_csv(args.outmatrix,sep='\t',header=False,index=False,mode='wb',compression='gzip',encoding='utf-8')


	# plot profile
	ts_mean = np.asarray(ts_df.apply(lambda x : x.mean(), axis=0),dtype=float)
	nts_mean = np.asarray(nts_df.apply(lambda x : x.mean(), axis=0),dtype=float)
	start_label = '+' + str(upstream / 1000) + 'kb'
	end_label = '-' + str(downstream / 1000) + 'kb'
	refsite = upstream//binsize
	xticks = [0,refsite,len(ts_mean)]
	center_label = 'CENTER'
	xlabels = [start_label,center_label,end_label]

	plot_g(args.outfigure, ts_mean, nts_mean, xticks, xlabels)




def plot_g(outfig, a, b, xticks, xlabels):

	x = list(range(len(a)))
	fig, ax = plt.subplots(figsize=(9, 6),dpi=300)
	ax.plot(x, a, ls='solid', lw=1.2, color='blue', label='TS')
	ax.plot(x, b, ls='solid', lw=1.2, color='red', label='NTS')
	#ax.set_ylim(0,1)
	#ax.set_yticks([0,.2,.4,.6,.8,1])
	#ax.set_yticklabels([0,.2,.4,.6,.8,1])
	ax.set_xlim(-1,len(x))
	ax.set_xticks(xticks)
	ax.set_xticklabels(xlabels,rotation=0,fontdict={'size':12,'family':'Verdana','weight':500})
	ax.set_ylabel('Background G content',fontdict={'size':15,'family':'DejaVu Serif','weight':600})
	ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5))

	fig.savefig(outfig,dpi=300,bbox_inches="tight")
	plt.close()

	return




def count_G_by_bin(chrom,up,down,binsize,strand):

	chrom = chrom
	binsize = binsize
	tb = py2bit.open(global_vars['2bit'])
	chrom_max_value = int(tb.chroms(chrom))
	if down > chrom_max_value:
		return 0, 0

	p,m = [],[]
	start_sites = list(range(up,down,binsize))
	for pos in start_sites:
		p.append(get_base_content(chrom,pos,pos+binsize,tb,fraction=True,base='G'))
		m.append(get_base_content(chrom,pos,pos+binsize,tb,fraction=True,base='C'))

	if None in p or None in m:
		return 0, 0

	if strand == '+':
		nts = p
		ts = m

	if strand == '-':
		nts = m[::-1]
		ts = p[::-1]

	tb.close()


	return ts,nts




def g_process_refsite(start,end,upstream,downstream,strand,refpoint='tss'):

	if refpoint == 'tss':
		if strand == '+':
			a = start - upstream
			b = start + downstream

		if strand == '-':
			b = end + upstream
			a = end - downstream

	elif refpoint == 'tes':
		if strand == '+':
			a = end - upstream
			b = end + downstream

		if strand == '-':
			b = start + upstream
			a = start - downstream

	else:
		c = int(round((start+end)/2,0))
		a = c - upstream
		b = c + downstream

	return a, b




def gc_process_refsite(start,end,upstream,downstream,rp='tss'):

	if rp == 'tss':
		a = start - upstream
		b = start + downstream

	elif rp == 'tes':
		a = end - upstream
		b = end + downstream

	else:
		midpoint = int(round((start+end)/2,0))
		a = midpoint - upstream
		b = midpoint + downstream

	return a, b


def count_GC_by_bin(chrom,up,down,binsize):
	chrom = chrom
	binsize = binsize
	tb = py2bit.open(global_vars['2bit'])
	chrom_max_value = int(tb.chroms(chrom))
	if down > chrom_max_value:
		return None

	gc = []
	start_sites = list(range(up,down,binsize))
	for pos in start_sites:
		gc.append(get_base_content(chrom,pos,pos+binsize,tb,fraction=True,base='GC'))
	tb.close()

	if None in gc:
		return None
	return gc


def run_gc(args):
	
	# reload parameters
	RegionFile = args.region
	GenomeFile = args.genome
	binsize = args.binsize
	upstream = args.upstream
	downstream = args.downstream
	refpoint = args.refpoint.lower()

	# create DataFrame to save values
	nbins = (upstream + downstream) // binsize
	combin_df = pd.DataFrame(0,index=[],columns=range(nbins),dtype='float64')


	global global_vars
	global_vars = {}
	global_vars['2bit'] = GenomeFile

	counter = 0
	with open(RegionFile) as fi:
		for line in fi:
			line = line.rstrip().split('\t')
			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			up,down = gc_process_refsite(start,end,upstream,downstream,refpoint)

			if up < 0 or down < 0:
				continue
				
			subGC = count_GC_by_bin(chrom,up,down,binsize)
			if subGC == None:
				continue

			subGC = np.asarray(subGC,dtype=float)
			sub_GC_df = pd.DataFrame(subGC.reshape(1,nbins),index=[counter])
			combin_df = combin_df.append(sub_GC_df)

			counter += 1


	# calculating score completed

	if args.outmatrix:
		combin_df.to_csv(args.outmatrix,sep='\t',header=False,index=False,mode='wb',compression='gzip',encoding='utf-8')


	# plot profile
	combin_mean = np.asarray(combin_df.apply(lambda x : x.mean(), axis=0),dtype=float)
	start_label = '+' + str(upstream / 1000) + 'kb'
	end_label = '-' + str(downstream / 1000) + 'kb'
	refsite = upstream//binsize
	xticks = [0,refsite,len(combin_mean)]
	center_label = 'CENTER'
	xlabels = [start_label,center_label,end_label]

	plot_gc(args.outfigure,combin_mean,xticks,xlabels)


def plot_gc(outfig, value, xticks, xlabels):

	x = list(range(len(value)))
	fig, ax = plt.subplots(figsize=(9, 6),dpi=300)
	ax.plot(x, value, ls='solid', lw=1.2, color='#32CD32')
	#ax.set_ylim(0,1)
	#ax.set_yticks([0,.2,.4,.6,.8,1])
	#ax.set_yticklabels([0,.2,.4,.6,.8,1])
	ax.set_xlim(-1,len(x))
	ax.set_xticks(xticks)
	ax.set_xticklabels(xlabels,rotation=0,fontdict={'size':12,'family':'Verdana','weight':500})
	ax.set_ylabel('Background GC content',fontdict={'size':15,'family':'DejaVu Serif','weight':600})

	fig.savefig(outfig,dpi=300,bbox_inches="tight")
	plt.close()

	return


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



def _check_if_strand(infile):

	with open(infile) as fi:
		line = fi.readline().strip().split('\t')
		try:
			strand = line[5]
		except IndexError:
			sys.stderr.write('Error IndexError:\nplease follow BED format standard, '
				'strand info should be put on 6th column.\n')
			sys.exit()
		except:
			sys.stderr.write('Error:\nplease follow BED format standard.\n')
			sys.exit()

		if not strand in ['+','-']:
			sys.stderr.write('Error:\nstrand info should be put on 6th column. [ +, - ]\n')
			sys.exit()

	return



def main():

	args = get_arguments().parse_args()

	binsize = args.binsize
	upstream = args.upstream
	downstream = args.downstream

	# check input parameters
	if not upstream % binsize == 0:
		raise ValueError('Error ValueError:\nupstream can not be divisible by binsize.\n')
		sys.exit()

	if not downstream % binsize == 0:
		raise ValueError('Error ValueError:\ndownstream can not be divisible by binsize.\n')
		sys.exit()

	if args.gc:
		# gc = False
		# check if region file contains strand infomation
		_check_if_strand(args.region)
		run_g(args)

	else:
		# gc = True
		run_gc(args)


	# All done.
	return


if __name__ == '__main__':
	main()
