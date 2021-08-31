#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''
Description:
Used to generate background G content profile along gene bodies
Calculating background G content in scaled region.

Note:
	Poor compatibility and generalization, for this study only.

'''


import sys, os
import gzip
import py2bit
from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import multiprocessing
import matplotlib.pyplot as plt


def get_arguments():

	parser = argparse.ArgumentParser(description=
		'Calculating background G/C content in sclaed region.')

	parser.add_argument('-r','--region', required = True, help=
		'Genome region file with strand-information in Bed-format.')

	parser.add_argument('-g','--genome', required = True, help=
		'Genome in two bit format. [.2bit]')

	parser.add_argument('-a', '--upstream', type=int, default=2000, help='Distance upstream of the reference-point [bp].')
	
	parser.add_argument('-d', '--downstream', type=int, default=2000, help='Distance downstream of the reference-point [bp].')

	parser.add_argument('-m','--bodylength',type=int,default=5000,help=
		'Distance in bases to which all regions will be fit [bp].'
		'default:5000')

	parser.add_argument('-bs', '--binsize', type=int, default=50, help='Length in bases of the binsize')

	parser.add_argument('-o', '--outfigure', help='File name to save the image. '
		'The image format is depend on the suffix of input name.')

	parser.add_argument('-j', '--cores', required = False, type=int, default=1, help=
		'Number of processor.\n\tdefault:1')

	parser.add_argument('--outmatrix', required = False, help=
		'If true, raw g/c content score matrix would be saved into this file.')


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



def plot_scale(ts, nts, xticks, xlabels, outfig):

	fig, ax = plt.subplots(figsize=(9,6),dpi=300)
	x = range(len(ts))

	ax.plot(x, ts, ls='solid', lw=1.2, color='blue', label='TS')
	ax.plot(x, nts, ls='solid', lw=1.2, color='red', label='NTS')
	ax.set_xlim(-1,len(x))
	ax.set_xticks(xticks)
	ax.set_xticklabels(xlabels,rotation=0,fontdict={'size':12,'family':'Verdana','weight':500})
	ax.set_ylabel('Background G content',fontdict={'size':15,'family':'DejaVu Serif','weight':600})
	ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5))

	fig.savefig(outfig, dpi=300, bbox_inches="tight")
	plt.close()

	return



def scale_wrapper(args):
	return scale_wrapper_worker(*args)


def scale_wrapper_worker(chrom,start,end,strand):

	# load global vars
	binsize = global_vars['binsize']
	upstream = global_vars['upstream']
	downstream = global_vars['downstream']
	bodylength = global_vars['bodylength']
	tb = py2bit.open(global_vars['tbit'])

	up_bins, mid_bins, down_bins = global_vars['bins']
	nbins = int((up_bins + mid_bins + down_bins)*2)
	ts, nts = [], []

	# for positive strand
	if strand == '+':
		up_region = [start - upstream, start]
		mid_region = [start, end]
		down_region = [end, end + downstream]

		if up_region[0] < 0:
			return None

		# upstream-region
		up_pos = np.linspace(up_region[0], up_region[1], up_bins, endpoint=False, dtype=int)
		up_pos = np.append(up_pos,up_region[1])
		for i in range(len(up_pos)-1):
			sub_start = up_pos[i]
			sub_end = up_pos[i+1]
			nts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True))
			ts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True,base='C'))

		# mid-region
		mid_pos = np.linspace(mid_region[0], mid_region[1], mid_bins, endpoint=False, dtype=int)
		mid_pos = np.append(mid_pos,mid_region[1])
		for i in range(len(mid_pos)-1):
			sub_start = mid_pos[i]
			sub_end = mid_pos[i+1]
			nts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True))
			ts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True,base='C'))

		# down-region
		down_pos = np.linspace(down_region[0], down_region[1], down_bins, endpoint=False, dtype=int)
		down_pos = np.append(down_pos,down_region[1])
		for i in range(len(down_pos)-1):
			sub_start = down_pos[i]
			sub_end = down_pos[i+1]
			nts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True))
			ts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True,base='C'))

		ts.extend(nts)
		merged = ts[:]
		assert len(merged) == nbins


	# for negtive strand
	elif strand == '-':
		up_region = [end, end + upstream]
		mid_region = [start, end]
		down_region = [start - downstream, start]

		if down_region[0] < 0:
			return None

		# down-region
		down_pos = np.linspace(down_region[0], down_region[1], down_bins, endpoint=False, dtype=int)
		down_pos = np.append(down_pos,down_region[1])
		for i in range(len(down_pos)-1):
			sub_start = down_pos[i]
			sub_end = down_pos[i+1]
			nts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True,base='C'))
			ts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True))

		# mid-region
		mid_pos = np.linspace(mid_region[0], mid_region[1], mid_bins, endpoint=False, dtype=int)
		mid_pos = np.append(mid_pos,mid_region[1])
		for i in range(len(mid_pos)-1):
			sub_start = mid_pos[i]
			sub_end = mid_pos[i+1]
			nts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True,base='C'))
			ts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True))

		# upstream-region
		up_pos = np.linspace(up_region[0], up_region[1], up_bins, endpoint=False, dtype=int)
		up_pos = np.append(up_pos,up_region[1])
		for i in range(len(up_pos)-1):
			sub_start = up_pos[i]
			sub_end = up_pos[i+1]
			nts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True,base='C'))
			ts.append(get_base_content(chrom,sub_start,sub_end,tb,fraction=True))
	
		ts = ts[::-1]
		nts = nts[::-1]

		ts.extend(nts)
		merged = ts[:]
		assert len(merged) == nbins
		
	tb.close()
	return merged



def main():

	# parse arguments
	args = get_arguments().parse_args()
	cores = args.cores

	up_bins = args.upstream // args.binsize
	mid_bins = args.bodylength // args.binsize
	down_bins = args.downstream // args.binsize
	nbins = int((up_bins + mid_bins + down_bins)*2)

	global global_vars
	global_vars = {}
	global_vars['tbit'] = args.genome
	global_vars['upstream'] = args.upstream
	global_vars['downstream'] = args.downstream
	global_vars['bodylength'] = args.bodylength
	global_vars['binsize'] = args.binsize
	global_vars['bins'] = [up_bins, mid_bins, down_bins]

	# adding header to prepared outfile
	if args.outmatrix:
		out_Mat = args.outmatrix
		outfile = gzip.open(out_Mat,'wb')
		header = '# RegionFileUsing:{}\tUpstream:{}\tDownstream:{}\tBodylength:{}\tBinsize:{}\tBoundary:{}\n'.format(
			os.path.basename(args.region), args.upstream, args.downstream, args.bodylength, args.binsize, int(up_bins + mid_bins + down_bins))
		outfile.write(header.encode('utf-8'))
		outfile.close()

	# create an empty dataframe
	DF = pd.DataFrame(0,index=[],columns=range(nbins),dtype=float)

	# create mutilthread pools and run
	TASKS = []
	line_counter = 0
	with open(args.region) as f:
		for line in f:
			line_counter += 1
			line = line.rstrip().split('\t')
			frag_chrom = line[0]
			frag_start = int(float(line[1]))
			frag_end = int(float(line[2]))
			frag_strand = line[5]
			arglist = [frag_chrom,frag_start,frag_end,frag_strand]
			TASKS.append(tuple(arglist))

	pool = multiprocessing.Pool(cores)
	res = pool.map_async(scale_wrapper,TASKS).get()
	pool.close()
	pool.join()


	idx_counter = 0
	skip_counter = 0
	for item in res:
		if item:
			tmp_df = pd.DataFrame(np.array(item).reshape(1,nbins),index=[idx_counter],dtype=float)
			DF = DF.append(tmp_df)
			idx_counter += 1
		else:
			skip_counter += 1
			continue

	# save results
	if 'out_Mat' in locals().keys():
		DF.to_csv(out_Mat,sep='\t',header=False,index=False,mode='ab',compression='gzip',encoding='utf-8')

	if skip_counter:
		sys.stderr.write('Warning: Of total {} regions, {} out of range.\n'.format(line_counter,skip_counter))

	# plot profile
	# prepare data for plotting
	TS_res, NTS_res = [], []
	cut = int(up_bins + mid_bins + down_bins)
	ts_df = DF.iloc[:, :cut]
	nts_df = DF.iloc[:, cut:]
	ts_mean = ts_df.apply(lambda x : x.mean(), axis=0).tolist()
	nts_mean = nts_df.apply(lambda x : x.mean(), axis=0).tolist()

	start_label = '+' + str(args.upstream / 1000) + 'kb'
	end_label = '-' + str(args.downstream / 1000) + 'kb'
	up_tick_loc = args.upstream / args.binsize
	down_tick_loc = args.downstream / args.binsize
	x_ticks = [0, up_tick_loc-1, cut-down_tick_loc-1, cut-1]
	xlabels = [start_label,'start','end',end_label]

	plot_scale(ts_mean, nts_mean, x_ticks, xlabels, args.outfigure)



if __name__ == '__main__':
	main()
