#!/usr/bin/env python
# -*- coding:utf-8 -*-

'''
Description:
	Rather than center on one specific point, this script will shrunk or strech total region into samle length.
	Also signal was divided into two parts according to template-strand or non-template-strand

'''


import sys, os, gzip
from collections import defaultdict
import pandas as pd
import numpy as np
import argparse
import pysam
import multiprocessing



def get_arguments():
	
	parser = argparse.ArgumentParser()

	required = parser.add_argument_group('Required arguments')
	required.add_argument('-b','--bam',required=True,help='Indexed Bam-file')
	required.add_argument('-r','--region',required=True,help=
		'Genome region file with strand-information in Bed-format.')

	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('-a','--upstream',type=int,default=2000,help=
		'Distance upstream of the reference-point [bp].'
		'default:2000')

	optional.add_argument('-d','--downstream',type=int,default=2000,help=
		'Distance downstream of the reference-point [bp].'
		'default:2000')

	optional.add_argument('-m','--bodylength',type=int,default=5000,help=
		'Distance in bases to which all regions will be fit [bp].'
		'default:5000')

	optional.add_argument('-bs','--binsize',type=int,default=50,
		help='Length in bases of the binsize' 
		'calculating scores of the genome region of binsize. [bp].'
		'default:50')
	
	optional.add_argument('--offset',type=int,default=0,help=
		'Uses this offset inside of each read as the signal.'
		'default:0')

	optional.add_argument('-j','--cores',type=int,default=1,help='number of processor.\n\tdefault:1')

	outopt = parser.add_argument_group('Output arguments')
	outopt.add_argument('-o','--outMatrix',help='File name to save the matrix. [Ends with .gz]')

	return parser


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


def scale_wrapper(args):
	return scale_wrapper_worker(*args)


def scale_wrapper_worker(chrom,start,end,strand):

	# load global vars
	bf = pysam.AlignmentFile(global_vars['bam'],'rb')
	binsize = global_vars['binsize']
	upstream = global_vars['upstream']
	downstream = global_vars['downstream']
	bodylength = global_vars['bodylength']
	offset = global_vars['offset']
	total_mapped = bf.mapped
	chrom_max_size = bf.lengths[bf.references.index(chrom)]

	up_bins, mid_bins, down_bins = global_vars['bins']
	nbins = int((up_bins + mid_bins + down_bins)*2)
	ts, nts = [], []

	# for positive strand
	if strand == '+':
		up_region = [start - upstream, start]
		mid_region = [start, end]
		down_region = [end, end + downstream]

		# skip regions out of range
		if up_region[0] < 0 or down_region[1] > chrom_max_size:
			return None

		# upstream-region
		up_pos = np.linspace(up_region[0], up_region[1], up_bins, endpoint=False, dtype=int)
		up_pos = np.append(up_pos,up_region[1])
		for i in range(len(up_pos)-1):
			sub_start = up_pos[i]
			sub_end = up_pos[i+1]
			nts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='+') / total_mapped / (sub_end - sub_start) * 1e9 )
			ts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='-') / total_mapped / (sub_end - sub_start) * 1e9 )

		# mid-region
		mid_pos = np.linspace(mid_region[0], mid_region[1], mid_bins, endpoint=False, dtype=int)
		mid_pos = np.append(mid_pos,mid_region[1])
		for i in range(len(mid_pos)-1):
			sub_start = mid_pos[i]
			sub_end = mid_pos[i+1]
			nts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='+') / total_mapped / (sub_end - sub_start) * 1e9 )
			ts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='-') / total_mapped / (sub_end - sub_start) * 1e9 )

		# down-region
		down_pos = np.linspace(down_region[0], down_region[1], down_bins, endpoint=False, dtype=int)
		down_pos = np.append(down_pos,down_region[1])
		for i in range(len(down_pos)-1):
			sub_start = down_pos[i]
			sub_end = down_pos[i+1]
			nts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='+') / total_mapped / (sub_end - sub_start) * 1e9 )
			ts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='-') / total_mapped / (sub_end - sub_start) * 1e9 )

		ts.extend(nts)
		merged = ts[:]
		assert len(merged) == nbins

	# for negtive strand
	elif strand == '-':
		up_region = [end, end + upstream]
		mid_region = [start, end]
		down_region = [start - downstream, start]

		if down_region[0] < 0 or up_region[1] > chrom_max_size:
			return None

		# down-region
		down_pos = np.linspace(down_region[0], down_region[1], down_bins, endpoint=False, dtype=int)
		down_pos = np.append(down_pos,down_region[1])
		for i in range(len(down_pos)-1):
			sub_start = down_pos[i]
			sub_end = down_pos[i+1]
			nts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='-') / total_mapped / (sub_end - sub_start) * 1e9 )
			ts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='+') / total_mapped / (sub_end - sub_start) * 1e9 )

		# mid-region
		mid_pos = np.linspace(mid_region[0], mid_region[1], mid_bins, endpoint=False, dtype=int)
		mid_pos = np.append(mid_pos,mid_region[1])
		for i in range(len(mid_pos)-1):
			sub_start = mid_pos[i]
			sub_end = mid_pos[i+1]
			nts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='-') / total_mapped / (sub_end - sub_start) * 1e9 )
			ts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='+') / total_mapped / (sub_end - sub_start) * 1e9 )

		# upstream-region
		up_pos = np.linspace(up_region[0], up_region[1], up_bins, endpoint=False, dtype=int)
		up_pos = np.append(up_pos,up_region[1])
		for i in range(len(up_pos)-1):
			sub_start = up_pos[i]
			sub_end = up_pos[i+1]
			nts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='-') / total_mapped / (sub_end - sub_start) * 1e9 )
			ts.append(bamCountinRegion(bf,chrom,sub_start,sub_end,offset=offset,strand='+') / total_mapped / (sub_end - sub_start) * 1e9 )

		# reverse values outside the loops	
		ts = ts[::-1]
		nts = nts[::-1]

		ts.extend(nts)
		merged = ts[:]
		

	bf.close()
	return merged


def main():

	# parse arguments
	args = get_arguments().parse_args()
	cores = args.cores

	# calculate bins number
	up_bins = args.upstream // args.binsize
	mid_bins = args.bodylength // args.binsize
	down_bins = args.downstream // args.binsize
	nbins = int((up_bins + mid_bins + down_bins) * 2)
	boundary = int(up_bins + mid_bins + down_bins)

	# create global variables
	global global_vars
	global_vars = {}
	global_vars['bam'] = args.bam
	global_vars['upstream'] = args.upstream
	global_vars['downstream'] = args.downstream
	global_vars['bodylength'] = args.bodylength
	global_vars['binsize'] = args.binsize
	global_vars['offset'] = args.offset
	global_vars['bins'] = [up_bins, mid_bins, down_bins]

	# create matrix file and print the header line
	out_Mat = args.outMatrix
	outfile = gzip.open(out_Mat,'wb')
	header = '# Mode:scale\tUpstream:{}\tDownstream:{}\tBinsize:{}\tBoundary:{}\tBodylength:{}\tRegionFileUsing:{}\tScoreFileUsing:{}\n'.format(
		args.upstream, args.downstream, args.binsize, boundary, args.bodylength, os.path.basename(args.region), os.path.basename(args.bam))
	outfile.write(header.encode('utf-8'))
	outfile.close()

	# create a dataframe to save values
	DF = pd.DataFrame(0,index=[],columns=range(nbins),dtype=float)

	# create empty list to save sub_task
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
			frag_name = line[4]
			arglist = [frag_chrom,frag_start,frag_end,frag_strand]
			TASKS.append(tuple(arglist))

	# run
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

	DF.to_csv(out_Mat,sep='\t',header=False,index=False,mode='ab',compression='gzip',encoding='utf-8')

	if skip_counter:
		sys.stderr.write('Warning: Of total {} regions, {} out of range.\n'.format(line_counter,skip_counter))


if __name__ == '__main__':
	main()
