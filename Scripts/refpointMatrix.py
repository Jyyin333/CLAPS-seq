#!/usr/bin/env python
# -*- coding:utf-8 -*-

'''
Description:
	This script calculates scores along genome per region,
	Create a matrix that can be used to plot distribution
	Importantly, different with deeptools computeMatrix, \
	here we divide signal into two parts according to template-strand or non-template-strand

Note:
	The input score files should be in BAM format
	You need to use another script if using bigWig file while the result is the same

'''

import numpy as np
import pandas as pd
import pysam
import argparse, gzip, sys, os
import multiprocessing
from collections import defaultdict


def get_arguments():
	parser = argparse.ArgumentParser()

	required = parser.add_argument_group('Required arguments')
	required.add_argument('-b','--bam',required=True,help='Indexed Bam-file')
	required.add_argument('-r','--region',required=True,help=
		'Genome region file with strand-information in Bed-format.')

	optional = parser.add_argument_group('Optional arguments')
	optional.add_argument('-refpoint',type=str,default='tss',help='Options: tss, tes, center '
		'corresponding to the 2nd column, 3rd column in the region file or their average.')
	
	optional.add_argument('-a','--upstream',type=int,default=5000,help=
		'Distance upstream of the reference-point [bp].'
		'default:5000')

	optional.add_argument('-d','--downstream',type=int,default=5000,help=
		'Distance downstream of the reference-point [bp].'
		'default:5000')

	optional.add_argument('-bs','--binsize',type=int,default=50,
		help='Length in bases of the binsize' 
		'calculating scores of the genome region of binsize. [bp].'
		'default:50')

	optional.add_argument('--stepsize',type=int,default=-1,help=
		'Length in bases of the distance between bins'
		'in case of low sequencing depth, repeat sampling to amplify the signal.'
		'default is None.')

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


def process_refsite(frag_start,frag_end,frag_strand,rp='tss'):
	ref_site = 0
	if rp == 'tss':
		if frag_strand == '+':
			ref_site = frag_start
		if frag_strand == '-':
			ref_site = frag_end
	elif rp == 'tes':
		if frag_strand == '+':
			ref_site = frag_end
		if frag_strand == '-':
			ref_site = frag_start
	elif rp == 'center':
		ref_site = int(round((frag_start+frag_end)/2,0))

	return ref_site


def count_reads_in_bin_wrapper(args):
	return count_reads_in_bin(*args)

def count_reads_in_bin(chrom,refpoint,strand,offset=0):
	global global_vars

	ts,nts = [],[]
	refpoint = refpoint
	upstream = global_vars['upstream']
	downstream = global_vars['downstream']
	stepsize = global_vars['stepsize']
	binsize = global_vars['binsize']

	bf = pysam.AlignmentFile(global_vars['bam'], 'rb')
	total_mapped = bf.mapped

	if strand == '+':
		region_up = refpoint - upstream
		region_down = refpoint + downstream
		if region_up < 0 or region_down < 0:
			return False
		sites = np.arange(region_up,region_down,stepsize,dtype=int)
		for idx in range(len(sites)):
			istart = sites[idx]
			iend = istart + binsize
			ts.append(bamCountinRegion(bf,chrom,istart,iend,offset=offset,strand='-'))
			nts.append(bamCountinRegion(bf,chrom,istart,iend,offset=offset,strand='+'))

	elif strand == '-':
		region_up = refpoint + upstream
		region_down = refpoint - downstream
		if region_up < 0 or region_down < 0:
			return False
		sites = np.arange(region_up,region_down,-stepsize,dtype=int)
		for idx in range(len(sites)):
			iend = sites[idx]
			istart = iend - binsize
			ts.append(bamCountinRegion(bf,chrom,istart,iend,offset=offset,strand='+'))
			nts.append(bamCountinRegion(bf,chrom,istart,iend,offset=offset,strand='-'))

	
	bf.close()

	ts.extend(nts)
	merged = [round(x / total_mapped / binsize * 1e9, 6) for x in ts] # similar to RPKM in RNA-seq
	return merged


def main():

	# parse arguments
	args = get_arguments().parse_args()
	binsize = args.binsize
	upstream = args.upstream
	downstream = args.downstream
	stepsize = args.stepsize
	offset = args.offset
	cores = args.cores
	refpoint = args.refpoint

	if stepsize == -1:
		stepsize = binsize
	if stepsize > binsize:
		stepsize = binsize
		sys.stderr.write("stepsize can't be greater than binsize, set to same as binsize\n")
	if not refpoint in ['tss','tes','center']:
		refpoint = 'tss'
		sys.stderr.write('Invalid params:refpoint, set to tss.\n')

	# create matrix file and output the header line
	outFile = gzip.open(args.outMatrix,'wb')
	boundary = len(np.arange(0,upstream+downstream,stepsize))
	header = '# Mode:refpoint\tUpstream:{}\tDownstream:{}\tBinsize:{}\tBoundary:{}\tStepsize:{}\tRegionFileUsing:{}\tScoreFileUsing:{}\n'.format(
		args.upstream, args.downstream, args.binsize, boundary, stepsize, os.path.basename(args.region), os.path.basename(args.bam))
	outFile.write(header.encode('utf-8'))
	outFile.close()

	# create a dataframe to save values
	nbins = (int(abs(downstream+upstream)) // stepsize) * 2
	combined_data = pd.DataFrame(0,index=[],columns=range(nbins),dtype=float)

	# create global variables
	global global_vars
	global_vars = {}
	global_vars['bam'] = args.bam
	global_vars['upstream'] = upstream
	global_vars['downstream'] = downstream
	global_vars['stepsize'] = stepsize
	global_vars['binsize'] = binsize
	global_vars['outFile'] = outFile


	# create empty list to save sub_task
	TASKS = []
	line_counter = 0

	with open(args.region,'r') as fi:
		for line in fi:
			line = line.rstrip().split('\t')
			frag_chrom = line[0]
			frag_start = int(float(line[1]))
			frag_end = int(float(line[2]))
			frag_strand = line[5]
			frag_name = line[4]
			ref_site = process_refsite(frag_start,frag_end,frag_strand,rp=refpoint)
			arglist = [frag_chrom,ref_site,frag_strand,offset]
			TASKS.append(tuple(arglist))
			line_counter += 1

	pool = multiprocessing.Pool(cores)
	res = pool.map_async(count_reads_in_bin_wrapper,TASKS).get()
	pool.close()
	pool.join()

	invalid_regions = 0
	idx_counter = 0
	for item in res:
		if item:
			tmp = pd.DataFrame(np.array(item).reshape(1,nbins),index=[idx_counter],dtype=float)
			combined_data = combined_data.append(tmp)
			idx_counter += 1
		# remove regions out of valid range
		else:
			invalid_regions += 1
			continue

	# output scores matrix
	combined_data.to_csv(args.outMatrix,sep='\t',header=False,index=False,mode='ab',compression='gzip',encoding='utf-8')

	if invalid_regions != 0:
		sys.stderr.write('Warning: Of total {} regions, {} out of valid range.\n'.format(line_counter,invalid_regions))


if __name__ == '__main__':
	main()
