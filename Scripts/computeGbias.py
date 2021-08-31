#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''

Description:
Be used to explore the relationship between OG distribution and background G/C content.
This script is a modified & simplified version of computeGCbias function in deeptools to get an appropriate output for downstream analysis.
You can view computeGCbias codes at https://github.com/deeptools/deepTools/tree/master/deeptools

'''

import sys, os
import argparse
import random
import numpy as np
import pysam
import multiprocessing
from pyfaidx import Fasta


def get_arguments():

	parser = argparse.ArgumentParser()

	parser.add_argument('-b','--bam', required = True, nargs = 2, help=
		'Bamfile with .bai index.  [ Treat.bam  Control.bam ]')

	parser.add_argument('-g','--genome', required=True, help=
		'Genome in fasta format.  e.g. hg19.fa')

	parser.add_argument('-s','--samplesize', type=int, default=10000000, help=
		'Number of sampling points to be considered. default:10,000,000')

	parser.add_argument('-l','--regionsize', type=int, default=1000, help=
		'To plot the reads per %%GC over a regionthe size of the region is required. default:1000')

	parser.add_argument('-o','--out', help='file name to save results.')

	parser.add_argument('-j','--cores', type=int, default=1, help=
		'Number of processors to use. default:1')

	return parser



def _cal_effective_genome_size(genome, keep_chrom = None):
	'''
	return:
			dict = {
				'total_effective_genomesize':10000,
				'total_genomesize':12000,
				'chr1' :{
					'A':1, 'C':1,'G':1,'T':1,'effective':1
				},
			}
	'''
	genome = Fasta(genome)
	chrom_list = genome.keys()

	if keep_chrom:
		chrom_list = keep_chrom
		for i in chrom_list:
			assert i in genome.keys()

	res = {}
	gsize = 0
	for chrom in chrom_list:
		if not chrom in res.keys():
			res[chrom] = {}

		max_size = len(genome[chrom])
		gsize += max_size

		for base in ['A','C','G','T','N']:
			res[chrom][base] = genome[chrom][:max_size].seq.upper().count(base)

		res[chrom]['effective'] = res[chrom]['A'] + res[chrom]['C'] + res[chrom]['G'] + res[chrom]['T']

	esize = 0
	for chrom in chrom_list:
		esize += res[chrom]['effective']

	res['total_effective_genomesize'] = esize
	res['total_genomesize'] = gsize

	genome.close()

	return res



def _get_base_content_worker(*vargs):

	chrom, a, b, ref = vargs
	res = {}
	fres = {}

	bases_seq = ref[chrom][a:b].seq.upper()
	for Base in ['A','C','G','T','N']:
		bc = bases_seq.count(Base)
		res[Base] = bc

	tc = sum(res.values())
	for k,v in res.items():
		nv = float(v / tc)
		fres[k] = nv
	
	return res



def get_base_content(chrom, start, end, fa, base='G', fraction=True):

	# chech if base valid
	base = base.upper().strip()
	for i in range(len(base)):
		try:
			assert base[i] in ['A','C','G','T']
		except AssertionError:
			sys.stderr.write('Invalid base type.\n')
			return None

	if end > len(fa[chrom]):
		end = len(fa[chrom])

	bases = _get_base_content_worker(chrom, start, end, fa)

	if len(base) == 1:
		if fraction:
			return (bases[base]) / float(end - start)
		return bases[base]

	else:
		if fraction:
			return sum([bases[base[i]] for i in range(len(base))]) / float(end - start)
		return sum([bases[base[i]] for i in range(len(base))])



def _getTotalReads(f):

	bf = pysam.AlignmentFile(f,'rb')
	mapped = bf.mapped
	bf.close()

	return mapped


def fast_merge(*ARGS):

	sub_list = ARGS[0]
	tmp = []
	for i in sub_list:
		tmp += i
	return tmp


def multiProc(staticArgs, func, chromSize, genomeChunkSize=None, numberofProcessors=1):

	if not genomeChunkSize:
		genomeChunkSize = 1e5
		genomeChunkSize = int(genomeChunkSize)

	TASKS = []

	for chrom, size in chromSize:
		start = 0
		for startPos in range(start, size, genomeChunkSize):
			endPos = min(startPos + genomeChunkSize, size)

			regions = [[startPos, endPos]]
			for reg in regions:
				argsList = []

				argsList.extend([chrom, reg[0], reg[1]])
				argsList.extend(staticArgs)
				# argsList : (chr1, 0, 0 + chunksize, stepsize, regionszie)

				TASKS.append(tuple(argsList))


	if len(TASKS) > 1 and numberofProcessors > 1:
		random.shuffle(TASKS)
		pool = multiprocessing.Pool(numberofProcessors)
		res = pool.map_async(func, TASKS).get()
		pool.close()
		pool.join()

	else:
		res = list(map(func, TASKS))

	return res


def _get_chromsize(bam):

	bf = pysam.AlignmentFile(bam)
	chromsize = [(x[0],x[1]) for x in zip(bf.references,bf.lengths)]
	bf.close()

	return chromsize


def cal_signal_per_GC_wrapper(args):
	return cal_signal_per_GC_worker(*args)


def cal_signal_per_GC_worker(chrom,start,end,stepSize,regionSize):

	ref_fa = global_vars['ref_fa']
	t_bf = pysam.AlignmentFile(global_vars['T'],'rb')
	c_bf = pysam.AlignmentFile(global_vars['C'],'rb')
	sub_signal_per_gc = []

	positions = np.arange(start,end,stepSize)

	for idx in range(len(positions)):
		i = positions[idx]
		if len(ref_fa[chrom]) < i + regionSize:
			break

		try:
			g_content = get_base_content(chrom=chrom, start=int(i), end=int(i + regionSize),
				fa=ref_fa, base='G', fraction=True)
			c_content = get_base_content(chrom=chrom, start=int(i), end=int(i + regionSize),
				fa=ref_fa, base='C', fraction=True)

		except Exception as detail:
			sys.stderr.write("{}:{}-{}\n".format(chrom, i, i + regionSize))
			sys.stderr.write(detail)
			continue

		tp, tm, cp, cm = 0, 0, 0, 0
		for read in t_bf.fetch(chrom, int(i), int(i+regionSize)):
			if read.is_reverse:
				tm += 1
			else:
				tp += 1

		for read in c_bf.fetch(chrom, int(i), int(i+regionSize)):
			if read.is_reverse:
				cm += 1
			else:
				cp += 1


		#tp, tm, cp, cm = _count_reads_by_strand(t_bf, c_bf, chrom, int(i), int(i+regionSize))
		tall = tp + tm
		call = cp + cm

		sub_signal_per_gc.append(tuple([tp, tm, cp, cm, tall, call, g_content, c_content]))

	t_bf.close()
	c_bf.close()
	return sub_signal_per_gc




def main():

	args = get_arguments().parse_args()

	t_bam, c_bam = args.bam
	genome = args.genome
	samplesize = args.samplesize
	regionsize = args.regionsize
	cores = args.cores

	# split genome into chunks
	fa_dict = _cal_effective_genome_size(genome)
	esize = fa_dict.get('total_effective_genomesize')
	gsize = fa_dict.get('total_genomesize')
	t_total_mapped = _getTotalReads(t_bam)

	stepsize = max(int(gsize / samplesize), 1)
	old = _get_chromsize(t_bam)
	keeps = ['chr'+str(i) for i in range(1,23)] + ['chrX']
	chromsize = [x for x in old if x[0] in keeps]
	reads_per_bp = float(t_total_mapped / esize)
	chunksize = int(min(2e6, 4e5 / reads_per_bp))

	#
	global global_vars
	global_vars = {}
	global_vars['T'] = t_bam
	global_vars['C'] = c_bam
	global_vars['ref_fa'] = Fasta(genome)

	
	# run cal_signal_per_GC_worker
	res = multiProc((stepsize, regionsize),
		func = cal_signal_per_GC_wrapper,
		chromSize = chromsize,
		genomeChunkSize = chunksize,
		numberofProcessors = cores)


	# accerlerate merge
	stasks = []
	nlen = len(res)
	div = 20 if len(res) < 1000 else 100
	sp = nlen // div

	for i in range(sp):
		ys = tuple(res[i*div:(i+1)*div])
		stasks.append(ys)
		if nlen % div and i+1 == sp:
			stasks.append(tuple(res[(i+1)*div:]))

	pool = multiprocessing.Pool(cores)
	st = pool.map_async(fast_merge,stasks).get()
	pool.close()
	pool.join()

	signal_per_gc = []
	for i in st:
		signal_per_gc += i

	signal_per_gc = np.asarray(signal_per_gc)


	# save results
	np.savetxt(args.out, signal_per_gc, delimiter='\t', fmt='%s',
		header='T_plus  T_minus C_plus  C_minus  T_all  C_all  G_ratio  C_ratio')


	# close files
	global_vars['ref_fa'].close()
	


if __name__ == '__main__':
	main()
