#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''
Description:
Used to generate signal tracks with bigwig format that can be imported into IGV or other visual tools


Note:
	This can be a time-consuming operation due to large number of repeated sampling, therefore we do not recommend the construction of genome-wide signal tracks

'''

import sys, argparse, re
import pysam
import pyBigWig
import numpy as np

# get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-b', required=True, help='Bam file with .bai .')
parser.add_argument('-o', required=True, help='File name to save results.')
parser.add_argument('-r', type=str, help='A genome region with following format:\n\tchr1:10000-20000 .')

args = parser.parse_args()


# parse region
pattern = re.compile(r'(chr[0-9X]\d?):(\d*)-(\d*)')

m = pattern.search(args.r)
chrom = str(m.group(1))
start = int(m.group(2))
end = int(m.group(3))

# check if chrom in bamfile
bamfile = pysam.AlignmentFile(args.b,'rb')
try:
	assert chrom in bamfile.references
except AssertionError:
	sys.stderr.write('Error:\n{} do not contains {} records.'.format(str(args.b), chrom))
	sys.exit()


total = bamfile.mapped
max_size = bamfile.lengths[bamfile.references.index(chrom)]
vs = []


for i in range(start,end,10):
	sub_start = i - 150
	sub_end = i + 150
	counts = bamfile.count(chrom, sub_start, sub_end, read_callback='nofilter')

	norm_v = counts / total / 300 * 1e10
	vs = np.append(vs, norm_v)

#
bamfile.close()


out_chroms = [chrom] * len(vs)
out_starts = list(range(start,end,10))
out_ends = list(range(start+10,end+10,10))
vs = list(vs)

assert len(out_chroms) == len(vs) == len(out_starts) == len(out_ends)

#
bw = pyBigWig.open(args.o,'w')
bw.addHeader([(chrom, max_size)])
bw.addEntries(out_chroms, out_starts, ends=out_ends, values=vs)

bw.close()
