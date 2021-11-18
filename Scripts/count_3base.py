#!/usr/bin/env python
# -*- coding:utf-8 -*-
# usage : count_3base.py xxx.bed(.gz) xxx.3bases
import sys, gzip

def count_multi_bases(seq_list, i=3, j=6):
    base_dict = {}
    for seq in seq_list:
        base = seq[i:j]
        k = base_dict.get(base, 0)
        base_dict[base] = k + 1
    return base_dict

bed = sys.argv[1]
if not bed.endswith('gz'):
    seq_list = [line.split("\t")[3] for line in open(bed)]
else:
    seq_list = [line.decode().split("\t")[3] for line in gzip.open(bed)]


outfile = open(sys.argv[2], 'w')
base_dict = count_multi_bases(seq_list, 3, 6)

total = 0
for key in base_dict.keys():
    count = base_dict.get(key, 0)
    total += count
    print(key + "\t" + str(count),file=outfile)
    
print("Total\t" + str(total),file=outfile)

outfile.close()
