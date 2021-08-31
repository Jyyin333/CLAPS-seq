#!/usr/bin/env python
# -*- coding:utf-8 -*-
import sys,gzip

bed = sys.argv[1]
fo = open(sys.argv[2],'w')

if str(bed).endswith('.gz'):
    with gzip.open(bed) as f:
        for line in f:
            line = line.decode().rstrip().split('\t')
            sequence = str(line[3])
            if not sequence[4] == 'G':
                continue
            else:
                print('\t'.join(map(str,line)),file = fo)

else:
    with open(bed) as f:
        for line in f:
            line = line.rstrip().split('\t')
            sequence = str(line[3])
            if not sequence[4] == 'G':
                continue
            else:
                print('\t'.join(map(str,line)),file = fo)

fo.close()
