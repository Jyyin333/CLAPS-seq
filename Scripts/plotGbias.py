#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''

Description:
Plot G/C bias using Readcounts file generated from computeGbias.py

'''


import sys, os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats



def bin_by(x, y, nbins=10):

	bins = np.linspace(0, 1, nbins+1)
	bins[-1] += 1
	indices = np.digitize(y, bins)

	res = []
	for i in range(1, len(bins)):
		res.append(x[indices == i])

	bins = bins[:-1]

	return res, bins


def get_ci(array):
	# 
	# calculate standard error and mean
	SEM = stats.sem(array)
	MEAN = np.mean(array)

	# default using 95% confidence interval
	# corresponding Z statistic: 1.96
	up = MEAN + (1.96 * SEM)
	low = MEAN - (1.96 * SEM)

	return low, up


def _transfer(fname):

	res = []

	with open(fname) as fi:
		for line in fi:
			if line.startswith('#'):
				continue

			else:
				line = line.rstrip().split('\t')
				tc = int(float(line[4]))
				cc = int(float(line[5]))
				if tc == 0 or cc == 0:
					continue

				tp = int(float(line[0]))
				tm = int(float(line[1]))
				cp = int(float(line[2]))
				cm = int(float(line[3]))
				pg = float(line[6])
				cg = float(line[7])

				if tp != 0 and cp != 0:
					res.append([tp, cp, pg])

				if tm != 0 and cm != 0:
					res.append([tm,cm,cg])

				else:
					continue

	return res


def plot_G(fname, outfig, se):

	dataname = str(os.path.basename(fname))
	title = dataname.split('.')[0].split('_')[0]

	# transfer data type
	data = _transfer(fname)

	tall = sum([x[0] for x in data])
	call = sum([x[1] for x in data])
	ratio = call / tall

	reads_per_gc = []

	for item in data:
		n_count = item[0] / item[1] * ratio
		reads_per_gc.append([n_count, item[2]])

	reads_per_gc = [[np.log2(x[0]),x[1]] for x in reads_per_gc]
	reads_per_gc = np.asarray(reads_per_gc, dtype=float)

	# prepare data for  plotting
	reads, GC = reads_per_gc.T
	signal_per_gc, bin_labels = bin_by(reads, GC, nbins=100)
	# set GC limit
	to_keep = [idx for idx, x in enumerate(bin_labels) if 0.1 <= round(x,3) <= 0.35]
	signal_per_gc = [signal_per_gc[x] for x in to_keep]
	bin_labels = [bin_labels[x] for x in to_keep]

	y,low,up = [],[],[]
	for i in range(len(signal_per_gc)):
		ci = get_ci(signal_per_gc[i])
		low.append(ci[0])
		up.append(ci[1])
		y.append(np.mean(signal_per_gc[i]))


	x = range(len(y))
	fig, ax = plt.subplots(figsize=(9,6),dpi=300)
	ax.plot(x,y,ls='-',lw=1.5,color='#9ACD32')
	if se:
		ax.fill_between(x,low,up,facecolor='#9ACD32',alpha=0.5)

	ax.set_xlabel('G ratio',fontsize=12, weight='bold')
	ax.set_ylabel('Normalized OG signal',fontsize=12,fontdict={'size':15,'family':'DejaVu Serif','weight':600})
	xticks = list(range(0,len(y)+1,5))
	ax.set_xticks(xticks)
	ax.set_xticklabels([0.1,0.15,0.2,0.25,0.3,0.35])
	ax.set_title(title, fontsize=15, weight='bold')
	fmt = str(outfig).split('.')[-1]
	fig.savefig(outfig, dpi=300, bbox_inches='tight', format=fmt)
	plt.close()

	return


def plot_GC(fname, outfig, se=False):

	dataname = str(os.path.basename(fname))
	title = dataname.split('.')[0].split('_')[0]

	reads_per_gc = []

	tall = sum([float(x.rstrip().split('\t')[4]) for x in open(fname) if not x.startswith('#')])
	call = sum([float(x.rstrip().split('\t')[5]) for x in open(fname) if not x.startswith('#')])

	with open(fname) as fi:
		quote = fi.readline()

		for line in fi:
			line = line.rstrip().split('\t')
			t_count = int(float(line[4]))
			c_count = int(float(line[5]))
			gcratio = float(line[-1]) + float(line[-2])

			# filter count==0
			if c_count == 0 and t_count != 0:
				tall -= t_count
				continue

			if c_count != 0 and t_count == 0:
				call -= c_count
				continue

			if c_count == 0 and t_count == 0:
				continue

			n_count = float(t_count / c_count)
			reads_per_gc.append([n_count,gcratio])


	ratio = call / tall
	# normalize 2 inpput and take base2 logarithm
	reads_per_gc = [[x[0]*ratio, x[1]] for x in reads_per_gc]
	reads_per_gc = [[np.log2(x[0]),x[1]] for x in reads_per_gc]
	reads_per_gc = np.asarray(reads_per_gc, dtype=float)

	# prepare data for  plotting
	reads, GC = reads_per_gc.T
	signal_per_gc, bin_labels = bin_by(reads, GC, nbins=100)
	# set GC limit
	to_keep = [idx for idx, x in enumerate(bin_labels) if 0.2 <= x <= 0.7]
	signal_per_gc = [signal_per_gc[x] for x in to_keep]
	bin_labels = [bin_labels[x] for x in to_keep]

	y,low,up = [],[],[]
	for i in range(len(signal_per_gc)):
		ci = get_ci(signal_per_gc[i])
		low.append(ci[0])
		up.append(ci[1])
		y.append(np.mean(signal_per_gc[i]))


	x = range(len(y))
	fig,ax = plt.subplots(figsize=(9,6),dpi=300)
	ax.plot(x,y,ls='-',lw=1.5,color='#9ACD32')
	if se:
		ax.fill_between(x,low,up,facecolor='#9ACD32',alpha=0.5)

	ax.set_xlabel('GC ratio',fontsize=12, weight='bold')
	ax.set_ylabel('Normalized OG signal',fontsize=12,fontdict={'size':15,'family':'DejaVu Serif','weight':600})
	xticks = list(range(0,len(y)+1,10))
	ax.set_xticks(xticks)
	ax.set_xticklabels([0.2,0.3,0.4,0.5,0.6,0.7])
	ax.set_title(title, fontsize=15, weight='bold')
	fmt = str(outfig).split('.')[-1]
	fig.savefig(outfig, dpi=300, bbox_inches='tight', format=fmt)
	plt.close()

	return



def main():

	parser = argparse.ArgumentParser()

	parser.add_argument('-i', help='ReadCounts file generated from computeGbias.py')

	parser.add_argument('-o', help='File name to save figure.')

	parser.add_argument('-se', type=str, default='false', help=
		'Whether the final figure contains 95%% confidence interval. Default is False')

	parser.add_argument('-gc', type=str, default='true', help=
		'By default, plotting relationship between OG and GC content. '
		'If false, plotting relationship between OG and G content. Default is True')

	args = parser.parse_args()


	infile = args.i
	outfile = args.o

	# check -gc
	gc = str(args.gc).lower()
	if gc in ['true','t','y','yes']:
		gc = True
	else:
		gc = False


	# check -se
	se = str(args.se).lower()
	if se in ['true','t','y','yes']:
		se = True
	else:
		se = False

	# into plot
	if gc:
		plot_GC(infile,  outfile, se)
	else:
		plot_G(infile,  outfile, se)


if __name__ == '__main__':
	main()
