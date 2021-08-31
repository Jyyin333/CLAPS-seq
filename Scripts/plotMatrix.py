#!/usr/bin/env python
# -*- coding:utf-8 -*-


'''
Description:
Be used to plot stranded signal distribution with matrix[.gz] file generated from refpointMatrix.py or scaleMatrix.py
Default output image only contains one subplot
If control sample was provided, then image would have two subplots.

'''

import pandas as pd
import numpy as np
import argparse
import sys, os
import gzip
import re
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt



def process_header(line, clabel = None):

	# parse matrix header line

	mode_pat = re.compile(r'Mode:(\w+)')
	boundary_pat = re.compile(r'Boundary:(\d+)')
	bs_pat = re.compile(r'Binsize:(\d+)')
	ups_pat = re.compile(r'Upstream:(\d+)')
	dos_pat = re.compile(r'Downstream:(\d+)')
	mode = str(mode_pat.search(line).group(1))
	# check Mode
	if mode not in ['scale','refpoint']:
		sys.stderr.write('\n{}\tInvalid input file, please check the header line, the Mode should be "scale" or "refpoint".\n'.format(filename))
		sys.exit()

	bs = int(bs_pat.search(line).group(1))
	up_tick_loc = int(ups_pat.search(line).group(1)) / bs
	boundary = int(boundary_pat.search(line).group(1))

	x_ticks = [0, up_tick_loc-1, boundary-1]
	if mode == 'scale':
		down_tick_loc = int(dos_pat.search(line).group(1)) / bs
		x_ticks = [0, up_tick_loc-1, boundary-down_tick_loc-1, boundary-1]


	# make xlabels
	upstream = int(ups_pat.search(line).group(1))
	downstream = int(dos_pat.search(line).group(1))
	start_label = '+' + str(upstream / 1000) + 'kb'
	end_label = '-' + str(downstream / 1000) + 'kb'

	if mode == 'refpoint' and clabel == None:
		x_labels = [start_label, 'TSS', end_label]

	if mode == 'refpoint' and clabel != None:
		x_labels = [start_label, clabel[0], end_label]

	if mode == 'scale' and clabel == None:
		x_labels = [start_label, 'TSS', 'TES', end_label]

	if mode == 'scale' and clabel != None:
		x_labels = [start_label, clabel[0], clabel[1], end_label]


	return mode, boundary, x_ticks, x_labels


def read_in_file(infile):

	filename = str(infile)

	# read in header line and get some infomation
	mode_pat = re.compile(r'Mode:(\w+)')
	boundary_pat = re.compile(r'Boundary:(\d+)')

	with gzip.open(filename) as fi:
		header = fi.readline().decode()
		mode, boundary, x_ticks, x_labels = process_header(header)

	# create a list to save averaged values
	TS_res, NTS_res = [],[]
	ts_df = pd.read_table(filename,sep='\t',header=None,index_col=None,usecols=range(boundary),comment='#',dtype=float)
	nts_df = pd.read_table(filename,sep='\t',header=None,index_col=None,usecols=range(boundary,2*boundary),comment='#',dtype=float)

	# calculate confidence interval for each bins under 95% CL
	# append lower limit
	for i in range(boundary):
		TS_res.append(get_ci(np.array(ts_df.iloc[:,i]))[0])
		NTS_res.append(get_ci(np.array(nts_df.iloc[:,i]))[0])

	# append mean
	TS_res.extend(ts_df.apply(lambda x : x.mean(), axis=0).tolist())
	NTS_res.extend(nts_df.apply(lambda x : x.mean(), axis=0).tolist())

	# append upper limit
	for i in range(boundary):
		TS_res.append(get_ci(np.array(ts_df.iloc[:,i]))[1])
		NTS_res.append(get_ci(np.array(nts_df.iloc[:,i]))[1])

	return mode, boundary, TS_res, NTS_res, x_ticks, x_labels



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


def plot_w2c(data, boundary, ticks, labels, title, outfig, se=False):

	fig, ax = plt.subplots(figsize=(9,6), dpi=300)

	cut = boundary
	x = list(range(cut))

	assert len(data) == 4
	a,b,c,d = data

	# normalize 2 input
	norm_ts_low = np.asarray(a[:cut],dtype=float) / np.asarray(c[:cut],dtype=float)
	norm_ts_mean = np.asarray(a[cut:2*cut],dtype=float) / np.asarray(c[cut:2*cut],dtype=float)
	norm_ts_up = np.asarray(a[2*cut:],dtype=float) / np.asarray(c[2*cut:],dtype=float)

	norm_nts_low = np.asarray(b[:cut],dtype=float) / np.asarray(d[:cut],dtype=float)
	norm_nts_mean = np.asarray(b[cut:2*cut],dtype=float) / np.asarray(d[cut:2*cut],dtype=float)
	norm_nts_up = np.asarray(b[2*cut:],dtype=float) / np.asarray(d[2*cut:],dtype=float)


	ax.plot(x,norm_ts_mean, linestyle='-',lw=1.5,color='blue',label='TS')
	ax.plot(x,norm_nts_mean, linestyle='-',lw=1.5,color='red',label='NTS')
	if se:
		ax.fill_between(x,norm_ts_low,norm_ts_up,facecolor='blue',alpha=0.5)
		ax.fill_between(x,norm_nts_low,norm_nts_up,facecolor='red',alpha=0.5)

	ax.set_xlabel('')
	ax.set_ylabel('Normalized OG signal',fontdict={'size':15,'family':'DejaVu Serif','weight':600})
	ax.set_xlim(-1,cut)
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels,fontdict={'size':12,'family':'Verdana','weight':500})
	ax.set_title(title)
	ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5))
	
	fmt = str(outfig).split('.')[-1]
	plt.savefig(outfig,dpi=300,format=fmt,bbox_inches='tight')
	plt.close()
	return



def plot_wo2c(data, boundary, ticks, labels, title, outfig, se=False):

	fig, ax = plt.subplots(figsize=(9,6), dpi=300)

	cut = boundary
	x = list(range(cut))

	assert len(data) == 2
	a,b = data

	ts_low = np.asarray(a[:cut],dtype=float)
	ts_mean = np.asarray(a[cut:2*cut],dtype=float)
	ts_up = np.asarray(a[2*cut:],dtype=float)

	nts_low = np.asarray(b[:cut],dtype=float)
	nts_mean = np.asarray(b[cut:2*cut],dtype=float)
	nts_up = np.asarray(b[2*cut:],dtype=float)


	ax.plot(x,ts_mean, linestyle='-',lw=1.5,color='blue',label='TS')
	ax.plot(x,nts_mean, linestyle='-',lw=1.5,color='red',label='NTS')
	if se:
		ax.fill_between(x,ts_low,ts_up,facecolor='blue',alpha=0.5)
		ax.fill_between(x,nts_low,nts_up,facecolor='red',alpha=0.5)

	ax.set_xlabel('')
	ax.set_ylabel('OG signal',fontdict={'size':15,'family':'DejaVu Serif','weight':600})
	ax.set_xlim(-1,cut)
	ax.set_xticks(ticks)
	ax.set_xticklabels(labels,fontdict={'size':12,'family':'Verdana','weight':500})
	ax.set_title(title)
	ax.legend(loc='center left',bbox_to_anchor=(1.01,0.5))

	fmt = str(outfig).split('.')[-1]
	plt.savefig(outfig,dpi=300,format=fmt,bbox_inches='tight')
	plt.close()
	return


def pre_plot(data, boundary, ticks, labels, title, outfig, se, mode, n2c):

	if n2c:
		plot_w2c(data, boundary, ticks, labels, title, outfig, se)
	else:
		plot_wo2c(data, boundary, ticks, labels, title, outfig, se)



def main():

	# parse arguments
	parser = argparse.ArgumentParser(description=
		'Be used to plot stranded signal distribution with matrix[.gz] file'
		' generated from profileMatrix.py or scaleMatrix.py')

	parser.add_argument('-t', required=True, help='Matrix file [.gz] generated from profileMatrix.py')

	parser.add_argument('-c', required=False, help=
		'If control sample exists, then Treat sample signals would be normalized by this value.')

	parser.add_argument('-o', help='File name to save the image. '
		'The image format is depend on the suffix of input name.')

	parser.add_argument('--labels', nargs='+', type=str, help=
		'X-axis tick labels for reference point or scale region.'
		'Should match the type of input matrix.'
		'Default: TSS for refpoint, TSS TES for scale region.')

	parser.add_argument('--title', type=str, help='Set title for the figure.'
		'Default is None.')

	parser.add_argument('--se', type=bool, help=
		'Whether the final figure contains 95%% confidence interval. Default is False')

	args = parser.parse_args()



	# parse arguments
	seplot = False
	if args.se:
		seplot = True

	title = ' '
	if args.title:
		title = args.title

	outfig = args.o



	# load treat file and choose mode
	t_mode, t_boundary, t_ts, t_nts, t_x_ticks, t_x_labels = read_in_file(args.t)
	data = [t_ts, t_nts]

	flabels = t_x_labels[:]
	if args.labels:
		user_defined_labels = args.labels

		if t_mode == 'refpoint' and len(user_defined_labels) != 1:
			sys.stderr.write('Error ValueError:')
			sys.stderr.write('\nInput matrix is in reference point mode, but the input labels has {} elements.\n'.format(len(user_defined_labels)))
			sys.exit()

		if t_mode == 'scale' and len(user_defined_labels) != 2:
			sys.stderr.write('Error ValueError:')
			sys.stderr.write('\nInput matrix is in scale region mode, but the input labels has {} elements.\n'.format(len(user_defined_labels)))
			sys.exit()

		else:
			flabels = [t_x_labels[0]]
			for i in range(len(user_defined_labels)):
				flabels.append(user_defined_labels[i])
			flabels.append(t_x_labels[1])



	# check if control sample exists
	n2c = False
	if args.c:
		n2c = True
		c_mode, c_boundary, c_ts, c_nts, c_x_ticks, c_x_labels = read_in_file(args.c)

		try:
			assert t_mode==c_mode
		except AssertionError:
			sys.stderr.write('Error AssertionError:')
			sys.stderr.write('\nInput control matrix Mode is not same with Treat matrix.\n')
			sys.exit()

		try:
			assert t_boundary==c_boundary
		except AssertionError:
			sys.stderr.write('Error AssertionError:')
			sys.stderr.write('\nInput control matrix Boundary is not same with Treat matrix.\n')
			sys.exit()

		data.extend([c_ts, c_nts])


	# into plotting
	pre_plot(data, t_boundary, t_x_ticks, flabels, title, outfig, seplot, t_mode, n2c)

	# Done


if __name__ == '__main__':
	main()
