#!/usr/bin/python env

#Python 3 standard library

from __future__ import print_function
import os
import sys
from collections import defaultdict
from shutil import which
import argparse
from argparse import HelpFormatter
import subprocess
import re

#additional modules

import pysam
import pyfaidx
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot
from plotly import tools


def main():


	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This script is included in VISOR and was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa', required=True)
	required.add_argument('-bam', '--bamfile', help='single-strand .bam file/s (crick or watson or crick + watson) to plot', metavar='.bam', nargs='*', action='append', required=True)
	required.add_argument('-O', '--output', help='where plots will be saved', metavar='folder', required=True)


	optional = parser.add_argument_group('Bin size')

	optional.add_argument('-bin', '--binsize', help= 'bin size for counting [200000]', metavar='', default=200000, type=int)
	optional.add_argument('-l', '--label', help= 'label/s to identify sample/s [None]', metavar='', nargs='*', action='append', default=None)


	args = parser.parse_args()


	if not os.path.exists(os.path.abspath(args.output + '/comparison')):

		try:

			os.makedirs(os.path.abspath(args.output + '/comparison'))

		except:

			print('It was not possible to create the results folder. Specify a path for which you have write permissions')
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.dirname(os.path.abspath(args.output + '/comparison')),os.W_OK): #path exists but no write permissions on that folder

			print('You do not have write permissions on the directory in which results will be stored. Specify a folder for which you have write permissions')
			sys.exit(1)



	if args.label is not None:

		try:
			
			assert(len(args.label[0]) == len(args.bamfile[0]))
			labels=args.label[0]

		except:

			print('Different number of samples and labels')
			sys.exit(1)

	else:

		labels=[]

		for i in range(len(args.bamfile[0])):

			labels.append('.bam ' + str(i+1))

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) 

	except:

		print('Reference file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)


	if which('samtools') is None:

		print('Samtools is required')
		sys.exit(1)

	else:

		for bamfiles in args.bamfile[0]:

			try:

				subprocess.check_call(['samtools', 'quickcheck', os.path.abspath(bamfiles)])

			except:

				print(os.path.abspath(bamfiles), 'does not exist, is not readable or is not a valid .bam file')
				sys.exit(1)


	reslist=[]

	progress = ProgressBar(len(args.bamfile[0]), fmt=ProgressBar.FULL)

	print('Counting')

	for bams in args.bamfile[0]:

		watscount, crickcount = Count(os.path.abspath(args.genome), os.path.abspath(bams), args.binsize)
		reslist.append((watscount, crickcount))
		
		progress.current += 1
		progress()

	progress.done()

	print('Done')
	print('Plotting')


	Plotter(reslist,labels, os.path.abspath(args.output))

	print('Done')


	sys.exit(0)



class CustomFormat(HelpFormatter):

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()



class AutoVivification(dict):

	def __getitem__(self, item):
		
		try:
			
			return dict.__getitem__(self, item)
		
		except KeyError:
			
			value = self[item] = type(self)()
			
			return value



class ProgressBar(object):

	DEFAULT = 'Progress: %(bar)s %(percent)3d%%'
	FULL = '%(bar)s %(current)d/%(total)d (%(percent)3d%%) %(remaining)d to go'

	def __init__(self, total, width=40, fmt=DEFAULT, symbol='=',output=sys.stderr):
		
		assert len(symbol) == 1

		self.total = total
		self.width = width
		self.symbol = symbol
		self.output = output
		self.fmt = re.sub(r'(?P<name>%\(.+?\))d',r'\g<name>%dd' % len(str(total)), fmt)
		self.current = 0

	def __call__(self):

		percent = self.current / float(self.total)
		size = int(self.width * percent)
		remaining = self.total - self.current
		bar = '[' + self.symbol * size + ' ' * (self.width - size) + ']'
		args = {
			'total': self.total,
			'bar': bar,
			'current': self.current,
			'percent': percent * 100,
			'remaining': remaining
		}

		print('\r' + self.fmt % args, file=self.output, end='')

	def done(self):

		self.current = self.total
		self()
		print('', file=self.output)




def Plotter(reslist, labels, output):

	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y', 'M']] #allowed chromosomes

	for chroms in classic_chrs:

		chromtraces=[]

		for i,res in enumerate(reslist):

			watson_trace=res[0] #watson trace for bam file
			crick_trace=res[1] #crick trace for bam file

			wats_chromdict=watson_trace[chroms]
			crick_chromdict=crick_trace[chroms]

			if i == 0:

				tracewatson=go.Bar(
					x = list(wats_chromdict.keys()),
					y = list(wats_chromdict.values()),
					marker = dict(color = 'rgb(55, 83, 109)'),
					text=list(wats_chromdict.values()),
					hoverinfo = 'text+name',
					name = 'watson',
					legendgroup='a')
				
				tracecrick= go.Bar(
					x = list(crick_chromdict.keys()),
					y = [-x for x in list(crick_chromdict.values())],
					marker = dict(color = '#B35F5F'),
					text=list(crick_chromdict.values()),
					hoverinfo = 'text+name',
					name = 'crick',
					legendgroup='b'
					)

			else:

				tracewatson=go.Bar(
					x = list(wats_chromdict.keys()),
					y = list(wats_chromdict.values()),
					marker = dict(color = 'rgb(55, 83, 109)'),
					text=list(wats_chromdict.values()),
					hoverinfo = 'text+name',
					name = 'watson',					
					showlegend = False,
					legendgroup='a')
				
				tracecrick= go.Bar(
					x = list(crick_chromdict.keys()),
					y = [-x for x in list(crick_chromdict.values())],
					marker = dict(color = '#B35F5F'),
					text=list(crick_chromdict.values()),
					hoverinfo = 'text+name',
					name = 'crick',					
					showlegend = False,
					legendgroup='b'
					)
			
			
			chromtraces.append((tracewatson,tracecrick))


		fig = tools.make_subplots(rows=len(chromtraces), cols=1, subplot_titles=labels,vertical_spacing=0.1)
		
		for i in range(len(chromtraces)):

			wats=chromtraces[i][0]
			crick=chromtraces[i][1]

			fig.append_trace(wats, i+1,1)
			fig.append_trace(crick, i+1, 1)

		fig['layout'].update(height=600*len(chromtraces), title=chroms)

		for i in range(len(chromtraces)):

			if i == 0:

				fig['layout']['yaxis']['title'] = 'Count'
				fig['layout']['xaxis']['title'] = 'Genomic position'

			else:

				fig['layout']['yaxis' + str(i+1)]['title'] = 'Count'
				fig['layout']['xaxis' + str(i+1)]['title'] = 'Genomic position'


		plot(fig, filename = os.path.abspath(output + '/' + chroms + '.html'), auto_open=False)




def Count(genomein, bamfilein, binsize):


	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y', 'M']] #allowed chromosomes


	watscount=AutoVivification()
	crickcount=AutoVivification()
	
	bam = pysam.AlignmentFile(os.path.abspath(bamfilein), "rb")
	genome=pyfaidx.Fasta(os.path.abspath(genomein))


	for chrs in classic_chrs:


		chrom=genome[chrs]
		start=0
		chrom_end=len(chrom[:len(chrom)].seq)
		end= start + binsize

		while start <= chrom_end:

			listw=watson(bam, chrs, start, end)
			listc=crick(bam, chrs, start, end)

			watscount[chrs][start] = len(list(listw))
			crickcount[chrs][start] = len(list(listc))

			start += binsize
			end += binsize

			if chrom_end - start < binsize:

				end == chrom_end


	bam.close()

	return watscount, crickcount




def watson(bam, chromosome, start, end):

	read_dict = defaultdict(lambda: [None, None])

	for read in bam.fetch(chromosome, start, end):

		if not read.is_proper_pair or read.is_secondary or read.is_supplementary:

			continue

		elif read.is_read1 and read.is_reverse: #if read1 is reverse skip

			continue

		elif read.is_read2 and not read.is_reverse: #if read2 is not reverse skip

			continue

		else: #read 1 is forward and read 2 is reverse

			qname = read.query_name

			if qname not in read_dict:

				if read.is_read1:

					read_dict[qname][0] = read

				else:

					read_dict[qname][1] = read
			else:

				if read.is_read1:

					yield read, read_dict[qname][1]
				
				else:

					yield read_dict[qname][0], read

				del read_dict[qname]





def crick(bam, chromosome, start, end):

	read_dict = defaultdict(lambda: [None, None])

	for read in bam.fetch(chromosome, start, end):

		if not read.is_proper_pair or read.is_secondary or read.is_supplementary:

			continue

		elif read.is_read1 and not read.is_reverse: #if read1 is not reverse skip

			continue

		elif read.is_read2 and read.is_reverse: #if read2 is reverse skip

			continue

		else:

			qname = read.query_name

			if qname not in read_dict:

				if read.is_read1:

					read_dict[qname][0] = read

				else:

					read_dict[qname][1] = read
			else:

				if read.is_read1:

					yield read, read_dict[qname][1]
				
				else:

					yield read_dict[qname][0], read

				del read_dict[qname]





if __name__ == '__main__':

	main()
