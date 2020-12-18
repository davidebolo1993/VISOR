#!/usr/bin/python env

#Python 3 standard library

import os
import sys
import re
import argparse
from argparse import HelpFormatter
from collections import defaultdict
from datetime import datetime

#additional modules

import pysam
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline import plot
from plotly import tools



class AutoVivification(dict):

	'''
	Fast way to create nested dictionaries
	'''

	def __getitem__(self, item):
		
		try:
			
			return dict.__getitem__(self, item)
		
		except KeyError:
			
			value = self[item] = type(self)()
			
			return value


class CustomFormat(HelpFormatter):

	'''
	Custom Help format
	'''

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



def WR(stranded_bam,chromosome,start,end):

	'''
	Extract Watson reads from strand-seq BAM
	'''

	WR = defaultdict(lambda: [None, None])

	for read in stranded_bam.fetch(chromosome,start,end):

		if read.is_proper_pair and not read.is_secondary and not read.is_supplementary:

			if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse): #read1 forward/ read2 reverse

				if read.query_name not in WR:

					if read.is_read1:

						WR[read.query_name][0] = read

					else:

						WR[read.query_name][1] = read

				else:

					if read.is_read1:

						yield read, WR[read.query_name][1]

					else:

						yield WR[read.query_name][0], read

					del WR[read.query_name]



def CR(stranded_bam,chromosome,start,end):

	'''
	Extract Crick reads from strand-seq BAM
	'''

	CR = defaultdict(lambda: [None, None])

	for read in stranded_bam.fetch(chromosome,start,end):

		if read.is_proper_pair and not read.is_secondary and not read.is_supplementary:

			if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse): #read2 forward and read1 reverse

				if read.query_name not in CR:

					if read.is_read1:

						CR[read.query_name][0] = read

					else:

						CR[read.query_name][1] = read

				else:

					if read.is_read1:

						yield read, CR[read.query_name][1]

					else:

						yield CR[read.query_name][0], read


					del CR[read.query_name]


def Count(bamfilein,binsize,chromosomes):

	'''
	Count Watons and Crick pairs
	'''

	bam = pysam.AlignmentFile(bamfilein, "rb", require_index=True)
	header=bam.header
	cdict=dict()
	chromosomes_info=list(header.items())[1][1]

	for i in range(len(chromosomes_info)):

		key,val=chromosomes_info[i].values()
		cdict[key]=val

	if chromosomes is not None: #in that case we want all the chromosomes in BAM header

		unwanted=set(cdict)-set(chromosomes[0])
		
		for unwanted_key in unwanted: 

			del cdict[unwanted_key]

	watsoncount=AutoVivification()
	crickcount=AutoVivification()
	
	for c,s in cdict.items():

		start=0
		chrom_end=s
		end=start+binsize

		while end <= chrom_end:

			listw=WR(bam,c,start,end)
			listc=CR(bam,c,start,end)

			watsoncount[c][start] = len(list(listw))
			crickcount[c][start] = len(list(listc))

			start+=binsize
			end+=binsize

	bam.close()

	return watsoncount,crickcount



def SSplot(wats_count,crick_count, output):

	'''
	Plot distribution of W/C reads along chromosome
	'''

	chromtraces=[]
	chromnames=list(wats_count.keys())

	for l,chroms in enumerate(wats_count.keys()):

		wats_chromdict=wats_count[chroms]
		crick_chromdict=crick_count[chroms]

		legend=True if l == 0 else False

		tracewatson=go.Bar(
			x = list(wats_chromdict.keys()),
			y = list(wats_chromdict.values()),
			marker = dict(color = '#F3A461'),
			text=list(wats_chromdict.values()),
			hoverinfo = 'text+name',
			name = 'watson',
			legendgroup='a',
			showlegend=legend)
				
		tracecrick= go.Bar(
			x = list(crick_chromdict.keys()),
			y = [-x for x in list(crick_chromdict.values())],
			marker = dict(color = '#668B8C'),
			text=list(crick_chromdict.values()),
			hoverinfo = 'text+name',
			name = 'crick',
			legendgroup='b',
			showlegend=legend
			)
			
		chromtraces.append((tracewatson,tracecrick))

	fig = tools.make_subplots(rows=len(chromtraces), cols=1,vertical_spacing=0.1,print_grid=False,subplot_titles=chromnames)
		
	for i in range(len(chromtraces)):


		fig.append_trace(chromtraces[i][0],i+1,1)
		fig.append_trace(chromtraces[i][1],i+1,1)

		if i == 0:

			fig['layout']['yaxis']['title'] = 'Count'
			#fig['layout']['xaxis']['title'] = 'Genomic position'

		elif 0 < i < len(chromtraces)-1 :

			fig['layout']['yaxis' + str(i+1)]['title'] = 'Count'
			#fig['layout']['xaxis' + str(i+1)]['title'] = 'Genomic position'

		else:

			fig['layout']['yaxis' + str(i+1)]['title'] = 'Count'
			fig['layout']['xaxis' + str(i+1)]['title'] = 'Genomic position'


	fig['layout'].update(height=600*len(chromtraces), title='Watson/Crick Counts')

	plot(fig, filename = output, auto_open=False)



def main():

	'''
	Execute the code and store the resulting plot to HTML

	'''

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This script is included in VISOR and was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-b','--bamfile', help='strand-seq BAM', metavar='BAM', required=True)
	required.add_argument('-o', '--output', help='output browsable HTML', metavar='HTML', required=True)

	optional = parser.add_argument_group('Additional parameters')

	optional.add_argument('--binsize', help= 'bin size for counting watson and crick reads [200000]', metavar='', default=200000, type=int)
	optional.add_argument('--chromosome', help= 'one or more chromosomes to plot. If None, plot all chromosomes from BAM header [None]', metavar='', nargs='+', action='append', default=None)

	args = parser.parse_args()


	BAM=os.path.abspath(args.bamfile)
	OUT=os.path.abspath(args.output)

	if not OUT.endswith('.html'):

		OUT = OUT+'.html'


	try:

		pysam.quickcheck(BAM)

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Errror] BAM ' + BAM + ' does not exist, is not readable or is not a valid BAM')
		sys.exit(1)


	if not os.path.exists(BAM + '.bai'):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Missing ' + BAM + ' index. Creating')

		try:

			pysam.index(BAM)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Errror] BAM ' + BAM + ' could not be indexed')
			sys.exit(1)	

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Counting Watson/Crick reads')

	watsoncount,crickcount = Count(BAM,args.binsize,args.chromosome)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Plotting')

	SSplot(watsoncount,crickcount,OUT)


	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)


if __name__ == '__main__':

	main()