#!/usr/bin/python env

#python 3 standard library

import os
import sys
import logging
from shutil import which
import subprocess
import timeit
import random
from collections import defaultdict
from multiprocessing import Process

#additional modules

import pybedtools
import pysam



def run(parser,args):

	#validate ouput

	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('It was not possible to create the results folder. Specify a path for which you have write permissions')
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK): #path exists but no write permissions on that folder

			print('You do not have write permissions on the directory in which results will be stored. Specify a folder for which you have write permissions')
			sys.exit(1)

	logging.basicConfig(filename=os.path.abspath(args.output + '/VISOR_SHORtS.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


	if not os.path.exists(os.path.abspath(args.output + '/simulations_haplotype1')):

		try:

			os.makedirs(os.path.abspath(args.output+ '/simulations_haplotype1'))

		except:

			logging.error('It was not possible to create haplotype 1 results folder. Specify a path for which you have write permissions')
			sys.exit(1)

			
	if not os.path.exists(os.path.abspath(args.output + '/simulations_haplotype2')):

		try:

			os.makedirs(os.path.abspath(args.output+ '/simulations_haplotype2'))

		except:

			logging.error('It was not possible to create haplotype 2 results folder. Specify a path for which you have write permissions')
			sys.exit(1)




	#check if external tools can be executed

	external_tools=['wgsim', 'bwa', 'samtools']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run VISOR SHORtS')
			sys.exit(1)



	#validate genome


	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) 

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)



	if not os.path.exists(os.path.abspath(args.genome + '.sa')):

		try:

			logging.info('Creating bwa index for reference genome')
			subprocess.check_call(['bwa', 'index', os.path.abspath(args.genome)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

		except:

			logging.error('It was not possible to generate bwa index for reference genome. Aborted')
			sys.exit(1)


	#validate h1.fa


	try:

		with open(os.path.abspath(args.hap1fa),'r') as file:

			assert(file.readline().startswith('>')) 

	except:

		logging.error('Specified haplotype 1 .fasta file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)


	if not os.path.exists(os.path.abspath(args.hap1fa + '.sa')) and args.type=='single-strand'

		try:

			runInParallel(BWA_Index,os.path.abspath(args.hap1fa),os.path.abspath(args.hap2fa))

		except:

			logging.error('Could not create bwa indexes. Aborted.')
			sys.exit(1)


	#validate h2.fa


	try:

		with open(os.path.abspath(args.hap2fa),'r') as file:

			assert(file.readline().startswith('>')) 

	except:

		logging.error('Specified haplotype 2 .fasta file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)




	#validate .bed1


	bedh1 = pybedtools.BedTool(os.path.abspath(args.hap1bed)) #this one is required

	
	try:
	
		srtbedh1 = bedh1.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		logging.error('Incorrect .bed format for haplotype 1')
		sys.exit(1)



	bedh2 = pybedtools.BedTool(os.path.abspath(args.hap2bed)) #this one is required


	#validate .bed2

	
	try:
	
		srtbedh2 = bedh1.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		logging.error('Incorrect .bed format for haplotype 2')
		sys.exit(1)




	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y', 'M']] #allowed chromosomes

	labels_seen_h1=set()
	labels_seen_h2=set()


	start=timeit.default_timer()

	logging.info('Simulating from haplotype 1')

	for entries in srtbedh1: #validate each entry

		if str(entries[0]) not in classic_chrs:

			logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file for haplotype 1. Allowed chromosomes are chr1-22, chrX, chrY and chrM')
			sys.exit(1)

		try:

			int(entries[1])

		except:

			logging.error('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file for haplotype 1. Start must be an integer')
			sys.exit(1)


		try:

			int(entries[2])

		except:

			logging.error('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file for haplotype 1. End must be an integer')
			sys.exit(1)


		if (int(entries[2]) - int(entries[1]) == 0):

			logging.error('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed for haplotype 1')
			sys.exit(1)


		if str(entries[3]) not in labels_seen_h1:


			labels_seen_h1.add(str(entries[3]))

			try:

				if args.type == 'double-strand':

					ClassicSimulate(os.path.abspath(args.genome), args.threads, os.path.abspath(args.hap1fa), str(entries[0]), int(entries[1]), int(entries[2]), str(entries[3]), args.error, args.coverage, args.length, args.indels, args.probability, os.path.abspath(args.output + '/simulations_haplotype1'))

				else:

					SSSimulate(args.threads, os.path.abspath(args.hap1fa), str(entries[0]), int(entries[1]), int(entries[2]), str(entries[3]), args.error, args.coverage, args.length, args.indels, args.probability, os.path.abspath(args.output + '/simulations_haplotype1'), args.noise)
					SingleStrand(os.path.abspath(args.genome), args.threads, os.path.abspath(args.output + '/simulations_haplotype1/' + str(entries[3]) + '.srt.bam'), str(entries[3]), args.noise, os.path.abspath(output + '/simulations_haplotype1'))

			except:

				logging.exception('Something went wrong during simulations for haplotype 1, ' + str(entries[0]) + ':' + str(entries[1]) + '-' + str(entries[2]) + '. Log is below.')



		else:

			logging.error('Multiple entries with same label in .bed file for haplotype 1')
			sys.exit(1)

	
	logging.info('Simulations from haplotype 1 completed')

	logging.info('Simulating from haplotype 2')

	for entries in srtbedh2: #validate each entry

		if str(entries[0]) not in classic_chrs:

			logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file for haplotype 2. Allowed chromosomes are chr1-22, chrX, chrY and chrM')
			sys.exit(1)

		try:

			int(entries[1])

		except:

			logging.error('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file for haplotype 2. Start must be an integer')
			sys.exit(1)


		try:

			int(entries[2])

		except:

			logging.error('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file for haplotype 2. End must be an integer')
			sys.exit(1)


		if (int(entries[2]) - int(entries[1]) == 0):

			logging.error('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed for haplotype 2')
			sys.exit(1)


		if str(entries[3]) not in labels_seen_h2:

			labels_seen_h2.add(str(entries[3]))

			try:

				if args.type == 'double-strand':

					ClassicSimulate(os.path.abspath(args.genome), args.threads, os.path.abspath(args.hap2fa), str(entries[0]), int(entries[1]), int(entries[2]), str(entries[3]), args.error, args.coverage, args.length, args.indels, args.probability, os.path.abspath(args.output + '/simulations_haplotype2'))

				else:

					SSSimulate(args.threads, os.path.abspath(args.hap2fa), str(entries[0]), int(entries[1]), int(entries[2]), str(entries[3]), args.error, args.coverage, args.length, args.indels, args.probability, os.path.abspath(args.output + '/simulations_haplotype2'))
					SingleStrand(os.path.abspath(args.genome), args.threads, os.path.abspath(args.output + '/simulations_haplotype2/' + str(entries[3]) + '.srt.bam'), str(entries[3]), args.noise, os.path.abspath(output + '/simulations_haplotype2'))


			except:

				logging.exception('Something went wrong during simulations for haplotype 1, ' + str(entries[0]) + ':' + str(entries[1]) + '-' + str(entries[2]) + '. Log is below.')

		else:

			logging.error('Multiple entries with same label in .bed file for haplotype 2')
			sys.exit(1)

	
	logging.info('Simulations from haplotype 2 completed')

	end=timeit.default_timer()
	elapsed=(end-start)/60
	logging.info('Simulations generated in ' + str(elapsed) + ' minutes')
	logging.info('Done')




def runInParallel(function, *arguments):

	proc = []

	for args in arguments:

		p = Process(target=function, args=args)
		p.start()
		proc.append(p)

	for p in proc:

		p.join()


def BWA_Index(fasta):

	subprocess.call(['bwa', 'index', os.path.abspath(fasta)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))



def ClassicSimulate(genome, cores, haplotype, chromosome, start, end, label, error, coverage, length, indels, probability, output):


	#prepare region

	with open(os.path.abspath(output + '/region.tmp.fa'), 'w') as regionout:

		subprocess.call(['samtools', 'faidx', haplotype, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))

	numreads= round((coverage*(end-start)) / length) 

	#simulate reads

	subprocess.call(['wgsim', '-e', str(error), '-N', str(numreads), '-1', str(length), '-2', str(length), '-R', str(indels), '-X', str(probability), os.path.abspath(output + '/region.tmp.fa'), os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.fa'))

	#align to modified reference

	with open(os.path.abspath(output + '/region.tmp.sam'), 'w') as samout:

		subprocess.call(['bwa', 'mem', '-t', str(cores), reference, os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stdout=samout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.1.fq'))
	os.remove(os.path.abspath(output + '/region.2.fq'))

	with open(os.path.abspath(output + '/region.tmp.bam'), 'w') as bamout:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/region.tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.sam'))

	with open(os.path.abspath(output + '/' + label + '.srt.bam'), 'w') as srtbamout:

		subprocess.call(['samtools', 'sort', os.path.abspath(output + '/region.tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.bam'))

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))



def SSSimulate(cores, haplotype, chromosome, start, end, label, error, coverage, length, indels, probability, output, noise):


	#prepare region

	with open(os.path.abspath(output + '/region.tmp.fa'), 'w') as regionout:

		subprocess.call(['samtools', 'faidx', haplotype, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))

	numreads= round((coverage*(end-start)) / length) 

	#simulate reads

	subprocess.call(['wgsim', '-e', str(error), '-N', str(numreads), '-1', str(length), '-2', str(length), '-R', str(indels), '-X', str(probability), os.path.abspath(output + '/region.tmp.fa'), os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.fa'))

	#align to modified reference


	with open(os.path.abspath(output + '/region.tmp.sam'), 'w') as samout:

		subprocess.call(['bwa', 'mem', '-t', str(cores), haplotype, os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stdout=samout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.1.fq'))
	os.remove(os.path.abspath(output + '/region.2.fq'))

	with open(os.path.abspath(output + '/region.tmp.bam'), 'w') as bamout:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/region.tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.sam'))

	with open(os.path.abspath(output + '/' + label + '.srt.bam'), 'w') as srtbamout:

		subprocess.call(['samtools', 'sort', os.path.abspath(output + '/region.tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.bam'))

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))




def SingleStrand(genome, cores, bamfilein, label, noisefraction, output):

	bam = pysam.AlignmentFile(bamfilein, "rb")

	watsonbam = pysam.AlignmentFile(os.path.abspath(output + '/' + label + '.watson.bam'), "wb", template=bam)

	watslist=list(watson_orientation(bam))

	for read1,read2 in watslist:
			
		watsonbam.write(read1)
		watsonbam.write(read2)

	
	crickbam = pysam.AlignmentFile(os.path.abspath(output + '/' + label + '.crick.bam'), "wb", template=bam)

	cricklist=list(crick_orientation(bam))

	for read1,read2 in cricklist:
			
		crickbam.write(read1)
		crickbam.write(read2)


	if noisefraction > 0:

		discordant_pair_watson= round((len(watslist)*noisefraction)/100)
		discordant_pair_crick= round((len(cricklist)*noisefraction)/100)

		sample_crick=random.sample(cricklist,discordant_pair_watson)

		for read1,read2 in sample_crick:

			watsonbam.write(read1)
			watsonbam.write(read2)

		sample_watson=random.sample(watslist,discordant_pair_crick)

		for read1,read2 in sample_watson:

			crickbam.write(read1)
			crickbam.write(read2)

	watsonbam.close()
	crickbam.close()
	bam.close()

	os.remove(bamfilein)
	os.remove(bamfilein + '.bai')


	with open(os.path.abspath(output + '/' + label + '.watson.fq'), 'w') as watsonfq:

		subprocess.call(['samtools', 'fastq', os.path.abspath(output + '/' + label + '.watson.bam')], stdout=watsonfq, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/' + label + '.watson.bam'))

	#subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.watson.srt.bam')],stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/' + label + '.crick.fq'), 'w') as crickfq:

		subprocess.call(['samtools', 'fastq', os.path.abspath(output + '/' + label + '.crick.bam')], stdout=watsonfq, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/' + label + '.crick.bam'))
	
	
	#subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.crick.srt.bam')],stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/watson.tmp.sam'), 'w') as watsonsam:

		subprocess.call(['bwa', 'mem', '-t', str(cores), genome, os.path.abspath(output + '/' + label + '.watson.fq')], stdout=watsonsam, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/crick.tmp.sam'), 'w') as cricksam:

		subprocess.call(['bwa', 'mem', '-t', str(cores), genome, os.path.abspath(output + '/' + label + '.crick.fq')], stdout=cricksam, stderr=open(os.devnull, 'wb'))


	os.remove(os.path.abspath(output + '/' + label + '.watson.fq'))
	os.remove(os.path.abspath(output + '/' + label + '.crick.fq'))

	with open(os.path.abspath(output + '/watson.tmp.bam'), 'w') as watsonbam:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/watson.tmp.sam')], stdout=watsonbam, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/crick.tmp.bam'), 'w') as crickbam:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/crick.tmp.sam')], stdout=crickbam, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/watson.tmp.sam'))
	os.remove(os.path.abspath(output + '/crick.tmp.sam'))


	with open(os.path.abspath(output + '/' + label + 'watson.srt.bam'), 'w') as watsonsort:

		subprocess.call(['samtools', 'sort', os.path.abspath(output + '/watson.tmp.bam')], stdout=watsonsort, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/' + label + 'crick.srt.bam'), 'w') as cricksort:

		subprocess.call(['samtools', 'sort', os.path.abspath(output + '/crick.tmp.bam')], stdout=cricksort, stderr=open(os.devnull, 'wb'))


	os.remove(os.path.abspath(output + '/watson.tmp.bam'))
	os.remove(os.path.abspath(output + '/crick.tmp.bam'))


	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + 'watson.srt.bam')],stderr=open(os.devnull, 'wb'))
	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + 'crick.srt.bam')],stderr=open(os.devnull, 'wb'))





def watson_orientation(bam):

	read_dict = defaultdict(lambda: [None, None])

	for read in bam.fetch():

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





def crick_orientation(bam):

	read_dict = defaultdict(lambda: [None, None])

	for read in bam.fetch():

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
