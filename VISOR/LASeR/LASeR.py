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

	logging.basicConfig(filename=os.path.abspath(args.output + '/VISOR_LASeR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


	if not os.access(os.path.dirname(os.path.abspath(args.genome)),os.W_OK):

		logging.error('It is required to have write access on reference folder to generate minimap2 .mmi indexes')
		sys.exit(1)


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

	external_tools=['pbsim', 'minimap2', 'samtools']

	for tools in external_tools:

		if which(tools) is None:

			logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run VISOR LASeR')
			sys.exit(1)



	#validate genome


	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) 

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)



	#validate h1.fa


	try:

		with open(os.path.abspath(args.hap1fa),'r') as file:

			assert(file.readline().startswith('>')) 

	except:

		logging.error('Specified haplotype 1 .fasta file does not exist, is not readable or is not a valid .fasta file')
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

	model_qc=os.path.abspath(os.path.dirname(__file__) + '/model_qc_clr')

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

				Simulate(os.path.abspath(args.genome), args.threads, os.path.abspath(args.hap1fa), str(entries[0]), int(entries[1]), int(entries[2]), str(entries[3]), model_qc, args.accuracy, args.coverage, args.length, args.ratio, os.path.abspath(args.output + '/simulations_haplotype1'))

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

				Simulate(os.path.abspath(args.genome), args.threads, os.path.abspath(args.hap2fa), str(entries[0]), int(entries[1]), int(entries[2]), str(entries[3]), model_qc, args.accuracy, args.coverage, args.length, args.ratio, os.path.abspath(args.output + '/simulations_haplotype2'))

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




def Simulate(genome, cores, haplotype, chromosome, start, end, label, model_qc, accuracy, coverage, length, ratio, output):


	#prepare region

	with open(os.path.abspath(output + '/region.tmp.fa'), 'w') as regionout:

		subprocess.call(['samtools', 'faidx', haplotype, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))


	#simulate reads

	subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', 'sim','--length-mean', str(length), '--accuracy-mean', str(accuracy), '--difference-ratio', ratio, '--depth', str(coverage), os.path.abspath(output + '/region.tmp.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))


	os.remove(os.path.abspath(output + '/region.tmp.fa'))
	os.remove(os.path.abspath(output + '/sim_0001.ref'))
	os.remove(os.path.abspath(output + '/sim_0001.maf'))


	if not os.path.exists(os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.mmi')):

		if not os.path.exists(os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.fa')): #checks for the presence of a .fa ref for the chromosome in the reference folder. If not present, creates it. Will save time during alignments.

			with open(os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa'),'w') as fout: 

				subprocess.call(['samtools', 'faidx', reference, chromosome], stdout=fout, stderr=open(os.devnull, 'wb'))

		subprocess.call(['minimap2', '-d', os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.mmi'), os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb')) #create .mmi: faster when it comes to simulate multiple times from same chromosome

	new_mmi=os.path.abspath(os.path.dirname(reference) +'/' + chromosome + '.mmi')

	#align to reference

	with open(os.path.abspath(output + '/region.tmp.sam'), 'w') as samout:

		subprocess.call(['minimap2', '-ax', 'map-ont', new_mmi, os.path.abspath(output + '/sim_0001.fastq')], stdout=samout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/sim_0001.fastq'))

	with open(os.path.abspath(output + '/region.tmp.bam'), 'w') as bamout:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/region.tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.sam'))

	with open(os.path.abspath(output + '/' + label + '.srt.bam'), 'w') as srtbamout:

		subprocess.call(['samtools', 'sort', os.path.abspath(output + '/region.tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.bam'))

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))
