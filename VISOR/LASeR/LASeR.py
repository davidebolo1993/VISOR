#!/usr/bin/python env

#python 3 standard library

import os
import sys
import logging
from shutil import which
import subprocess
import glob

#additional modules

import pybedtools
import pyfaidx



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

		if os.listdir(os.path.abspath(args.output)):

			print('Specified output directory is not empty. Specify another directory or clean the chosen one')
			sys.exit(1)


	logging.basicConfig(filename=os.path.abspath(args.output + '/VISOR_LASeR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


	if not os.access(os.path.dirname(os.path.abspath(args.genome)),os.W_OK):

		logging.error('It is required to have write access on reference folder to generate minimap2 .mmi indexes')
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



	fastaslist=args.haplotypefasta[0]


	for fastas in fastaslist:

		try:

			with open(os.path.abspath(fastas),'r') as file:

				assert(file.readline().startswith('>')) 

		except:

			logging.error(os.path.abspath(fastas) + ' file does not exist, is not readable or is not a valid .fasta file')
			sys.exit(1)



	#validate .bed

	bed = pybedtools.BedTool(os.path.abspath(args.bedfile)) #this one is required

	
	try:

		srtbed = bed.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		logging.error('Incorrect .bed format for -bed/--bedfile')
		sys.exit(1)


	allelic=args.allelicfraction
	fa=pyfaidx.Fasta(os.path.abspath(args.genome))
	classic_chrs = fa.keys() #allowed chromosomes

	model_qc=os.path.abspath(os.path.dirname(__file__) + '/model_qc_clr')

	logging.info('Simulations starts now')


	for fastas in fastaslist:

		haploname=os.path.basename(os.path.abspath(fastas)).split('.')[0]

		counter = 0

		for entries in srtbed: #validate each entry

			counter +=1
			
			if str(entries[0]) not in classic_chrs:

				logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file.')
				sys.exit(1)

			try:

				int(entries[1])

			except:

				logging.error('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file. Start must be an integer')
				sys.exit(1)


			try:

				int(entries[2])

			except:

				logging.error('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file. End must be an integer')
				sys.exit(1)


			if (int(entries[2]) - int(entries[1]) == 0):

				logging.error('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed file')
				sys.exit(1)


			try:

				float(entries[3])

			except:

				logging.error('Cannot convert ' + str(entries[3]) + ' to float number in .bed file. Coverage bias must be a float')
				sys.exit(1)

			try:

				Simulate(os.path.abspath(args.genome), args.threads, os.path.abspath(args.hap1fa), str(entries[0]), int(entries[1]), int(entries[2]), args.identifier + '.' + str(counter), model_qc, args.accuracy, (args.coverage / 100 * float(entries[3])), allelic, args.length, args.ratio, os.path.abspath(args.output + '/simulations_' +haploname))

			except:

				logging.exception('Something went wrong during simulations for ' + os.path.abspath(fastas) + '. Log is below.')

	

		if counter == 1:


			os.rename(os.path.abspath(args.output + '/simulations_' + haploname + '/' + args.identifier + '.' + str(counter) + '.srt.bam'), os.path.abspath(args.output + '/simulations_' + haploname + '/' +  args.identifier + '.srt.bam'))
			os.rename(os.path.abspath(args.output + '/simulations_' + haploname + '/' + args.identifier + '.' + str(counter) + '.srt.bam.bai'), os.path.abspath(args.output + '/simulations_' + haploname + '/' +  args.identifier + '.srt.bam.bai'))

		else:


			bams = glob.glob(os.path.abspath(args.output + '/simulations_' + haploname + '/' + '*.srt.bam'))

			with open(os.path.abspath(args.output + '/simulations_' + haploname + '/bamtomerge.txt'), 'a') as bamout:

				for file in bams:

					bamout.write(file + '\n')

			subprocess.call(['samtools', 'merge', '-b', os.path.abspath(args.output + '/simulations_' + haploname + '/bamtomerge.txt'), os.path.abspath(args.output + '/simulations_' + haploname + '/' + args.identifier + '.srt.bam')])
			os.remove(os.path.abspath(args.output + '/simulations_' + haploname + '/bamtomerge.txt'))
			subprocess.call(['samtools', 'index', os.path.abspath(args.output + '/simulations_' + haploname + '/' + args.identifier + '.srt.bam')])


			for b in bams:

				os.remove(b)
				os.remove(b + '.bai')


	logging.info('Done')




def Simulate(genome, cores, haplotype, chromosome, start, end, label, model_qc, accuracy, coverage, allelic, length, ratio, output):



	#prepare region

	with open(os.path.abspath(output + '/region.tmp.fa'), 'w') as regionout:

		subprocess.call(['samtools', 'faidx', haplotype, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))



	if not allelic == 100: #simulate part from modified and part from reference

		with open(os.path.abspath(output + '/reference.region.tmp.fa'), 'w') as regionout:

			subprocess.call(['samtools', 'faidx', genome, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))


		coveragevar = (coverage/100)*allelic
		coverageref = coverage - coveragevar


		subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', output + '/simref','--length-mean', str(length), '--accuracy-mean', str(accuracy), '--difference-ratio', ratio, '--depth', str(coverageref), os.path.abspath(output + '/reference.region.tmp.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(output + '/reference.region.tmp.fa'))
		os.remove(os.path.abspath(output + '/simref_0001.ref'))
		os.remove(os.path.abspath(output + '/simref_0001.maf'))


		subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', output + '/simvar','--length-mean', str(length), '--accuracy-mean', str(accuracy), '--difference-ratio', ratio, '--depth', str(coveragevar), os.path.abspath(output + '/region.tmp.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(output + '/region.tmp.fa'))
		os.remove(os.path.abspath(output + '/simvar_0001.ref'))
		os.remove(os.path.abspath(output + '/simvar_0001.maf'))


		with open (os.path.abspath(output + '/sim_0001.fastq'), 'w') as regionout:

			subprocess.call(['cat', os.path.abspath(output + '/simref_0001.fastq'), os.path.abspath(output + '/simvar_0001.fastq')], stdout=regionout, stderr=open(os.devnull, 'wb'))


		os.remove(os.path.abspath(output + '/simref_0001.fastq'))
		os.remove(os.path.abspath(output + '/simvar_0001.fastq'))


	else:

		subprocess.call(['pbsim', '--model_qc', model_qc, '--prefix', output + '/sim','--length-mean', str(length), '--accuracy-mean', str(accuracy), '--difference-ratio', ratio, '--depth', str(coverage), os.path.abspath(output + '/region.tmp.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(output + '/region.tmp.fa'))
		os.remove(os.path.abspath(output + '/sim_0001.ref'))
		os.remove(os.path.abspath(output + '/sim_0001.maf'))



	if not os.path.exists(os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.mmi')):

		if not os.path.exists(os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.fa')): #checks for the presence of a .fa ref for the chromosome in the reference folder. If not present, creates it. Will save time during alignments.

			with open(os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.fa'),'w') as fout: 

				subprocess.call(['samtools', 'faidx', genome, chromosome], stdout=fout, stderr=open(os.devnull, 'wb'))

		subprocess.call(['minimap2', '-d', os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.mmi'), os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.fa')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb')) #create .mmi: faster when it comes to simulate multiple times from same chromosome

	new_mmi=os.path.abspath(os.path.dirname(genome) +'/' + chromosome + '.mmi')

	#align to reference

	with open(os.path.abspath(output + '/region.tmp.sam'), 'w') as samout:

		subprocess.call(['minimap2', '-ax', 'map-ont', '-t', str(cores), new_mmi, os.path.abspath(output + '/sim_0001.fastq')], stdout=samout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/sim_0001.fastq'))

	with open(os.path.abspath(output + '/region.tmp.bam'), 'w') as bamout:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/region.tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.sam'))

	with open(os.path.abspath(output + '/' + label + '.srt.bam'), 'w') as srtbamout:

		subprocess.call(['samtools', 'sort', os.path.abspath(output + '/region.tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.bam'))

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))
