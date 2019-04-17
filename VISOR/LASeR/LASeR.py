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

			print('It was not possible to create the output folder. Specify a path for which you have write permissions')
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.dirname(os.path.abspath(args.output)),os.W_OK): #path exists but no write permissions on that folder

			print('You do not have write permissions on the output folder. Specify a folder for which you have write permissions')
			sys.exit(1)


		if os.listdir(os.path.abspath(args.output)):

			print('Specified output folder is not empty. Specify another directory or clean the chosen one')
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



	#validate .bed

	bed = pybedtools.BedTool(os.path.abspath(args.bedfile)) #this one is required

	
	try:

		srtbed = bed.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		logging.error('Incorrect .bed format for -bed/--bedfile')
		sys.exit(1)


	inputs=args.sample[0]

	if len(inputs) > 1:

		if args.clonefraction is None:

			logging.error('When specifying multiple -s/--sample, multiple -cf/--clonefraction percentages must be specified')
			sys.exit(1)

		else: #something has been specified

			fractions=args.clonefraction[0]

			if len(fractions) != len(inputs):

				logging.error('When specifying multiple -s/--sample, the same number of -cf/--clonefraction percentages must be specified')
				sys.exit(1)

			for fraction in fractions:

				try:

					float(fraction)

				except:

					logging.error('Each fraction percentage in -cf/--clonefraction must be float')
					sys.exit(1)


			totalfraction = sum(map(float,fractions))

			if totalfraction > 100:

				logging.error('Sum of fractions percentages in -cf/--clonefraction cannot exceed 100.0')
				sys.exit(1)


	fa=pyfaidx.Fasta(os.path.abspath(args.genome))
	classic_chrs = fa.keys() #allowed chromosomes
	model_qc=os.path.abspath(os.path.dirname(__file__) + '/model_qc_clr')

	logging.info('Running simulations')

	if len(inputs) == 1: #just one folder, use a classic simulation

		logging.info('Single input for -s/--sample')

		#find .fasta in folder

		fastas = glob.glob(os.path.abspath(inputs[0] + '/*.fa'))

		if fastas == []:

			logging.error('Given folder ' + inputs[0] + ' does not contain any valid .fasta inputs')
			sys.exit(1)


		for folder,fasta in enumerate(fastas):

			os.makedirs(os.path.abspath(args.output + '/' + str(folder))) #create directory for this haplotype

			counter=0

			for entries in srtbed: #validate each entry

				counter+=1
				
				if str(entries[0]) not in classic_chrs:

					logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file')
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

					logging.error('Cannot convert ' + str(entries[3]) + ' to float number in .bed file. Capture bias must be a float')
					sys.exit(1)


				try:

					allelic=float(entries[4])

				except:

					logging.error('Cannot convert ' + str(entries[4]) + ' to float number in .bed file. Sample fraction must be a float percentage')
					sys.exit(1)

				try:

					Simulate(os.path.abspath(args.genome), args.threads, os.path.abspath(fasta), str(entries[0]), int(entries[1]), int(entries[2]), str(counter), model_qc, args.accuracy, (args.coverage / 100 * float(entries[3]))/len(fastas), allelic, args.length, args.ratio, os.path.abspath(args.output + '/' + str(folder)),folder +1, 1)
				
				except:
					
					logging.exception('Something went wrong during simulations for ' + os.path.abspath(fasta) + '. Log is below.')


		subdirs=[os.path.join(os.path.abspath(args.output), o) for o in os.listdir(os.path.abspath(args.output)) if os.path.isdir(os.path.join(os.path.abspath(args.output),o))]
		bams = [y for x in os.walk(os.path.abspath(args.output)) for y in glob.glob(os.path.join(x[0], '*.srt.bam'))]

		with open(os.path.abspath(args.output + '/bamstomerge.txt'), 'w') as bamstomerge:

			for bam in bams:

				bamstomerge.write(bam + '\n')

		subprocess.call(['samtools', 'merge', '-b', os.path.abspath(args.output + '/bamstomerge.txt'), os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')], stderr=open(os.devnull, 'wb'))
		subprocess.call(['samtools', 'index', os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')],stderr=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(args.output + '/bamstomerge.txt'))

		for bam in bams:

			os.remove(bam)
			os.remove(bam + '.bai')

		for dirs in subdirs: #now they are empty and vcan be removed safely

			os.rmdir(dirs)

	else: # simulate subclones

		logging.info('Multiple inputs for -s/--sample. Assuming each input is a subclone')

		for fract,inp in enumerate(inputs): #each input is now a subclone

			os.makedirs(os.path.abspath(args.output + '/' + str(fract)))

			subfastas=glob.glob(os.path.abspath(inp) + '/*.fa')
			subfastasfraction= float(fractions[fract]) #percentage of this clone in total in the final .bam
			eachhaplofraction=subfastasfraction/len(subfastas)

			for folder,subfasta in enumerate(subfastas):

				os.makedirs(os.path.abspath(args.output + '/' + str(fract) + '/' + str(folder)))

				counter=0

				for entries in srtbed: #validate each entry

					counter +=1
					
					if str(entries[0]) not in classic_chrs:

						logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file')
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

						logging.error('Cannot convert ' + str(entries[3]) + ' to float number in .bed file. Capture bias must be a float')
						sys.exit(1)


					try:		

						Simulate(os.path.abspath(args.genome), args.threads, os.path.abspath(subfasta), str(entries[0]), int(entries[1]), int(entries[2]), str(counter), model_qc, args.accuracy, ((args.coverage / 100 * float(entries[3]))/100)*eachhaplofraction, 100.0, args.length, args.ratio, os.path.abspath(args.output + '/' + str(fract) + '/' + str(folder)), folder+1, fract+1)

					except:

						logging.exception('Something went wrong during simulations for ' + os.path.abspath(subfasta) + '. Log is below.')

		subdirs=[os.path.join(os.path.abspath(args.output), o) for o in os.listdir(os.path.abspath(args.output)) if os.path.isdir(os.path.join(os.path.abspath(args.output),o))]
		bams = [y for x in os.walk(os.path.abspath(args.output)) for y in glob.glob(os.path.join(x[0], '*.srt.bam'))]


		with open(os.path.abspath(args.output + '/bamstomerge.txt'), 'w') as bamstomerge:

			for bam in bams:

				bamstomerge.write(bam + '\n')


		subprocess.call(['samtools', 'merge', '-b', os.path.abspath(args.output + '/bamstomerge.txt'), os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')], stderr=open(os.devnull, 'wb'))
		subprocess.call(['samtools', 'index', os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')],stderr=open(os.devnull, 'wb'))

		os.remove(os.path.abspath(args.output + '/bamstomerge.txt'))

		for bam in bams:

			os.remove(bam)
			os.remove(bam + '.bai')


		for subs in subdirs:

			subsub=[os.path.join(os.path.abspath(subs), o) for o in os.listdir(os.path.abspath(subs)) if os.path.isdir(os.path.join(os.path.abspath(subs),o))]

			for s in subsub: #more safe than other solutions

				os.rmdir(s)

			os.rmdir(subs)

	logging.info('Done')




def ModifyReadTags(inbam, haplonum, clone):

	bam=pysam.AlignmentFile(os.path.abspath(inbam), 'rb')

	outbam=pysam.AlignmentFile(os.path.abspath(inbam + '.tmp'), "wb", template=bam)

	for reads in bam.fetch():

		new_tags = reads.tags
		new_tags.append(('HP', haplonum))
		new_tags.append(('CL', clone))
		reads.tags = new_tags
		outbam.write(reads)

	bam.close()
	outbam.close()

	os.remove(os.path.abspath(inbam))
	os.remove(os.path.abspath(inbam + '.bai'))

	os.rename(os.path.abspath(inbam + '.tmp'), os.path.abspath(inbam))




def Simulate(genome, cores, haplotype, chromosome, start, end, label, model_qc, accuracy, coverage, allelic, length, ratio, output, haplonum, clone):


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

	ModifyReadTags(os.path.abspath(output + '/' + label + '.srt.bam'), haplonum, clone)

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))
