#!/usr/bin/python env

#python 3 standard library

import os
import sys
import glob
import logging
from shutil import which
import subprocess
import random
from collections import defaultdict

#additional modules

import pybedtools
import pysam
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


	logging.basicConfig(filename=os.path.abspath(args.output + '/VISOR_SHORtS.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


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
			BWA_Index(os.path.abspath(args.genome))

		except:

			logging.error('It was not possible to generate bwa index for reference genome')
			sys.exit(1)


	#validate .bed with regions to simulate


	bed = pybedtools.BedTool(os.path.abspath(args.bedfile))

	
	try:

		srtbed = bed.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		logging.error('Incorrect .bed format for -bed/--bedfile')
		sys.exit(1)



	#validate SCEBED if single-strand


	if args.type=='single-strand':

		if args.scebedfile is not None:

			scebed=pybedtools.BedTool(os.path.abspath(args.scebedfile))

			try:

				srtscebed = scebed.sort()

			except:

				logging.error('Incorrect .bed format for -scebed/--scebedfile')
				sys.exit(1)


		else:

			srtscebed = None


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

					logging.error('Each fraction percentage in -cf/--clonefraction must be a float')
					sys.exit(1)


			totalfraction = sum(map(float,fractions))

			if totalfraction > 100:

				logging.error('Sum of fractions percentages in -cf/--clonefraction cannot exceed 100.0')
				sys.exit(1)



	fa=pyfaidx.Fasta(os.path.abspath(args.genome))
	generate=os.path.abspath(os.path.dirname(__file__) + '/generate.sh')
	classic_chrs = fa.keys() #allowed chromosomes
	tag=args.noaddtag
	logging.info('Running simulations')

	if len(inputs) == 1: #just one folder, use a classic doulbe-strand or single-strand simulations

		logging.info('Single input for -s/--sample')

		if args.type == 'double-strand':

			logging.info('Double-strand data')

		else:

			logging.info('Single-strand data')

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

					logging.error('Cannot convert ' + str(entries[3]) + ' to float number in .bed file. Capture bias must be a float')
					sys.exit(1)


				try:

					allelic=float(entries[4])

				except:

					logging.error('Cannot convert ' + str(entries[4]) + ' to float number in .bed file. Sample fraction must be a float percentage')
					sys.exit(1)


				if alellic > 100:

					logging.error('Purity value cannot exceed 100.0 in .bed file.')
					sys.exit(1)

				try:

					if args.type == 'double-strand':

						ClassicSimulate(tag, os.path.abspath(args.genome), args.threads, os.path.abspath(fasta), str(entries[0]), int(entries[1]), int(entries[2]), str(counter), allelic, args.error, (args.coverage / 100 * float(entries[3]))/len(fastas), args.length, args.indels, args.probability, args.insertsize, args.standardev, os.path.abspath(args.output + '/' + str(folder)), folder+1, 1)

					else: #generate single-stranded data


						if not os.path.exists(os.path.abspath(fasta + '.sa')):

							try:

								logging.info('Creating bwa index for ' + os.path.abspath(fasta))
								BWA_Index(os.path.abspath(fasta))

							except:

								logging.error('It was not possible to generate bwa index for ' + os.path.abspath(fasta))
								sys.exit(1)



						haploname = os.path.basename(os.path.abspath(fasta)).split('.')[0] #this is important only if scebed is given

						SSSimulate(args.threads, os.path.abspath(fasta), str(entries[0]), int(entries[1]), int(entries[2]), args.error, (args.coverage / 100 * float(entries[3]))/len(fastas), args.length, args.indels, args.probability, args.insertsize, args.standardev, os.path.abspath(args.output + '/' + str(folder)))
						SingleStrand(haploname, str(entries[0]), generate, os.path.abspath(args.genome), args.threads, os.path.abspath(args.output + '/' + str(folder) + '/region.tmp.srt.bam'), str(counter), args.noise, os.path.abspath(args.output + '/' + str(folder)), srtscebed)

				except:

					
					logging.exception('Something went wrong during simulations for ' + os.path.abspath(fasta) + '. Log is below.')


		
		subdirs=[os.path.join(os.path.abspath(args.output), o) for o in os.listdir(os.path.abspath(args.output)) if os.path.isdir(os.path.join(os.path.abspath(args.output),o))]

		if args.type == 'double-strand':

			logging.info('Merging double-strand data')

			bams = [y for x in os.walk(os.path.abspath(args.output)) for y in glob.glob(os.path.join(x[0], '*.srt.bam'))]

			with open(os.path.abspath(args.output + '/bamstomerge.txt'), 'w') as bamstomerge:

				for bam in bams:

					bamstomerge.write(bam + '\n')

			subprocess.call(['samtools', 'merge', '-@', str(args.threads-1), '-b', os.path.abspath(args.output + '/bamstomerge.txt'), os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')], stderr=open(os.devnull, 'wb'))
			subprocess.call(['samtools', 'index', os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')],stderr=open(os.devnull, 'wb'))

			os.remove(os.path.abspath(args.output + '/bamstomerge.txt'))

			for bam in bams:

				os.remove(bam)
				os.remove(bam + '.bai')

			for dirs in subdirs: #now they are empty and vcan be removed safely

				os.rmdir(dirs)

		else: # args.type is single-strand, merge watsons from same haplos and cricks from same haplos. Additional script to merge between haplos.

			logging.info('Merging watson and crick .bam files for each haplotype')

			for dirs in subdirs:

				watsons=[]
				cricks=[]

				watsons.extend(glob.glob(os.path.abspath(dirs) + '/*watson.srt.bam'))
				cricks.extend(glob.glob(os.path.abspath(dirs) + '/*crick.srt.bam'))

				with open(os.path.abspath(dirs + '/watsonstomerge.txt'),'w') as watsonstomerge:

					for wats in watsons:

						watsonstomerge.write(wats + '\n')

				subprocess.call(['samtools', 'merge', '-@', str(args.threads-1), '-b', os.path.abspath(dirs + '/watsonstomerge.txt'), os.path.abspath(dirs + '/' + args.identifier + '.watson.srt.bam')], stderr=open(os.devnull, 'wb'))
				subprocess.call(['samtools', 'index', os.path.abspath(dirs + '/' + args.identifier + '.watson.srt.bam')],stderr=open(os.devnull, 'wb'))

				os.remove(os.path.abspath(dirs + '/watsonstomerge.txt'))

				for wats in watsons:

					os.remove(wats)
					os.remove(wats + '.bai')


				with open(os.path.abspath(dirs + '/crickstomerge.txt'),'w') as crickstomerge:

					for cri in cricks:

						crickstomerge.write(cri + '\n')

				subprocess.call(['samtools', 'merge', '-@', str(args.threads-1), '-b', os.path.abspath(dirs + '/crickstomerge.txt'), os.path.abspath(dirs + '/' + args.identifier + '.crick.srt.bam')], stderr=open(os.devnull, 'wb'))
				subprocess.call(['samtools', 'index', os.path.abspath(dirs + '/' + args.identifier + '.crick.srt.bam')],stderr=open(os.devnull, 'wb'))

				os.remove(os.path.abspath(dirs + '/crickstomerge.txt'))

				for cri in cricks:

					os.remove(cri)
					os.remove(cri + '.bai') #there will be an additional script to generate WC .bam files


	else: # simulate subclones with double-stranded data

		logging.info('Multiple inputs for -s/--sample. Assuming each input is a subclone')
		logging.info('Double-strand data')

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

						logging.error('Cannot convert ' + str(entries[3]) + ' to float number in .bed file. Capture bias must be a float')
						sys.exit(1)


					try:

						ClassicSimulate(tag,os.path.abspath(args.genome), args.threads, os.path.abspath(subfasta), str(entries[0]), int(entries[1]), int(entries[2]), str(counter),100.0, args.error, ((args.coverage / 100 * float(entries[3]))/100)*eachhaplofraction, args.length, args.indels, args.probability, args.insertsize, args.standardev, os.path.abspath(args.output + '/' + str(fract) + '/' + str(folder)), folder+1, fract+1)

					except:

						logging.exception('Something went wrong during simulations for ' + os.path.abspath(subfasta) + '. Log is below.')

		subdirs=[os.path.join(os.path.abspath(args.output), o) for o in os.listdir(os.path.abspath(args.output)) if os.path.isdir(os.path.join(os.path.abspath(args.output),o))]
		bams = [y for x in os.walk(os.path.abspath(args.output)) for y in glob.glob(os.path.join(x[0], '*.srt.bam'))]


		with open(os.path.abspath(args.output + '/bamstomerge.txt'), 'w') as bamstomerge:

			for bam in bams:

				bamstomerge.write(bam + '\n')


		subprocess.call(['samtools', 'merge', '-@', str(args.threads-1), '-b', os.path.abspath(args.output + '/bamstomerge.txt'), os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')], stderr=open(os.devnull, 'wb'))
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




def BWA_Index(fasta):

	subprocess.call(['bwa', 'index', os.path.abspath(fasta)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))



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



def ClassicSimulate(tag,genome, cores, haplotype, chromosome, start, end, label, allelic, error, coverage, length, indels, probability, insertsize, standarddev, output, haplonum, clone):

	#prepare region

	fa=pyfaidx.Fasta(os.path.abspath(haplotype))

	if chromosome not in fa.keys():

		return

	chr_= fa[chromosome]
	seq = chr_[:len(chr_)].seq

	with open(os.path.abspath(output + '/region.tmp.fa'), 'w') as regionout:

		subprocess.call(['samtools', 'faidx', haplotype, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))

	if len(seq) < end-start:

		logging.warning(str(chromosome) + ' in haplotype ' + os.path.abspath(haplotype) + ' is shorter than region to simulate.')
		numreads= round((coverage*len(seq)) / length)/2 #calculate chosen coverage and divide by 2 beacuse they are pairs

	else:

		numreads=round((coverage*(end-start)) / length)/2 

	
	if not allelic == 100:

		with open(os.path.abspath(output + '/reference.region.tmp.fa'), 'w') as regionout:

			subprocess.call(['samtools', 'faidx', genome, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))

		numreads1 = round((numreads/100)*allelic)		
		numreads2 = numreads - numreads1

		subprocess.call(['wgsim', '-e', str(error), '-d', str(insertsize), '-s', str(standarddev), '-N', str(numreads1), '-1', str(length), '-2', str(length), '-R', str(indels), '-X', str(probability), os.path.abspath(output + '/region.tmp.fa'), os.path.abspath(output + '/region.region.1.fq'), os.path.abspath(output + '/region.region.2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))
		subprocess.call(['wgsim', '-e', str(error), '-d', str(insertsize), '-s', str(standarddev), '-N', str(numreads2), '-1', str(length), '-2', str(length), '-R', str(indels), '-X', str(probability), os.path.abspath(output + '/reference.region.tmp.fa'), os.path.abspath(output + '/reference.region.1.fq'), os.path.abspath(output + '/reference.region.2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

		with open (os.path.abspath(output + '/region.1.fq'), 'w') as regionout:

			subprocess.call(['cat', os.path.abspath(output + '/region.region.1.fq'), os.path.abspath(output + '/reference.region.1.fq')], stdout=regionout, stderr=open(os.devnull, 'wb'))


		os.remove(os.path.abspath(output + '/region.region.1.fq'))
		os.remove(os.path.abspath(output + '/reference.region.1.fq'))


		with open (os.path.abspath(output + '/region.2.fq'), 'w') as regionout:

			subprocess.call(['cat', os.path.abspath(output + '/region.region.2.fq'), os.path.abspath(output + '/reference.region.2.fq')], stdout=regionout, stderr=open(os.devnull, 'wb'))


		os.remove(os.path.abspath(output + '/region.region.2.fq'))
		os.remove(os.path.abspath(output + '/reference.region.2.fq'))

		os.remove(os.path.abspath(output + '/reference.region.tmp.fa'))


	else:

		subprocess.call(['wgsim', '-e', str(error), '-d', str(insertsize), '-s', str(standarddev),'-N', str(numreads), '-1', str(length), '-2', str(length), '-R', str(indels), '-X', str(probability), os.path.abspath(output + '/region.tmp.fa'), os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

	
	os.remove(os.path.abspath(output + '/region.tmp.fa'))

	with open(os.path.abspath(output + '/region.tmp.sam'), 'w') as samout:

		subprocess.call(['bwa', 'mem', '-t', str(cores), genome, os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stdout=samout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.1.fq'))
	os.remove(os.path.abspath(output + '/region.2.fq'))

	with open(os.path.abspath(output + '/region.tmp.bam'), 'w') as bamout:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/region.tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.sam'))

	with open(os.path.abspath(output + '/' + label + '.srt.bam'), 'w') as srtbamout:

		subprocess.call(['samtools', 'sort', '-@', str(cores-1), os.path.abspath(output + '/region.tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.bam'))

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))
	
	if tag:

		ModifyReadTags(os.path.abspath(output + '/' + label + '.srt.bam'), haplonum, clone)

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.srt.bam')],stderr=open(os.devnull, 'wb'))



def SSSimulate(cores, haplotype, chromosome, start, end, error, coverage, length, indels, probability, insertsize, standarddev, output):

	#prepare region

	fa=pyfaidx.Fasta(os.path.abspath(haplotype))

	if chromosome not in fa.keys():

		return

	chr_= fa[chromosome]
	seq = chr_[:len(chr_)].seq

	with open(os.path.abspath(output + '/region.tmp.fa'), 'w') as regionout:

		subprocess.call(['samtools', 'faidx', haplotype, chromosome + ':' + str(start) +  '-' +str(end)], stdout=regionout, stderr=open(os.devnull, 'wb'))

	if len(seq) < end-start:

		logging.warning(str(chromosome) + ' in haplotype ' + os.path.abspath(haplotype) + ' is shorter than region to simulate.')
		numreads= round((coverage*len(seq)) / length)/2 #calculate chosen coverage and divide by 2 beacuse they are pairs

	else:

		numreads= round((coverage*(end-start)) / length)/2 
	

	#simulate reads

	subprocess.call(['wgsim', '-e', str(error), '-N', str(numreads), '-1', str(length), '-2', str(length), '-R', str(indels), '-X', str(probability), os.path.abspath(output + '/region.tmp.fa'), os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.fa'))

	#align to modified reference

	with open(os.path.abspath(output + '/region.tmp.sam'), 'w') as samout:

		subprocess.call(['bwa', 'mem', '-t', str(cores), haplotype, os.path.abspath(output + '/region.1.fq'), os.path.abspath(output + '/region.2.fq')], stdout=samout, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/region.tmp.bam'), 'w') as bamout:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/region.tmp.sam')], stdout=bamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.sam'))

	with open(os.path.abspath(output + '/region.tmp.srt.bam'), 'w') as srtbamout:

		subprocess.call(['samtools', 'sort', '-@', str(cores-1), os.path.abspath(output + '/region.tmp.bam')], stdout=srtbamout, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/region.tmp.bam'))

	subprocess.call(['samtools', 'index', os.path.abspath(output + '/region.tmp.srt.bam')],stderr=open(os.devnull, 'wb'))




def SingleStrand(haploname, chromosome, generate, genome, cores, bamfilein, label, noisefraction, output, scebed):

	bam = pysam.AlignmentFile(bamfilein, "rb")

	watslist=list(watson_orientation(bam))	
	cricklist=list(crick_orientation(bam))

	if scebed is not None:

		for entries in scebed:

			if str(entries[3]) == haploname: #perform SCE only on wanted haplotypes

				if str(entries[0]) == chromosome:

					watsinregion = list(watson_orientation_inregion(bam, chromosome, int(entries[1]), int(entries[2])))
					crickinregion = list(crick_orientation_inregion(bam, chromosome, int(entries[1]), int(entries[2])))

					watslist=list(set(watslist)-set(watsinregion)) + crickinregion
					cricklist = list(set(cricklist)-set(crickinregion)) + watsinregion


	with open(os.path.abspath(output + '/watsonreads.txt'), 'w') as watsonreads:

		for read1,read2 in watslist:
				
			watsonreads.write(read1.query_name + '\n')
			watsonreads.write(read2.query_name + '\n')	



	with open(os.path.abspath(output + '/crickreads.txt'), 'w') as crickreads:

		for read1,read2 in cricklist:
				
			crickreads.write(read1.query_name + '\n')
			crickreads.write(read2.query_name + '\n')	



	if noisefraction > 0:

		discordant_pair_watson= round((len(watslist)*noisefraction)/100)
		discordant_pair_crick= round((len(cricklist)*noisefraction)/100)

		sample_crick=random.sample(cricklist,discordant_pair_watson)

		with open(os.path.abspath(output + '/watsonreads.txt'), 'a') as watsonreads:

			for read1,read2 in sample_crick:

				watsonreads.write(read1.query_name + '\n')
				watsonreads.write(read2.query_name + '\n')	


		sample_watson=random.sample(watslist,discordant_pair_crick)

		with open(os.path.abspath(output + '/crickreads.txt'), 'a') as crickreads:

			for read1,read2 in sample_watson:

				crickreads.write(read1.query_name + '\n')
				crickreads.write(read2.query_name + '\n')	



	bam.close()

	os.remove(bamfilein)
	os.remove(bamfilein + '.bai')



	subprocess.call(['bash', generate, os.path.abspath(output)])

	os.remove(os.path.abspath(output + '/watsonreads.txt'))
	os.remove(os.path.abspath(output + '/crickreads.txt'))


	with open(os.path.abspath(output + '/watson.tmp.sam'), 'w') as watsonsam:

		subprocess.call(['bwa', 'mem', '-t', str(cores), genome, os.path.abspath(output + '/watson.1.fq'), os.path.abspath(output + '/watson.2.fq')], stdout=watsonsam, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/crick.tmp.sam'), 'w') as cricksam:

		subprocess.call(['bwa', 'mem', '-t', str(cores), genome, os.path.abspath(output + '/crick.1.fq'), os.path.abspath(output + '/crick.2.fq')], stdout=cricksam, stderr=open(os.devnull, 'wb'))


	os.remove(os.path.abspath(output + '/watson.1.fq'))
	os.remove(os.path.abspath(output + '/watson.2.fq'))

	os.remove(os.path.abspath(output + '/crick.1.fq'))
	os.remove(os.path.abspath(output + '/crick.2.fq'))


	with open(os.path.abspath(output + '/watson.tmp.bam'), 'w') as watsonbam:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/watson.tmp.sam')], stdout=watsonbam, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/crick.tmp.bam'), 'w') as crickbam:

		subprocess.call(['samtools', 'view', '-b', os.path.abspath(output + '/crick.tmp.sam')], stdout=crickbam, stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(output + '/watson.tmp.sam'))
	os.remove(os.path.abspath(output + '/crick.tmp.sam'))


	with open(os.path.abspath(output + '/' + label + '.watson.srt.bam'), 'w') as watsonsort:

		subprocess.call(['samtools', 'sort', '-@', str(cores-1), os.path.abspath(output + '/watson.tmp.bam')], stdout=watsonsort, stderr=open(os.devnull, 'wb'))


	with open(os.path.abspath(output + '/' + label + '.crick.srt.bam'), 'w') as cricksort:

		subprocess.call(['samtools', 'sort', '-@', str(cores-1), os.path.abspath(output + '/crick.tmp.bam')], stdout=cricksort, stderr=open(os.devnull, 'wb'))


	os.remove(os.path.abspath(output + '/watson.tmp.bam'))
	os.remove(os.path.abspath(output + '/crick.tmp.bam'))


	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.watson.srt.bam')],stderr=open(os.devnull, 'wb'))
	subprocess.call(['samtools', 'index', os.path.abspath(output + '/' + label + '.crick.srt.bam')],stderr=open(os.devnull, 'wb'))



def watson_orientation_inregion(bam, chromosome, start, end):

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




def crick_orientation_inregion(bam, chromosome, start, end):

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
