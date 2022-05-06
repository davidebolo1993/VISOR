#!/usr/bin/python3 env

#python 3 standard library

import os
import sys
import glob
import re
import subprocess
import math
import multiprocessing
import resource
import random
from datetime import datetime
from collections import defaultdict
from shutil import which,copyfileobj

#additional modules

import pybedtools
import pysam
import pyfaidx
#import mappy as mp #not calling mp.fastx_read 'cause is removing part of the read name (same behaviour of readfq from Heng)

from VISOR import __version__

class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	OUT = ''
	REF = ''
	BED = ''
	SAMPLES = list()

	#other variables

	refall=None
	threads=0
	mmpreset='map-ont'
	fastq=False
	compress=False

	#BadRead

	coverage=0
	regioncoverage=0
	error='nanopore2020'
	quality='nanopore2020'
	length=tuple()
	identity=tuple()
	junk=0
	random=0
	chimeras=0
	glitch=tuple()

	#bulk

	clonefraction=None
	cperc=0
	fperc=0
	ffiles=None
	ffile=None
	sampledir=''
	haplodir=''
	clonenumber=0
	hapnumber=0
	r_number=0
	tag=False


def readfq(fp): # this is a fast generator function

	'''
	Yield FASTQ record
	'''

	read=[]

	for line in fp:

		read.append(line.rstrip())

		if len(read) == 4:

			yield read[0][1:], read[1], read[3]

			read=[]


def Chunks(l,n):

	'''
	Split list in chunks based on number of threads
	'''

	return [l[i:i+n] for i in range(0, len(l), n)]


def atoi(text):

	'''
	Convert text to integers
	'''

	return int(text) if text.isdigit() else text


def natural_keys(text):

	'''
	Natural sort
	'''
	
	return [ atoi(c) for c in re.split(r'(\d+)', text)]


def RTag(sli,c):

	'''
	Add CL/HP-tag to BAM upon request (slows a bit)

	'''

	for s in sli:

		save = pysam.set_verbosity(0)
		bamfilein=pysam.AlignmentFile(s, mode='rb', require_index=False)
		pysam.set_verbosity(save)

		with pysam.AlignmentFile(s+'.tmp', mode='wb', template=bamfilein) as bamfileout:

			for r in bamfilein.fetch(until_eof=True):

				r.set_tag('CL', c.clonenumber, 'i')
				r.set_tag('HP', c.hapnumber, 'i')
				bamfileout.write(r)

		bamfilein.close()
		os.remove(s)
		os.rename(s+'.tmp', s)


def MultiBadRead(processor,c, tmpfa, hapcov, mateh):

	'''
	Write multiple FASTQ files that will be concatenated in the end
	'''

	badread_hap_cmd=['badread', 'simulate', '--reference', tmpfa, '--quantity', str(hapcov/c.threads)+'X', '--length', ','.join(str(x) for x in c.length), '--identity', ','.join(str(x) for x in c.identity), '--start_adapter_seq', '', '--end_adapter_seq', '', '--error_model', c.error, '--qscore_model', c.quality, '--junk_reads', str(c.junk), '--random_reads', str(c.random), '--chimeras', str(c.chimeras), '--glitches', ','.join(str(x) for x in c.glitch)]

	tmpout=os.path.abspath(mateh+'.'+processor)

	with open(tmpout, 'w') as mateout:

		subprocess.call(badread_hap_cmd, stdout=mateout, stderr=open(os.devnull, 'wb')) 


def BulkSim(w,c):

	'''
	Perform bulk simulations and re-align to the un-modified reference
	'''

	hfa=pyfaidx.Fasta(c.ffile)

	if w.chrom not in hfa.keys():

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Chromosome ' + w.chrom + ' not found in ' + c.ffile + '. Skipped simulation')

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Preparing simulation from ' + c.ffile + '. Clone ' + str(c.clonenumber) + '. Haplotype ' + str(c.hapnumber))

		chr_= hfa[w.chrom]
		seq_ = chr_[w.start-1:w.end].seq
		tmpfa=os.path.abspath(c.haplodir + '/' + 'htmp.fa')
		region=w.chrom+'_'+str(w.start)+'_'+str(w.end)

		with open(tmpfa, 'w') as tmpfout: #write temporary fa for sampling reads

			tmpfout.write('>' + region + '\n' + '\n'.join(re.findall('.{1,60}', seq_)) + '\n')

		mateh=os.path.abspath(c.haplodir + '/hr.tmp.fq')

		if float(w[4]) < 100.0:

			tmpref=os.path.abspath(c.haplodir + '/' + 'rtmp.fa')
			seq__=c.refall[w.chrom][w.start-1:w.end].seq
			
			with open(tmpref, 'w') as tmpfout: #write temporary fa for sampling reads

				tmpfout.write('>' + region + '\n' + '\n'.join(re.findall('.{1,60}', seq__)) + '\n')

			#simulate part from reference and part from haplotype
			
			hapcov=c.regioncoverage/100*float(w[4])
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulated coverage for this region will be ' + str(hapcov))

			refcov=c.regioncoverage-hapcov
			print('[' + now + '][Message] Simulated coverage for the corresponding reference region will be ' + str(refcov))
			
			mater=os.path.abspath(c.haplodir + '/rr.tmp.fq')

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulating')
				

			#Simulation with BadRead is great but slow. Implement multi-threading as in https://github.com/rrwick/Badread/wiki/Running-in-parallel. Keep an eye to RAM usage

			processes=[]

			for i in range(c.threads):

				processor='p' + str(i+1)
				p=multiprocessing.Process(target=MultiBadRead, args=(processor,c, tmpfa, hapcov, mateh))
				p.start()
				processes.append(p)

			for p in processes:

				p.join()

			multifiles=glob.glob(os.path.abspath(c.haplodir) + '/hr.*.p*')

			os.remove(tmpfa)

			matehnew=os.path.abspath(c.haplodir + '/hr.fq')

			with open(matehnew, 'w') as outfile:
	
				for fname in multifiles:
		
					with open(fname, 'r') as infile:
		
						for name,seq,qual in readfq(infile):
				
							newname='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fh_' + '|'.join(name.split(' '))
							newname='|'.join(newname.split('|')[:2])
							
							read=[newname, seq, '+', qual]

							outfile.write('\n'.join(read) + '\n')

					os.remove(fname)

			processes=[]

			for i in range(c.threads):

				processor='p' + str(i+1)
				p=multiprocessing.Process(target=MultiBadRead, args=(processor,c, tmpref, refcov, mater))
				p.start()
				processes.append(p)

			for p in processes:

				p.join()

			multifiles=glob.glob(os.path.abspath(c.haplodir) + '/rr.*.p*')

			os.remove(tmpref)

			with open(matehnew, 'a') as outfile:
	
				for fname in multifiles:
		
					with open(fname, 'r') as infile:
		
						for name,seq,qual in readfq(infile):
				
							newname='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fr_' + '|'.join(name.split(' '))
							newname='|'.join(newname.split('|')[:2])
							
							read=[newname, seq, '+', qual]

							outfile.write('\n'.join(read) + '\n')

					os.remove(fname)

		else:

			hapcov=c.regioncoverage
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulated coverage for this region will be ' + str(hapcov))

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulating')

			#Simulation with BadRead is great but slow. Implement multi-threading as in https://github.com/rrwick/Badread/wiki/Running-in-parallel. Keep an eye to RAM usage

			processes=[]

			for i in range(c.threads):

				processor='p' + str(i+1)
				p=multiprocessing.Process(target=MultiBadRead, args=(processor,c, tmpfa, hapcov, mateh))
				p.start()
				processes.append(p)

			for p in processes:

				p.join()

			multifiles=glob.glob(os.path.abspath(c.haplodir) + '/*.p*')

			os.remove(tmpfa)

			matehnew=os.path.abspath(c.haplodir + '/hr.fq')

			with open(matehnew, 'w') as outfile:
	
				for fname in multifiles:
		
					with open(fname, 'r') as infile:
		
						for name,seq,qual in readfq(infile):
				
							newname='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fh_' + '|'.join(name.split(' '))
							newname='|'.join(newname.split('|')[:2])
							
							read=[newname, seq, '+', qual]

							outfile.write('\n'.join(read) + '\n')

					os.remove(fname)

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Mapping simulated reads to the reference genome')

		RGid='nanopore' if c.mmpreset== 'map-ont' else 'pacbio'
		RGstring='@RG\\tID:'+RGid+'\\tSM:bulk'

		BAM=os.path.abspath(c.haplodir+'/'+str(c.r_number)+'.srt.bam')

		sam_cmd = ['minimap2', '-ax', c.mmpreset , '--MD', '--cs', '-Y', '-t', str(c.threads), '-R', RGstring, c.REF, matehnew]
		bam_cmd = ['samtools', 'sort', '-@', str(round(c.threads/2)), '-o', BAM]
		
		p1=subprocess.Popen(sam_cmd, stderr=open(os.devnull, 'wb'), stdout=subprocess.PIPE)
		bout=open(BAM, 'wb')
		p2=subprocess.run(bam_cmd, stdin=p1.stdout, stderr=open(os.devnull, 'wb'), stdout=bout)
		bout.close()

		if c.fastq:

			fastq1=os.path.abspath(c.OUT + '/r.fq')

			if c.compress:
				import gzip
				wfd=gzip.open(fastq1 + '.gz','ab')
				fd=open(matehnew,'rb')
			else:
				wfd=open(fastq1,'a')
				fd=open(matehnew,'r')

			copyfileobj(fd, wfd)

			fd.close()
			wfd.close()

		os.remove(matehnew)


def run(parser,args):

	'''
	Check arguments, run functions
	'''

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] VISOR LASeR v' + __version__)

	#fill container

	c.OUT=os.path.abspath(args.output)
	c.REF=os.path.abspath(args.genome)
	c.BED=os.path.abspath(args.bedfile)
	c.SAMPLES=[os.path.abspath(x) for x in args.sample[0]]
	c.fastq=args.fastq
	c.compress=args.compress

	#main

	if not os.path.exists(c.OUT):

		try:

			os.makedirs(c.OUT)

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Cannot create the output folder')
			sys.exit(1)

	else:

		if not os.access(os.path.abspath(c.OUT),os.W_OK):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Missing write permissions on the output folder')
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(c.OUT)):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] The output folder is not empty: specify another output folder or clean the current one')
			sys.exit(1)


	required=['bedtools', 'minimap2', 'samtools', 'badread']

	for x in required:

		if which(x) is None:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] ' + x +' must be in PATH')
			sys.exit(1)

		else:

			#check version (only for samtools?)
			if x == 'samtools':

				major,minor=subprocess.check_output(['samtools', '--version']).decode('utf-8').split('\n')[0].split('samtools ')[1].split('.')[:2]

				if int(major) < 1 or int(minor) < 9:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Required samtools version >= 1.9')
					sys.exit(1)

			elif x == 'minimap2':

				major,minor=subprocess.check_output(['minimap2', '--version']).decode('utf-8').rstrip().split('-')[0].split('.')

				if int(major) < 2 or int(minor) < 17:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Required minimap2 version >= 2.17')
					sys.exit(1)

	try:

		c.refall=pyfaidx.Fasta(c.REF)

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] Reference file does not exist, is not readable or is not a valid FASTA')
		sys.exit(1)

	try:

		bedfile=pybedtools.BedTool(c.BED)
		bedsrtd=bedfile.sort()

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] BED ' + c.BED + ' does not exist, is not readable or is not a valid BED')
		sys.exit(1)

	#validate all the fields in the input BED

	for j,x in enumerate(bedsrtd):

		if x.chrom not in c.refall.keys(): #not a chromosome in reference file

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Line ' + str(j+1) + ': column 1 (chromosome name) contains an invalid chromosome (not included in the reference provided)')
			sys.exit(1)

		try:

			assert(float(x[3]) <= 100.0) #check value and type same time

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Line ' + str(j+1) + ': column 4 (capture bias percentage) must be a float not greater than 100.0')
			sys.exit(1)

		try:

			float(x[4])

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (purity percentage) must be a float')
			sys.exit(1)

	c.tag=args.tag
	c.threads=args.threads

	if c.threads > multiprocessing.cpu_count():

		c.threads=multiprocessing.cpu_count()-1
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Specified number of cores exceeds the number of cores available. Using all but one')

	#fill c with wgsim parameters for simulations 

	c.coverage=args.coverage
	c.identity=(args.identity_min,args.identity_max,args.identity_stdev)
	c.length=(args.length_mean, args.length_stdev)
	c.error=args.error_model

	if c.error != 'nanopore2020' and c.error != 'nanopore2018' and c.error != 'pacbio2016': #then, should be a file

		c.error=os.path.abspath(args.error_model)

		if not os.path.exists(c.error):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Specified error model is not one of the accepted preset ("nanopore2020", "nanopore2018" or "pacbio2016") neither is a path to a trained error model for badread')
			sys.exit(1)

	c.quality=args.qscore_model

	if c.quality != 'nanopore2020' and c.quality != 'nanopore2018' and c.quality != 'pacbio2016': #then, should be a file

		c.quality=os.path.abspath(args.qscore_model)

		if not os.path.exists(c.quality):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Specified quality model is not one of the accepted preset ("nanopore2020", "nanopore2018" or "pacbio2016") neither is a path to a trained quality model for badread')
			sys.exit(1)

	c.junk=args.junk_reads
	c.random=args.random_reads
	c.chimeras=args.chimera_reads
	c.glitch=(args.glitches_rate, args.glitches_size, args.glitches_skip)

	if args.read_type == 'pacbio':

		c.mmpreset = 'map-pb' #else keep the standard "nanopore"

	if len(c.SAMPLES) > 1:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Preparing for bulk simulations with multiple clones')

		if args.clonefraction is None:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] Multiple clones but their percentages are not specified. Provide coherent clone percentages')	
			sys.exit(1)

		else:		

			c.clonefraction=[float(x) for x in args.clonefraction[0]]

			if len(c.clonefraction) != len(c.SAMPLES):

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] ' + str(len(c.SAMPLES)) + ' clones but ' + str(len(c.clonefraction)) + ' percentages specified')
				sys.exit(1)

			if sum(x for x in c.clonefraction) != 100.0:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] Sum of percentages must equal 100.0')
				sys.exit(1)

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Preparing for bulk simulations with a single clone')
		c.clonefraction=[100.0]

	allbams=[]
	hapdirs=[]
	clonedirs=[]

	for i,x in enumerate(c.SAMPLES):

		clonebams=[]
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Processing clone ' + str(i+1))
		c.sampledir=os.path.abspath(c.OUT + '/clone' + str(i+1))
		clonedirs.append(c.sampledir)
		os.makedirs(c.sampledir)
		c.ffiles=sorted(glob.glob(os.path.abspath(x) + '/*.fa'), key=natural_keys)
		c.cperc=c.clonefraction[i]
		c.fperc=c.cperc/len(c.ffiles) #sample/clone percentage divided equally for each FASTA
		c.clonenumber=i+1

		for k,s in enumerate(c.ffiles):

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Processing haplotype ' + str(k+1))
			c.haplodir=os.path.abspath(c.sampledir + '/h' + str(k+1))
			hapdirs.append(c.haplodir)
			os.makedirs(c.haplodir)
			c.ffile=c.ffiles[k]
			c.hapnumber=k+1

			for w in bedsrtd: #do not use multi-processing on this as minimap2 may require too much memory

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Message] Simulating from region ' + w.chrom + ':' + str(w.start) + '-' + str(w.end))
				c.r_number+=1
				c.regioncoverage=(c.coverage/100*float(w[3]))/100*c.fperc #this takes into account also number of fasta files to split coverage by
				BulkSim(w,c)

			c.r_number=0
			hsam=glob.glob(os.path.abspath(c.haplodir) + '/*.bam')

			if len(hsam) == 0:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Warning] No BAM generated for the current haplotype')

			else:

				if c.tag:

					reg=len(hsam)
					chunk_size=reg/c.threads
					slices=Chunks(hsam,math.ceil(chunk_size))
					processes = []

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Message] Adding CL- and HP-tags to BAM')

					for sli in slices:

						p=multiprocessing.Process(target=RTag, args=(sli,c))
						p.start()
						processes.append(p)

					for p in processes:

						p.join()

			clonebams.extend(hsam)

		if len(clonebams) == 0:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Warning] No BAM generated for the current clone')

		allbams.extend(clonebams)

	if len(allbams) == 0:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] No BAM generated for the current simulation')

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Merging simulated BAM')

		maxfile=round(resource.getrlimit(resource.RLIMIT_NOFILE)[0]/2) #use half of the lower bound
			
		if len(allbams) >= maxfile:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Lots of BAM here! Merging progressively')

			chunk_size=len(allbams)/maxfile
			slices=Chunks(allbams,math.ceil(chunk_size))
			BAML=os.path.abspath(c.OUT+'/bamlist.txt')

			for v,sli in enumerate(slices):

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Message] Round ' + str(v+1) + ' of merging')

				BAMO=os.path.abspath(c.OUT+'/' + str(v+1) + '.sim.srt.bam')

				with open(BAML, 'w') as bamlist:

					for s in sli:

						bamlist.write(s+'\n')

				merge_command = ['samtools', 'merge', '-@', str(c.threads), '-c', '-b', BAML, BAMO]
				subprocess.call(merge_command, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Final round of merging')

			finalbams=glob.glob(os.path.abspath(c.OUT) + '/*.bam')
			BAMO=os.path.abspath(c.OUT+'/sim.srt.bam')
					
			with open(BAML, 'w') as bamlist:

				for s in finalbams:

					bamlist.write(s+'\n')

			merge_command = ['samtools', 'merge', '-@', str(c.threads), '-c', '-b', BAML, BAMO]
			subprocess.call(merge_command, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
			pysam.index(BAMO)

				#remove sub-merged BAM files

			for b in finalbams:

				os.remove(b)				

		else:

			BAML=os.path.abspath(c.OUT+'/bamlist.txt')
			BAMO=os.path.abspath(c.OUT+'/sim.srt.bam')

			with open(BAML, 'w') as bamlist:

				for s in allbams:

					bamlist.write(s+'\n')

			merge_command = ['samtools', 'merge', '-@', str(c.threads), '-c', '-b', BAML, BAMO] #this can ideally be called also with '--write-index' but from 1.10 on
			subprocess.call(merge_command, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
			pysam.index(BAMO)

		os.remove(BAML) #remove BAM list

		for b in allbams: #remove all the initial BAM files

			os.remove(b)

		#clean directories
		#haplotpye-directories

		for d in hapdirs:

			os.rmdir(d)

		#clone-directories

		for d in clonedirs:

			os.rmdir(d)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)
