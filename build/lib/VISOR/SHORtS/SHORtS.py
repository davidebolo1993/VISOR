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
from pywgsim import wgsim
import mappy as mp


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
	fastq=False

	#pywgsim

	coverage=0
	regioncoverage=0
	error=0
	distance=0
	stdev=0
	length=0
	mutation=0
	indels=0
	extindels=0

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

	#strand-seq parameters

	strandseq=False
	sce_bed=None
	sce_bedregion=[]
	noise=0.00
	hapid=''
	cellid=''
	cellnum=0
	cellref=0
	cellhap=0
	singlecellnum=0


def redirect_stdout():

	'''
	Suppress c/c++/fortran stdout, keep python print calls
	'''

	sys.stdout.flush() # <--- important when redirecting to files
	newstdout = os.dup(1)
	devnull = os.open(os.devnull, os.O_WRONLY)
	os.dup2(devnull, 1)
	os.close(devnull)
	sys.stdout = os.fdopen(newstdout, 'w')


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


def MergeAll(c):

	'''
	Merge all the Watson/Crick BAM files in directory
	'''

	TypeW=glob.glob(os.path.abspath(c.haplodir) + '/*.W.srt.bam')


	if len(TypeW) == 0:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] No BAM generated for the current simulation')

	elif len(TypeW) > 1:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Merging Watson BAM for current cell/haplotype')

		BAML=os.path.abspath(c.haplodir+'/bamlist.txt')
		BAMO=os.path.abspath(c.haplodir+'/W.srt.bam')

		with open(BAML, 'w') as bamlist:

			for s in TypeW:

				bamlist.write(s+'\n')

		merge_command = ['samtools', 'merge', '-@', str(c.threads), '-c', '-b', BAML, BAMO]
		subprocess.call(merge_command, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
		pysam.index(BAMO)

		for b in TypeW:

			os.remove(b)

		os.remove(BAML)

	else:

		newbam=os.path.abspath(os.path.dirname(TypeW[0])+ '/W.srt.bam')
		os.rename(TypeW[0],newbam)
		pysam.index(newbam)

	TypeC=glob.glob(os.path.abspath(c.haplodir) + '/*.C.srt.bam')

	if len(TypeC) == 0:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] No BAM generated for the current simulation')

	elif len(TypeC) > 1:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Merging Crick BAM for current cell/haplotype')

		BAML=os.path.abspath(c.haplodir+'/bamlist.txt')
		BAMO=os.path.abspath(c.haplodir+'/C.srt.bam')

		with open(BAML, 'w') as bamlist:

			for s in TypeC:

				bamlist.write(s+'\n')

		merge_command = ['samtools', 'merge', '-@', str(c.threads), '-c', '-b', BAML, BAMO]
		subprocess.call(merge_command, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
		pysam.index(BAMO)

		for b in TypeC:

			os.remove(b)

		os.remove(BAML)

	else:

		newbam=os.path.abspath(os.path.dirname(TypeC[0])+ '/C.srt.bam')
		os.rename(TypeC[0],newbam)
		pysam.index(newbam)


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

		Ns=seq_.count('N') #normalize coverage on Ns
		Nreads=round(((c.regioncoverage*(len(seq_)-Ns))/c.length)/2) #for paired-end sequencing

		mate1h=os.path.abspath(c.haplodir + '/hr1.tmp.fq')
		mate2h=os.path.abspath(c.haplodir + '/hr2.tmp.fq')

		if float(w[4]) < 100.0:

			tmpref=os.path.abspath(c.haplodir + '/' + 'rtmp.fa')
			seq__=c.refall[w.chrom][w.start-1:w.end].seq
			
			with open(tmpref, 'w') as tmpfout: #write temporary fa for sampling reads

				tmpfout.write('>' + region + '\n' + '\n'.join(re.findall('.{1,60}', seq__)) + '\n')

			#simulate part from reference and part from haplotype

			haploreadsN=round(Nreads/100*float(w[4]))

			hapcov=haploreadsN*c.length*2/((w.end-w.start)-Ns)
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulated coverage for this region will be ' + str(hapcov))

			refreadsN=Nreads-haploreadsN
			refcov=refreadsN*c.length*2/((w.end-w.start)-Ns)
			print('[' + now + '][Message] Simulated coverage for the corresponding reference region will be ' + str(refcov))
			
			mate1r=os.path.abspath(c.haplodir + '/rr1.tmp.fq')
			mate2r=os.path.abspath(c.haplodir + '/rr2.tmp.fq')

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulating')

			wgsim.core(r1=mate1h, r2=mate2h, ref=tmpfa, err_rate=c.error, mut_rate=c.mutation, indel_frac=c.indels, indel_ext=c.extindels, N=haploreadsN, dist=c.distance, stdev=c.stdev, size_l=c.length, size_r=c.length,max_n=0.05, is_hap=0, is_fixed=0, seed=0)
			wgsim.core(r1=mate1r, r2=mate2r, ref=tmpref, err_rate=c.error, mut_rate=c.mutation, indel_frac=c.indels, indel_ext=c.extindels, N=refreadsN, dist=c.distance, stdev=c.stdev, size_l=c.length, size_r=c.length,max_n=0.05, is_hap=0, is_fixed=0, seed=0)

			os.remove(tmpfa)
			os.remove(tmpref)

			mate1hnew=os.path.abspath(c.haplodir + '/hr1.fq')
			mate2hnew=os.path.abspath(c.haplodir + '/hr2.fq')

			with open(mate1hnew,'w') as out1, open(mate2hnew,'w') as out2:

				for (name1,seq1,qual1),(name2,seq2,qual2) in zip(mp.fastx_read(mate1h),mp.fastx_read(mate2h)):

					#change name1/name2

					newname1='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fh_' + name1 
					newname2='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fh_' + name2

					read1=[newname1, seq1, '+', qual1]
					read2=[newname2, seq2, '+', qual2]

					out1.write('\n'.join(x for x in read1) + '\n')
					out2.write('\n'.join(x for x in read2) + '\n')

			os.remove(mate1h)
			os.remove(mate2h)

			with open(mate1hnew,'a') as out1, open(mate2hnew,'a') as out2:

				for (name1,seq1,qual1),(name2,seq2,qual2) in zip(mp.fastx_read(mate1r),mp.fastx_read(mate2r)):

					#change name1/name2

					newname1='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fr_' + name1 
					newname2='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fr_' + name2

					read1=[newname1, seq1, '+', qual1]
					read2=[newname2, seq2, '+', qual2]

					out1.write('\n'.join(read1) + '\n')
					out2.write('\n'.join(read2) + '\n')

			os.remove(mate1r)
			os.remove(mate2r)

			#split in chunks for multiprocessing

		else:

			hapcov=Nreads*c.length*2/((w.end-w.start)-Ns)
			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulated coverage for this region will be ' + str(hapcov))

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Simulating')

			wgsim.core(r1=mate1h, r2=mate2h, ref=tmpfa, err_rate=c.error, mut_rate=c.mutation, indel_frac=c.indels, indel_ext=c.extindels, N=Nreads, dist=c.distance, stdev=c.stdev, size_l=c.length, size_r=c.length,max_n=0.05, is_hap=0,  is_fixed=0, seed=0)

			os.remove(tmpfa)

			mate1hnew=os.path.abspath(c.haplodir + '/hr1.fq')
			mate2hnew=os.path.abspath(c.haplodir + '/hr2.fq')

			with open(mate1hnew,'w') as out1, open(mate2hnew,'w') as out2:

				for (name1,seq1,qual1),(name2,seq2,qual2) in zip(mp.fastx_read(mate1h),mp.fastx_read(mate2h)):

					#change name1/name2

					newname1='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fh_' + name1 
					newname2='@c' + str(c.clonenumber) + 'h' + str(c.hapnumber) + 'fh_' + name2

					read1=[newname1, seq1, '+', qual1]
					read2=[newname2, seq2, '+', qual2]

					out1.write('\n'.join(read1) +'\n')
					out2.write('\n'.join(read2) + '\n')

			os.remove(mate1h)
			os.remove(mate2h)

		
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Mapping simulated reads to the reference genome')

		BAM=os.path.abspath(c.haplodir+'/'+str(c.r_number)+'.srt.bam')

		sam_cmd = ['minimap2', '-ax', 'sr', '--MD', '--cs', '-Y', '--sam-hit-only', '-t', str(c.threads), '-R', '@RG\\tID:illumina\\tSM:bulk', c.REF, mate1hnew, mate2hnew]
		bam_cmd = ['samtools', 'sort', '-@', str(round(c.threads/2)), '-o', BAM]
		
		p1=subprocess.Popen(sam_cmd, stderr=open(os.devnull, 'wb'), stdout=subprocess.PIPE)
		bout=open(BAM, 'wb')
		p2=subprocess.run(bam_cmd, stdin=p1.stdout, stderr=open(os.devnull, 'wb'), stdout=bout)
		bout.close()

		if c.fastq:

			fastq1=os.path.abspath(c.OUT + '/r1.fq')
			fastq2=os.path.abspath(c.OUT + '/r2.fq')
			
			with open(fastq1,'a') as wfd:

				with open(mate1hnew,'r') as fd:

					copyfileobj(fd, wfd)

			with open(fastq2,'a') as wfd:

				with open(mate2hnew,'r') as fd:
					
					copyfileobj(fd, wfd)

		os.remove(mate1hnew)
		os.remove(mate2hnew)


def StrandSim(w,c):

	'''
	Perform first part of strand-seq simulations and re-align to the original haplotype
	'''

	hfa=pyfaidx.Fasta(c.ffile)

	if w.chrom not in hfa.keys():

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Chromosome ' + w.chrom + ' not found in ' + c.ffile + '. Skipped simulation')

	else:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Preparing simulation from ' + c.ffile + '. Haplotype ' + str(c.hapnumber))

		chr_= hfa[w.chrom]
		seq_ = chr_[w.start-1:w.end].seq
		tmpfa=os.path.abspath(c.haplodir + '/' + 'htmp.fa')
		region=w.chrom+'_'+str(w.start)+'_'+str(w.end)

		with open(tmpfa, 'w') as tmpfout: #write temporary fa for sampling reads

			tmpfout.write('>' + region + '\n' + '\n'.join(re.findall('.{1,60}', seq_)) + '\n')

		Ns=seq_.count('N') #normalize coverage on Ns
		Nreads=round(((c.regioncoverage*(len(seq_)-Ns))/c.length)/2) #for paired-end sequencing

		mate1h=os.path.abspath(c.haplodir + '/hr1.tmp.fq')
		mate2h=os.path.abspath(c.haplodir + '/hr2.tmp.fq')

		hapcov=Nreads*c.length*2/((w.end-w.start)-Ns)
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Simulated coverage for this region will be ' + str(hapcov))

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Simulating')

		wgsim.core(r1=mate1h, r2=mate2h, ref=tmpfa, err_rate=c.error, mut_rate=c.mutation, indel_frac=c.indels, indel_ext=c.extindels, N=Nreads, dist=c.distance, stdev=c.stdev, size_l=c.length, size_r=c.length,max_n=0.05, is_hap=0,  is_fixed=0, seed=0)

		os.remove(tmpfa)

		mate1hnew=os.path.abspath(c.haplodir + '/hr1.fq')
		mate2hnew=os.path.abspath(c.haplodir + '/hr2.fq')

		with open(mate1hnew,'w') as out1, open(mate2hnew,'w') as out2:

			for (name1,seq1,qual1),(name2,seq2,qual2) in zip(mp.fastx_read(mate1h),mp.fastx_read(mate2h)):

				#change name1/name2

				newname1='@c' + str(c.singlecellnum) + 'h' + str(c.hapnumber) + 'fh_' + name1
				newname2='@c' + str(c.singlecellnum) + 'h' + str(c.hapnumber) + 'fh_' + name2

				read1=[newname1, seq1, '+', qual1]
				read2=[newname2, seq2, '+', qual2]

				out1.write('\n'.join(read1) +'\n')
				out2.write('\n'.join(read2) + '\n')

		os.remove(mate1h)
		os.remove(mate2h)
	
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Mapping simulated reads to the corresponding haplotype')

		BAM=os.path.abspath(c.haplodir+'/'+str(c.r_number)+'.srt.bam')

		sam_cmd = ['minimap2', '-ax', 'sr', '--MD', '--cs', '-Y','--sam-hit-only','-t', str(c.threads),c.ffile, mate1hnew, mate2hnew]
		bam_cmd = ['samtools', 'sort', '-@', str(round(c.threads/2)), '-o', BAM]
		
		p1=subprocess.Popen(sam_cmd, stderr=open(os.devnull, 'wb'), stdout=subprocess.PIPE)
		bout=open(BAM, 'wb')
		p2=subprocess.run(bam_cmd, stdin=p1.stdout, stderr=open(os.devnull, 'wb'), stdout=bout)
		bout.close()

		os.remove(mate1hnew)
		os.remove(mate2hnew)

		#now re-parse BAM file to keep only Watson/Crick reads
		#Watson reads: read1 forward, read2 reverse
		#Crick reads: read2 forward, read1 reverse

		ivf=None

		if len(c.sce_bedregion) != 0:

			sce_string=''

			for s in c.sce_bedregion:

				if s[3] == c.cellid and s[4] == c.hapid:

					sce_string+=s.chrom+'\t'+str(s.start)+'\t'+str(s.end)+'\n'

			if sce_string != '':

				sce_fromscratch=pybedtools.BedTool(sce_string.rstrip(),from_string=True)
				ivf=sce_fromscratch.as_intervalfile() #intervals where to perform SCE events

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Message] Detected one ore more SCE event for current cell/haplotype')


		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Extracting Watson (R1F,R2R) and Crick (R1R,R2F) reads')

		save = pysam.set_verbosity(0)
		bamstrand=pysam.AlignmentFile(BAM, 'rb', require_index=False) #until-eof consumes the bamfile
		pysam.set_verbosity(save)
		Wreads=list(WR(bamstrand,ivf))
		bamstrand.close()


		save = pysam.set_verbosity(0)
		bamstrand=pysam.AlignmentFile(BAM, 'rb', require_index=False) #re-open for second round
		pysam.set_verbosity(save)
		Creads=list(CR(bamstrand,ivf))
		bamstrand.close()

		os.remove(BAM)

		if c.noise > 0:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Message] Adding noise to strands')

			CtoW=random.sample(Creads,round(len(Wreads)/100*c.noise))
			Wreads+=CtoW

			WtoC=random.sample(Wreads,round(len(Creads)/100*c.noise))
			Creads+=WtoC

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Writing Watson and Crick FASTQ')

		w1=os.path.abspath(c.haplodir + '/' + str(c.r_number)+'.w1.fq')
		w2=os.path.abspath(c.haplodir + '/' + str(c.r_number)+'.w2.fq')

		c1=os.path.abspath(c.haplodir + '/' + str(c.r_number)+'.c1.fq')
		c2=os.path.abspath(c.haplodir + '/' + str(c.r_number)+'.c2.fq')

		with open(w1, 'w') as wout1, open(w2, 'w') as wout2:

			for r1,r2 in Wreads:

				if r1.get_tag('OS') == 'W': #this is true W

					read1=['@'+r1.query_name, r1.query_sequence, '+', '2'*c.length]
					read2=['@'+r2.query_name, mp.revcomp(r2.query_sequence), '+', '2'*c.length]

				else: #write to Watson, but is Crick

					read1=['@'+r1.query_name, mp.revcomp(r1.query_sequence), '+', '2'*c.length]
					read2=['@'+r2.query_name, r2.query_sequence, '+', '2'*c.length]

				wout1.write('\n'.join(read1) +'\n')
				wout2.write('\n'.join(read2) +'\n')

		with open(c1, 'w') as cout1, open(c2, 'w') as cout2:

			for r1,r2 in Creads:

				if r1.get_tag('OS') == 'C': #this is true C

					read1=['@'+r1.query_name, mp.revcomp(r1.query_sequence), '+', '2'*c.length]
					read2=['@'+r2.query_name, r2.query_sequence, '+', '2'*c.length]

				else: #write to Crick, but is Watson

					read1=['@'+r1.query_name, r1.query_sequence, '+', '2'*c.length]
					read2=['@'+r2.query_name, mp.revcomp(r2.query_sequence), '+', '2'*c.length]

				cout1.write('\n'.join(read1) +'\n')
				cout2.write('\n'.join(read2) +'\n')

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Mapping Watson and Crick reads to the original reference')

		BAM=os.path.abspath(c.haplodir+'/'+str(c.r_number)+'.W.srt.bam')

		sam_cmd = ['minimap2', '-ax', 'sr', '--MD', '--cs', '-Y', '--sam-hit-only','-t', str(c.threads), '-R', '@RG\\tID:illumina\\tSM:strand', c.REF, w1, w2]
		bam_cmd = ['samtools', 'sort', '-@', str(round(c.threads/2)), '-o', BAM]

		p1=subprocess.Popen(sam_cmd, stderr=open(os.devnull, 'wb'), stdout=subprocess.PIPE)
		bout=open(BAM, 'wb')
		p2=subprocess.run(bam_cmd, stdin=p1.stdout, stderr=open(os.devnull, 'wb'), stdout=bout)
		bout.close()

		os.remove(w1)
		os.remove(w2)

		BAM=os.path.abspath(c.haplodir+'/'+str(c.r_number)+'.C.srt.bam')

		sam_cmd = ['minimap2', '-ax', 'sr', '--MD', '--cs', '-Y', '--sam-hit-only','-t', str(c.threads),'-R', '@RG\\tID:illumina\\tSM:strand',c.REF, c1, c2]
		bam_cmd = ['samtools', 'sort', '-@', str(round(c.threads/2)), '-o', BAM]

		p1=subprocess.Popen(sam_cmd, stderr=open(os.devnull, 'wb'), stdout=subprocess.PIPE)
		bout=open(BAM, 'wb')
		p2=subprocess.run(bam_cmd, stdin=p1.stdout, stderr=open(os.devnull, 'wb'), stdout=bout)
		bout.close()

		os.remove(c1)
		os.remove(c2)


def WR(stranded_bam, ivf):

	'''
	Extract Watson reads from strand-seq BAM. Switch Crick and Watson if hit in ivf
	'''

	WR = defaultdict(lambda: [None, None])

	for read in stranded_bam.fetch(until_eof=True):

		if read.is_proper_pair and not read.is_secondary and not read.is_supplementary:

			if ivf is None: #no need to check for intervals match

				if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse): #read1 forward/ read2 reverse

					read.set_tag('OS', 'W', 'Z') #used for debugging

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

			else: #there is a region to perform W-C switch in

				query=pybedtools.Interval(read.reference_name,read.reference_start,read.reference_end)

				if ivf.any_hits(query) >=1: #yeld Crick as Watson

					if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse): #read2 forward and read1 reverse

						read.set_tag('OS', 'C', 'Z') #used for debugging

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

				else: #Classic Watson Reads

					if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse): #read1 forward and read2 reverse

						read.set_tag('OS', 'W', 'Z') #used for debugging

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

def CR(stranded_bam, ivf):

	'''
	Extract Crick reads from strand-seq BAM. Switch Watson and Crick if hit in ivf
	'''

	CR = defaultdict(lambda: [None, None])

	for read in stranded_bam.fetch(until_eof=True):

		if read.is_proper_pair and not read.is_secondary and not read.is_supplementary:

			if ivf is None: #no need to check for intervals match

				if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse): #read2 forward and read1 reverse

					read.set_tag('OS', 'C', 'Z') #used for debugging

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

			else: #there is a region to perform W-C switch in

				query=pybedtools.Interval(read.reference_name,read.reference_start,read.reference_end)

				if ivf.any_hits(query) >=1: #yeld Watson as Crick

					if (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):

						read.set_tag('OS', 'W', 'Z') #used for debugging

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

				else: #Classic Crick Read

					if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse): #read2 forward and read1 reverse

						read.set_tag('OS', 'C', 'Z') #used for debugging

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


def run(parser,args):

	'''
	Check arguments, run functions
	'''

	redirect_stdout()# block pywgsim stdout

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] VISOR SHORtS v1.1')

	#fill container

	c.OUT=os.path.abspath(args.output)
	c.REF=os.path.abspath(args.genome)
	c.BED=os.path.abspath(args.bedfile)
	c.SAMPLES=[os.path.abspath(x) for x in args.sample[0]]
	c.strandseq=args.strandseq
	c.fastq=args.fastq

	if c.strandseq and len(c.SAMPLES) > 1:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] When performing strand-seq simulations, only 1 sample must be given as input')
		sys.exit(1)

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


	required=['bedtools', 'minimap2', 'samtools']

	for x in required:

		if which(x) is None:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] ' + x +' must be in PATH')
			sys.exit(1)

		else:

			#check version (only for samtools?)
			if x == 'samtools':

				major,minor=subprocess.check_output(['samtools', '--version']).decode('utf-8').split('\n')[0].split('samtools ')[1].split('.')

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
	c.noise=args.noise

	if c.threads > multiprocessing.cpu_count():

		c.threads=multiprocessing.cpu_count()-1
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Warning] Specified number of cores exceeds the number of cores available. Using all but one')

	#fill c with wgsim parameters for simulations 

	c.coverage=args.coverage
	c.error=args.error
	c.distance=args.distance
	c.stdev=args.stdev
	c.length=args.length
	c.mutation=args.mutation
	c.indels=args.indels
	c.extindels=args.extindels

	if not c.strandseq:

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

	else:

		#debug errors for wrong arguments as input

		c.cellnum=args.cells
		c.cellref=round(c.cellnum/100*args.refcells) #round to integer
		c.cellhap=c.cellnum-c.cellref

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] ' +str(c.cellref) + ' cells will be simulated from reference')
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] ' +str(c.cellhap) + ' cells will be simulated from sample')

		if args.sce is not None:

			c.sce_bed=os.path.abspath(args.sce)

			try:

				sce_bedfile=pybedtools.BedTool(c.sce_bed)
				sce_bedsrtd=sce_bedfile.sort()

			except:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] BED ' + c.sce_bed + ' does not exist, is not readable or is not a valid BED')
				sys.exit(1)

			for s in sce_bedsrtd:

				c.sce_bedregion.append(s)

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Preparing for strand-seq simulation with ' + str(c.cellnum) + ' single cells')

		if c.cellref > 0:

			counterref=0

			for cr in range(c.cellref):

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Message] Processing reference cell ' + str(cr+1))
				c.singlecellnum+=1
				c.sampledir=os.path.abspath(c.OUT + '/reference_cell' + str(cr+1))
				os.makedirs(c.sampledir)
				c.ffiles=[c.REF, c.REF]
				c.cperc=100.0 #all from diploid reference
				c.fperc=c.cperc/len(c.ffiles) #sample/clone percentage divided equally for each FASTA
				c.cellid='reference_cell'+str(cr+1) #for SCE events, if SCE BED provided

				for k,s in enumerate(c.ffiles):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Message] Processing haplotype ' + str(k+1))
					c.haplodir=os.path.abspath(c.sampledir + '/h' + str(k+1))
					os.makedirs(c.haplodir)
					c.ffile=c.ffiles[k]
					c.hapid='haplotype' + str(k+1) #for SCE events, if SCE BED provided
					c.hapnumber=k+1

					for w in bedsrtd:

						now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
						print('[' + now + '][Message] Simulating from region ' + w.chrom + ':' + str(w.start) + '-' + str(w.end))
						c.r_number+=1
						c.regioncoverage=c.coverage/100*c.fperc #this takes into account just number of fasta files to split coverage by
						StrandSim(w,c)

					#merge all watson/crick in haplotype. Assume they cannot exceed resource limit for a single haplotype of a single cell
					MergeAll(c)
					c.r_number=0

				counterref+=1

				if counterref == c.cellref: #break on occurrence

					break

		if c.cellhap > 0:

			counterhap=0

			for cr in range(c.cellhap):

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Message] Processing sample cell ' + str(cr+1))
				c.singlecellnum+=1
				c.sampledir=os.path.abspath(c.OUT + '/sample_cell' + str(cr+1))
				os.makedirs(c.sampledir)
				c.ffiles=sorted(glob.glob(os.path.abspath(c.SAMPLES[0]) + '/*.fa'), key=natural_keys)
				c.cperc=100.0 #all from sample (can have any ploidy)
				c.fperc=c.cperc/len(c.ffiles) #sample/clone percentage divided equally for each FASTA
				c.cellid='sample_cell'+str(cr+1) #for SCE events, if SCE BED provided

				for k,s in enumerate(c.ffiles):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Message] Processing haplotype ' + str(k+1))
					c.haplodir=os.path.abspath(c.sampledir + '/h' + str(k+1))
					os.makedirs(c.haplodir)
					c.ffile=c.ffiles[k]
					c.hapid='haplotype' + str(k+1) #for SCE events, if SCE BED provided
					c.hapnumber=k+1

					for w in bedsrtd:

						now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
						print('[' + now + '][Message] Simulating from region ' + w.chrom + ':' + str(w.start) + '-' + str(w.end))
						c.r_number+=1
						c.regioncoverage=c.coverage/100*c.fperc #this takes into account just number of fasta files to split coverage by
						StrandSim(w,c)

					MergeAll(c)
					c.r_number=0

				counterhap+=1

				if counterhap == c.cellhap: #break on occurrence

					break

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)