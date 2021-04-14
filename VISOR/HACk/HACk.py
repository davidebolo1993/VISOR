#!/usr/bin/python3 env

#python 3 standard library

import os
import sys
import random
import re
import bisect
from datetime import datetime
from operator import itemgetter
from shutil import which

#additional modules

import pybedtools
import pyfaidx

class c():

	'''
	Container. This stores argparser parameters. Used to pass multiple parameters at once.
	'''

	OUT = ''
	REF = ''
	BED = list()
	store = False

class Overlap():

	'''
	Exclude overlapping variants
	'''

	def __init__(self):
		
		self._intervals = []

	def intervals(self):
		
		return self._intervals

	def put(self, interval):
		
		istart, iend, itype= interval
		# Ignoring intervals that start after the window.                                       
		i = bisect.bisect_right(self._intervals, (iend, sys.maxsize))

		# Look at remaining intervals to find overlap.                                          
		for start, end, typ in self._intervals[:i]:
			
			if end >= istart:
				
				return False
		
		bisect.insort(self._intervals, interval)
		
		return True

def RandomI(nucs):

	'''
	Insert random nucleotide in nucleotide sequence
	'''

	index = random.randint(0, len(nucs)-1)
	nucs = nucs[:index] + random.choice(['A','T','C','G']) + nucs[index:]

	return nucs


def RandomD(nucs):

	'''
	Delete random nucleotide in nucleotide sequence
	'''

	index = random.randint(0, len(nucs)-1)
	nucs = nucs[:index] + nucs[index+1:]

	return nucs


def RandomX(nucs):

	'''
	Substitute random nucleotide in nucleotide sequence
	'''

	allnucs=['A','T','C','G']
	index = random.randint(0, len(nucs)-1)
	subchar=nucs[index]
	allnucs.remove(subchar)
	nucs=nucs[:index]+random.choice(allnucs)+nucs[index+1:] 

	return nucs


def write_chrom(chrs,pyseq_seq,hapfout):

	'''
	Write FASTA chromosome to file
	'''

	with open (hapfout, 'a') as faout:

		faout.write('>' + chrs + '\n' + '\n'.join(re.findall('.{1,60}', pyseq_seq)) + '\n') #write 60-chars FASTA


def HapMaker(pyref,pychroms,hapdict,hapfout):

	'''
	Create variant strings and insert them into FASTA haplotype
	'''

	for chrs in pychroms:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Processing chromosome ' + chrs)

		chrom=pyref[chrs]
		pyseq=chrom[:len(chrom)] #do not use .seq yet
		final_seq=''

		if chrs not in hapdict.keys(): #this means that it has not to be modified

			final_seq+=pyseq.seq

		else:

			#remove duplicates and overlapping features, if any

			alts=sorted(hapdict[chrs], key=itemgetter(0,1))
			ranges=[(x[0],x[1],x[2]) for x in alts]
			ov=Overlap()

			for el in ranges:

				ov.put(el)

			if len(ov.intervals()) < len(ranges):

				diff=len(ranges)-len(ov.intervals())
				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Warning] ' + str(diff) + ' variants will be skipped for current chromosome as they overlap others (only the first occuring among overlapping variants will be kept)')

			altsfltrd=[x for x in alts if (x[0],x[1],x[2]) in ov.intervals()] #exclude also different types

			i=0

			while i < len(altsfltrd):

				s,e,ty,inf,sb=altsfltrd[i]
				
				if s==0: #but this is unlikely to happen, no variants at the start/end of a chromosome 

					s+=1

				if i == 0: #write until first variant

					final_seq+=pyseq[:s-1].seq

				#now arrange variants

				if ty == 'inversion': #create reverse complement, as it has not been created previously

					final_seq+=pyseq[s-1:e].reverse.complement.seq+sb

				elif ty == 'deletion': #add nothing

					if e-s == 1 and inf == '1bp':

						final_seq+=pyseq[e-1].seq #ignore sb for 1 bp DEL

					else:

						final_seq+=sb

				elif ty == 'insertion': #for translocations, insertions are already in the correct f/r orientation

					final_seq+=pyseq[s-1:e].seq+inf+sb

				elif ty == 'deletion-insertion':

					final_seq+=inf+sb

				elif ty == 'tandem duplication':

					final_seq+=pyseq[s-1:e].seq*inf+sb

				elif ty == 'inverted tandem duplication':

					final_seq+=pyseq[s-1:e].seq+pyseq[s-1:e].reverse.complement.seq*(inf-1)+sb

				elif ty == 'SNP':

					final_seq+=pyseq[s-1:e-1].seq+inf #ignore sb for SNP

				elif ty == 'MNP': #ignore sb for MNP

					final_seq+=inf

				elif ty == 'perfect tandem repetition': #ignore sb for PTR

					motif,number=inf.split(':')
					final_seq+=pyseq[s-1:e].seq+motif*int(number)

				elif ty == 'approximate tandem repetition': #ignore sb for ATR

					motif,number,error=inf.split(':')
					app_rep=motif*int(number)
					allalts=['I','D','X']

					for a in range(int(error)):

						#pick random alt
						atype=random.choice(allalts)

						if atype == 'I':

							app_rep=RandomI(app_rep)

						elif atype == 'D':

							app_rep=RandomD(app_rep)

						else:

							app_rep=RandomX(app_rep)

					final_seq+=pyseq[s-1:e].seq+app_rep

				elif ty == 'tandem repeat expansion': #ignore sb for TRE

					motif,number=inf.split(':')
					final_seq+=pyseq[s-1:e].seq+motif*int(number)

				elif ty == 'tandem repeat contraction': #ignore sb for TRC

					motif,number=inf.split(':')
					rep=pyseq[s:e].seq
					index=len(motif)*int(number)
					final_seq+=pyseq[s-1].seq+rep[index:]

				if i < len(altsfltrd)-1: #till next start

					final_seq+=pyseq[e:altsfltrd[i+1][0]-1].seq

				elif i == len(altsfltrd)-1: #till chromosome end

					final_seq+=pyseq[e:].seq

				i+=1

		write_chrom(chrs,final_seq,hapfout)


def run(parser,args):

	'''
	Check arguments, run functions
	'''

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] VISOR HACk v1.1')

	#fill container

	c.OUT=os.path.abspath(args.output)
	c.REF=os.path.abspath(args.genome)
	c.BED=[os.path.abspath(x) for x in args.bedfile[0]]
	c.store=args.vcf #but not yet used. Just for future reference
	
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


	if which('bedtools') is None:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] bedtools must be in PATH')
		sys.exit(1)

	try:

		ref=pyfaidx.Fasta(c.REF)
		chrs=ref.keys()

	except:

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Error] Reference file does not exist, is not readable or is not a valid FASTA')
		sys.exit(1)


	#accepted variants
	possible_variants = ['SNP', 'MNP', 'inversion', 'deletion', 'insertion', 'tandem duplication', 'inverted tandem duplication', 'perfect tandem repetition', 'approximate tandem repetition', 'tandem repeat expansion', 'tandem repeat contraction', 'reciprocal translocation', 'translocation cut-paste', 'translocation copy-paste', 'interspersed duplication']
	valid_dna = 'ACGT'
	haplopattern=re.compile("^h[0-9]+$") #allowed haplotypes for inter-haplotype variants (h1,h2,...)

	#this will contain variants for each haplotype
	d=dict()
	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Validating variants in BED')

	for i,bed in enumerate(c.BED):

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Validating BED ' + bed)

		try:

			bedfile=pybedtools.BedTool(bed)
			bedsrtd=bedfile.sort()

		except:

			now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
			print('[' + now + '][Error] BED ' + bed + ' does not exist, is not readable or is not a valid BED')
			sys.exit(1)

		d["h{0}".format(i+1)]=dict() #one sub-dict for each BED/haplotype. This way of specifying different haplotypes works perfectly

		for j,x in enumerate(bedsrtd):

			if x.chrom not in chrs:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] Line ' + str(j+1) + ': column 1 (chromosome name) contains an invalid chromosome (not included in the reference provided)')
				sys.exit(1)

			if x.start <= 0 or x.start > len(ref[x.chrom]):

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] Line ' + str(j+1) + ': column 2 (chromosome start) contains an invalid coordinate (lower than chromosome start or greater than chromosome end)')
				sys.exit(1)

			if x.end > len(ref[x.chrom]): #this can't be 0 I guess

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] Line ' + str(j+1) + ': column 3 (chromosome end) contains an invalid coordinate (greater than chromosome end)')				
				sys.exit(1)

			#check if 4th field is a valid/supported variant

			if x[3] not in possible_variants:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] Line ' + str(j+1) + ': column 4 (variant type) contains an unsupported variant type')
				sys.exit(1)

			#check if 6th field can be converted to integer

			try:

				int(x[5])

			except:

				now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
				print('[' + now + '][Error] Line ' + str(j+1) + ': column 6 (length of random sequence to insert at breakpoint) must contain an integer')
				sys.exit(1)

			#now validate infos (5th field) for each variant

			if x[3] == 'SNP':

				if x[4] not in list(valid_dna): #single base must be a valid base. This also checks for length greater than 1

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid DNA base')
					sys.exit(1)

				if int(x[5]) != 0: #no random sequence at breakpoint: not a SV

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Line ' + str(j+1) + ': column 6 (length of random sequence to insert at breakpoint) coherced to 0')

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],'')]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],''))

			elif x[3] == 'MNP':

				if not all(y in list(valid_dna) for y in str(x[4])): #check that every base is a valid base

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ':  column 5 (variant information) contains an invalid DNA sequence')
					sys.exit(1)

				if len(str(x[4])) != x.end-x.start +1: #check that the length of the sequence that has to be replaced mathces the length of a user-defined sequence

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ':  column 5 (variant information) contains a sequence shorter/longer than region')
					sys.exit(1)

				if int(x[5]) != 0: #no random sequence at breakpoint: not a SV

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Line ' + str(j+1) + ': column 6 (length of random sequence to insert at breakpoint) coherced to 0')

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],'')]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],''))

			elif x[3] == 'inversion':

				if x[4] != 'None': #tolerate not-None, as there is only one chance

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Line ' + str(j+1) + ':  column 5 (variant information) coherced to None')

				randomseq=''.join(random.choices(valid_dna, k=int(x[5]))) #if x[5] equal to 0, this is empty anyway

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],randomseq)]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],randomseq))

			elif x[3] == 'deletion':

				if x[4] not in {'None', '1bp'}: #these are the accepted possibilities by now

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid istruction. Must be either "None" or "1bp"')
					sys.exit[1]

				randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],randomseq)]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],randomseq))

			elif x[3] == 'insertion':

				if not (all(y in list(valid_dna) for y in x[4].upper())): #sequence to insert must be a valid DNA string

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid DNA sequence')
					sys.exit(1)

				randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4].upper(),randomseq)]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4].upper(),randomseq))

			elif x[3] == 'tandem duplication' or x[3] == 'inverted tandem duplication': #same checks for these variants

				try:

					int(x[4])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer')
					sys.exit(1)

				randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], int(x[4]),randomseq)]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], int(x[4]),randomseq))

			elif x[3] == 'perfect tandem repetition':

				column5=x[4].split(':') #info must be in this format

				if len(column5) != 2:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must be in string:integer format')
					sys.exit(1)

				if not (all(y in list(valid_dna) for y in column5[0].upper())): #sequence must be a valid DNA string

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid DNA motif')
					sys.exit(1)

				try:

					int(column5[1])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer specifying the number of repetitions')
					sys.exit(1)

				if int(x[5]) != 0: #no random sequence at breakpoint: not a SV

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Line ' + str(j+1) + ': column 6 (length of random sequence to insert at breakpoint) coherced to 0')

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],'')]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],''))

			elif x[3] == 'approximate tandem repetition':

				column5=x[4].split(':') #info must be in this format

				if len(column5) != 3:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must be in string:integer:integer format')
					sys.exit(1)

				if not (all(y in list(valid_dna) for y in column5[0].upper())): #sequence must be a valid DNA string

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid DNA motif')
					sys.exit(1)

				try:

					int(column5[1])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer specifying the number of repetitions')
					sys.exit(1)

				try:

					int(column5[2])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer specifying the number of errors in repetition')
					sys.exit(1)

				if int(x[5]) != 0: #no random sequence at breakpoint: not a SV

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Line ' + str(j+1) + ': column 6 (length of random sequence to insert at breakpoint) coherced to 0')

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],'')]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],''))

			elif x[3] == 'tandem repeat expansion' or x[3] == 'tandem repeat contraction': #same checks for these variants

				column5=x[4].split(':') #info must be in this format

				if len(column5) != 2:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must be in string:integer format')
					sys.exit(1)

				if not (all(y in list(valid_dna) for y in column5[0].upper())): #sequence must be a valid DNA string

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid DNA motif')
					sys.exit(1)

				try:

					int(column5[1])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer specifying the number of repetitions to add or subtract')
					sys.exit(1)

				if int(x[5]) != 0: #no random sequence at breakpoint: not a SV

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Warning] Line ' + str(j+1) + ': column 6 (length of random sequence to insert at breakpoint) coherced to 0')

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, x[3], x[4],'')]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, x[3], x[4],''))

			elif x[3] == 'reciprocal translocation':

				column5=x[4].split(':') #info must be in this format

				if len(column5) != 5:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must be in string:string:integer:string:string format')
					sys.exit(1)

				if not haplopattern.match(column5[0]):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain a valid haplotype string')
					sys.exit(1)

				if column5[1] not in chrs:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid chromosome string (chromosome is not included in the reference provided)')
					sys.exit(1)

				try:

					int(column5[2])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer specifying the breakpoint coordinate on the second chromosome')
					sys.exit(1)

				if column5[3] not in {'forward', 'reverse'} or column5[4] not in {'forward', 'reverse'}:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain valid orientations (either "forward" or "reverse")')
					sys.exit(1)

				#translocate second to first

				newtype='deletion-insertion'

				#get second sequence[2]
				firstbase=int(column5[2])
				lastbase=int(column5[2])+(x.end-x.start+1)

				if firstbase <= 0 or firstbase > len(ref[column5[1]]):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) breakpoint coordinate lies outside chromosome (lower than chromosome start or greater than chromosome end)')
					sys.exit(1)

				if lastbase > len(ref[column5[1]]): #last base can't be 0 as it is start (that can't be lower than 0) + something.

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) breakpoint coordinate lies outside chromosome (greater than chromosome end)')					
					sys.exit(1)
				
				if column5[4] == 'reverse':

					transeq=ref[column5[1]][firstbase:lastbase].reverse.complement.seq

				else:

					transeq=ref[column5[1]][firstbase:lastbase].seq

				randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

				if x.chrom not in d["h{0}".format(i+1)].keys(): #store

					d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, newtype, transeq,randomseq)]

				else:

					d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, newtype, transeq,randomseq))

				#translocate first to second

				if column5[0] not in d.keys():

					d[column5[0]] = dict() #initialize haplotype dict if not present

				#get first sequence
				
				if column5[3] == 'reverse':

					transeq=ref[x.chrom][x.start-1:x.end].reverse.complement.seq

				else:

					transeq=ref[x.chrom][x.start-1:x.end].seq

				randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

				if column5[1] not in d[column5[0]].keys(): #store

					d[column5[0]][column5[1]] = [(firstbase+1, lastbase, newtype, transeq,randomseq)]

				else:

					d[column5[0]][column5[1]].append((firstbase+1, lastbase, newtype, transeq,randomseq))

			elif x[3] == 'translocation cut-paste' or x[3] == 'translocation copy-paste' or x[3] == 'interspersed duplication': #same info for these 2

				column5=x[4].split(':') #info must be in this format

				if len(column5) != 4:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must be in string:string:integer:string format')
					sys.exit(1)

				if not haplopattern.match(column5[0]):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain a valid haplotype string')
					sys.exit(1)

				if column5[1] not in chrs:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) contains an invalid chromosome string (chromosome is not included in the reference provided)')
					sys.exit(1)

				try:

					int(column5[2])

				except:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain an integer specifying the breakpoint coordinate on the second chromosome')
					sys.exit(1)

				if int(column5[2]) <= 0 or int(column5[2]) > len(ref[column5[1]]):

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) breakpoint coordinate lies outside chromosome (lower than chromosome start or greater than chromosome end)')
					sys.exit(1)				

				if column5[3] not in {'forward', 'reverse'}:

					now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
					print('[' + now + '][Error] Line ' + str(j+1) + ': column 5 (variant information) must contain a valid orientation (either "forward" or "reverse")')
					sys.exit(1)

				if column5[0] not in d.keys():

					d[column5[0]] = dict() #initialize haplotype dict if not present

				if x[3] == 'translocation cut-paste':

					newtype1='deletion'
					newtype2='insertion'

					if column5[3] == 'reverse':

						transeq=ref[x.chrom][x.start-1:x.end].reverse.complement.seq

					else:

						transeq=ref[x.chrom][x.start-1:x.end].seq

					randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

					#delete first

					if x.chrom not in d["h{0}".format(i+1)].keys(): #store

						d["h{0}".format(i+1)][x.chrom] = [(x.start, x.end, newtype1, 'None',randomseq)]

					else:

						d["h{0}".format(i+1)][x.chrom].append((x.start, x.end, newtype1, 'None',randomseq))

					#insert in second

					randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

					if column5[1] not in d[column5[0]].keys(): #store

						d[column5[0]][column5[1]] = [(int(column5[2])-1,int(column5[2]),newtype2, transeq,randomseq)]

					else:

						d[column5[0]][column5[1]].append((int(column5[2])-1,int(column5[2]),newtype2, transeq,randomseq))

				else: #is copy-paste/intersperded dup

					newtype='insertion'

					if column5[3] == 'reverse':

						transeq=ref[x.chrom][x.start-1:x.end].reverse.complement.seq

					else:

						transeq=ref[x.chrom][x.start-1:x.end].seq

					randomseq=''.join(random.choices(valid_dna, k=int(x[5])))

					#only insert

					if column5[1] not in d[column5[0]].keys(): #store

						d[column5[0]][column5[1]] = [(int(column5[2])-1,int(column5[2]),newtype, transeq,randomseq)]

					else:

						d[column5[0]][column5[1]].append((int(column5[2])-1,int(column5[2]),newtype, transeq,randomseq))

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] BED validated and variants organized')
	print('[' + now + '][Message] Generating modified FASTA haplotypes')

	for dicts in d.keys():

		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Generating haplotype ' + dicts)
		hapout=os.path.abspath(c.OUT + '/' + dicts + '.fa')
		HapMaker(ref,chrs,d[dicts],hapout)
		now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
		print('[' + now + '][Message] Indexing FASTA haplotype')
		pyfaidx.Faidx(hapout)

	now=datetime.now().strftime('%d/%m/%Y %H:%M:%S')
	print('[' + now + '][Message] Done')
	sys.exit(0)

