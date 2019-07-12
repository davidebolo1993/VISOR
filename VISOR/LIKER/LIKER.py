#!/usr/bin/python env

#python 3 standard library

import os
import sys
import glob
import re
import math
from shutil import which
import multiprocessing
import subprocess, shlex
import gzip

#additional modules

import pybedtools
import pyfaidx
import numpy as np

global barcodes

def run(parser,args):

	global barcodes

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

	
	external_tools=['wgsim']

	for tools in external_tools:

		if which(tools) is None:

			print(tools + ' was not found as an executable command. Install ' + tools + ' and re-run VISOR LIKER')
			sys.exit(1)

	bed = pybedtools.BedTool(os.path.abspath(args.bedfile))

	try:

		srtbed = bed.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		print('Incorrect .bed format for -bed/--bedfile')
		sys.exit(1)


	fastas = glob.glob(os.path.abspath(args.sample + '/*.fa'))

	if fastas == []:

		print('Given folder ' + os.path.abspath(args.sample) + ' does not contain any valid .fasta inputs')
		sys.exit(1)


	barcodepath=os.path.abspath(os.path.dirname(__file__) + '/4M-with-alts-february-2016.txt.gz')

	with gzip.open(barcodepath, 'rt') as fin:

		barcodes = fin.read().splitlines()


	Par=Container()

	Par.NMOL=args.molecules_number
	Par.LMOL=args.molecules_length
	Par.SRL=args.length
	Par.SRE=args.error
	Par.INS=args.insertsize
	Par.STDEV=args.standardev
	Par.INDELP=args.indels
	Par.INDELEXT=args.probability
	Par.NAME=args.identifier
	Par.PROCS=args.threads
	Par.COV=args.coverage/len(fastas)


	for folder,fasta in enumerate(fastas):

		print('Simulating from haplotype ' + os.path.abspath(fasta))

		Fa=pyfaidx.Fasta(fasta)
		allchrs=Fa.keys()

		os.makedirs(os.path.abspath(args.output + '/h' + str(folder+1))) #create directory for this haplotype

		Par.HAPLONUM=folder+1

		counter=0

		for entries in srtbed: #validate each entry

			counter+=1

			if str(entries[0]) not in allchrs:

				print(str(entries[0]) + ' is not a valid chromosome in .bed file. Skipped')
				continue					

			try:

				int(entries[1])

			except:

				print('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file. Start must be an integer')
				sys.exit(1)


			try:

				int(entries[2])

			except:

				print('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file. End must be an integer')
				sys.exit(1)


			if (int(entries[2]) - int(entries[1]) == 0):

				print('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed file')
				sys.exit(1)



			with open(os.path.abspath(args.output + '/h' + str(folder+1) + '/' + str(counter) +'.region.tmp.fa'), 'w') as regionout:

				subprocess.call(['samtools', 'faidx', os.path.abspath(fasta), str(entries[0]) + ':' + str(entries[1]) +  '-' +str(entries[2])], stdout=regionout, stderr=open(os.devnull, 'wb'))


			Par.COUNT=counter

			FAREGION=pyfaidx.Fasta(os.path.abspath(args.output + '/h' + str(folder+1) + '/' + str(counter) +'.region.tmp.fa'))
			CHROMREGION=FAREGION[str(entries[0]) + ':' + str(entries[1]) +  '-' +str(entries[2])]
			SEQREGION=CHROMREGION[:len(CHROMREGION)].seq
			REGIONSTART=int(entries[1])

			print(Par)

			LinkedSim(Par, str(entries[0]), SEQREGION, REGIONSTART, os.path.abspath(args.output + '/h' + str(folder+1)))


	#organize output in output folder

	print('Storing in proper format for longranger')

	subdirs=sorted([os.path.join(os.path.abspath(args.output), o) for o in os.listdir(os.path.abspath(args.output)) if os.path.isdir(os.path.join(os.path.abspath(args.output),o))], key=natural_keys)

	counter=0

	for dirs in subdirs:

		counter +=1

		R1=[]
		R2=[]

		R1.extend(glob.glob(os.path.abspath(dirs) + '/*_R1_*'))
		R2.extend(glob.glob(os.path.abspath(dirs) + '/*_R2_*'))

		command1 = 'cat ' + ' '.join(x for x in R1)
		command2 = 'cat ' + ' '.join(x for x in R2)

		with open(os.path.abspath(args.output + '/' + Par.NAME + '_S1_L' + str(counter).zfill(3) + '_R1_001.fastq'), 'w') as fout:

			subprocess.call(shlex.split(command1), stdout=fout)

		subprocess.call(['gzip', os.path.abspath(args.output + '/' + Par.NAME + '_S1_L' + str(counter).zfill(3) + '_R1_001.fastq')])

		with open(os.path.abspath(args.output + '/' + Par.NAME + '_S1_L' + str(counter).zfill(3) + '_R2_001.fastq'), 'w') as fout:

			subprocess.call(shlex.split(command2), stdout=fout)

		subprocess.call(['gzip', os.path.abspath(args.output + '/' + Par.NAME + '_S1_L' + str(counter).zfill(3) + '_R2_001.fastq')])

		for x in R1:

			os.remove(x)


		for y in R2:

			os.remove(y)

		os.rmdir(dirs)

	print('Done')


#CLASSES


class Molecule(object):
		
	def __init__(self,length,start,end,index):
				
		self.seqidx=index
		self.length=length
		self.index_droplet=0
		self.barcode=None
		self.start=start
		self.end=end

	def __str__(self):
		
		return str(self.__class__) + ": " + str(self.__dict__)


class Container(object):
		
	def __init__(self):
		
		self.HAPLONUM=0
		self.COV=0
		self.NMOL=0
		self.LMOL=0
		self.SRL=0
		self.SRE=0
		self.INS=0
		self.STDEV=0
		self.INDELP=0
		self.INDELEXT=0
		self.NAME=None
		self.PROCS=0
		self.COUNT=0

	def __str__(self):
		
		return str(self.__class__) + ": " + str(self.__dict__)



#FUNCTIONS


def Chunks(l,n):

	return [l[i:i+n] for i in range(0, len(l), n)]


def atoi(text):

	return int(text) if text.isdigit() else text

def natural_keys(text):
	
	return [ atoi(c) for c in re.split(r'(\d+)', text)]



def randomlong(Par,refseq):

	index=0
			
	lensingle=len(refseq)
	N_frag=round(lensingle*Par.COV/(Par.LMOL)) #how many fragment we generate for a certain region 
				
	for i in range(N_frag):
						
		start=int(np.random.uniform(low=0,high=lensingle))
		length=int(np.random.exponential(scale=Par.LMOL)) #molecules length is randomly distributed according to exponential distribution

		if length==0:
					
			continue

		end=start+length-1

		if end>lensingle:
							 
			Molseq=refseq[start:lensingle]
			lengthnew=lensingle-start
			NewMol=Molecule(lengthnew,start,lensingle,index)
			MolSet.append(NewMol)
						

		else:
							 
			Molseq=refseq[start:end]
			NewMol=Molecule(length-1,start,end,index)
			MolSet.append(NewMol)
		
		index+=1	



def deternumdroplet(NMOL_SET,NMOL_MEAN):

	large_droplet=4000000

	frag_drop = np.random.poisson(NMOL_MEAN,large_droplet)
	totalfrag=0
	
	for i in range(large_droplet):
				
		totalfrag=totalfrag+frag_drop[i]
		
		if totalfrag<=len(NMOL_SET):
		
			assign_drop.append(frag_drop[i])
		
		else:
			
			last=len(NMOL_SET)-(totalfrag-frag_drop[i])
			assign_drop.append(last)
			break



def selectbarcode(assign_drop,MolSet,droplet_container):

	permutnum=np.random.permutation(len(MolSet))
	N_droplet=len(assign_drop)

	start=0
	
	for i in range(N_droplet):
				
		num_molecule_per_partition=assign_drop[i]
		index_molecule=permutnum[start:start+num_molecule_per_partition]
		#print(index_molecule)
		totalseqlen=0
		temp=[]
		start=start+num_molecule_per_partition
	
		for j in range(num_molecule_per_partition):
						
			index=index_molecule[j]
			temp.append(index)
			MolSet[index].index_droplet=i
			MolSet[index].barcode=barcodes[i]
			totalseqlen=totalseqlen+MolSet[index].length
		
		assigned_barcodes.add(barcodes[i])
		droplet_container.append(temp)




def Runner(processor,molecule,refseq,reftitle, refstart, Par, output):

	#remained=0

	for mol in molecule:

		moleculenumber=str(mol.seqidx+1)
		moleculedroplet=str(mol.index_droplet+1)
		barcodestring=str(mol.barcode)
		chromstart=str(refstart+mol.start-1)
		chromend=str(refstart+mol.end)

		header='>MOL:' + moleculenumber + '_GEM:' + moleculedroplet + '_BAR:' + barcodestring + '_CHROM:' + reftitle + '_START:' + chromstart + '_END:' + chromend
		seq=refseq[mol.start:mol.end+1]

		with open (os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.fa'), 'w') as fout:

			fout.write(header + '\n' + seq + '\n')

		NUM_READS=int(((mol.length-seq.count('N')))/(Par.SRL*2))


		if NUM_READS == 0: #got a region with only Ns

			#remained+=mol.length
			continue

		#remained=0

		subprocess.call(['wgsim', '-e', str(Par.SRE), '-d', str(Par.INS), '-s', str(Par.STDEV), '-N', str(NUM_READS), '-1', str(Par.SRL-(16+6)), '-2', str(Par.SRL), '-R', str(Par.INDELP), '-X', str(Par.INDELEXT), os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.fa'), os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.R1.fq.tmp'), os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.R2.fq')], stderr=open(os.devnull, 'wb'), stdout=open(os.devnull, 'wb'))

		#reopen FASTQ and modify adding barcode and random 6MER, with relative quality (default)

		RANDOM6MER=''.join(np.random.choice(['A','T','G','C','N'],6, replace=True))

		with open(os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.R1.fq.tmp'), 'r') as fin, open(os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.R1.fq'), 'w') as fout:

			for line in fin:

				if line.startswith(str(5)):

					fout.write(str(5)*(16+6) + line)

				elif line.startswith('@') or line.startswith('+'):

					fout.write(line)

				else:

					fout.write(barcodestring + RANDOM6MER + line)

		os.remove(os.path.abspath(output + '/' + processor + '.' + moleculenumber + '.R1.fq.tmp'))




def LinkedSim(Par,reftitle,refseq,refstart, output):

	global MolSet    
	global assign_drop
	global assigned_barcodes
	global droplet_container
	global barcodes
	
	MolSet=[]
	assign_drop=[]
	assigned_barcodes=set()
	droplet_container=[]

	print(str(len(barcodes)) + ' barcodes available')
	
	randomlong(Par,refseq)
		
	print('Generated ' + str(len(MolSet)) + ' molecules')

	deternumdroplet(MolSet,Par.NMOL)

	print('Assigned molecules to GEMs')
		
	selectbarcode(assign_drop,MolSet,droplet_container)
		
	print('Assigned a barcode to each molecule')

	print(str(len(assigned_barcodes)) + ' barcodes assigned')

	barcodes=[e for e in barcodes if e not in assigned_barcodes]

	chunk_size=len(MolSet)/Par.PROCS
	slices=Chunks(MolSet,math.ceil(chunk_size))

	#parallelize FASTQ creation

	processes=[]

	for i,molecule in enumerate(slices):

		processor='p'+str(i+1)
		print('processor ' + processor + ' is ready to fly')
		
		p=multiprocessing.Process(target=Runner, args=(processor,molecule,refseq,reftitle, refstart, Par, output))
		
		p.start()
		processes.append(p)
		
	for p in processes:
		
		p.join()

	R1=glob.glob(os.path.abspath(output +'/p*R1.fq'))
	R2=glob.glob(os.path.abspath(output +'/p*R2.fq'))
	R3=glob.glob(os.path.abspath(output +'/*.fa*'))

	command1 = 'cat ' + ' '.join(x for x in R1)
	command2 = 'cat ' + ' '.join(x for x in R2)

	with open(os.path.abspath(output + '/' + str(Par.COUNT) + '.' + Par.NAME + '_S1_L' + str(Par.HAPLONUM).zfill(3) + '_R1_001.fastq'), 'w') as fout:

		subprocess.call(shlex.split(command1), stdout=fout)

	with open(os.path.abspath(output + '/' + str(Par.COUNT) + '.' + Par.NAME + '_S1_L' + str(Par.HAPLONUM).zfill(3) + '_R2_001.fastq'), 'w') as fout:

		subprocess.call(shlex.split(command2), stdout=fout)

	for f1,f2 in zip(R1,R2):

		os.remove(f1)
		os.remove(f2)

	for f3 in R3:

		os.remove(f3)
