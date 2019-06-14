#!/usr/bin/python env

#python 3 standard library

import os
import sys
import random
import logging
from operator import itemgetter
import string
import shutil
import re

#additional modules

import pybedtools #useful to sort inside python
import pyfaidx # fastest way to deal with .fasta in python


def run(parser,args):

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

		
	logging.basicConfig(filename=os.path.abspath(args.output + '/VISOR_HACk.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) #genome .file starts with '>'

	except:

		logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)



	immutable_ref=pyfaidx.Fasta(os.path.abspath(args.genome)) #load referene, that will be used to modify real .fasta
	classic_chrs = immutable_ref.keys() #allowed chromosomes
	possible_variants = ['deletion', 'insertion', 'inversion', 'tandem duplication', 'inverted tandem duplication', 'SNP', 'tandem repeat expansion', 'tandem repeat contraction', 'perfect tandem repetition', 'approximate tandem repetition', 'translocation cut-paste', 'translocation copy-paste', 'interspersed duplication', 'reciprocal translocation'] #allowed variants
	valid_dna = 'ACGT' #allowed nucleotides
	haplopattern=re.compile("^h[0-9]+$") #allowed haplotypes for inter-haplotypes SVs are h1,h2,h3 ...

	bedlist=[]

	for bed in args.bedfile[0]:

		if bed not in bedlist:

			bedlist.append(bed) #remove possible duplicates but mantain input order that should be lost in set


	d=dict() #initialize one empty dict for each haplotype


	logging.info('Organizing SVs')


	for i,bed in enumerate(bedlist):

		if not os.path.exists(os.path.abspath(bed)):

			logging.error('.bed file ' + os.path.abspath(bed) + ' does not exist')
			sys.exit(1)

		else:

			bedh = pybedtools.BedTool(os.path.abspath(bed))

			try:

				srtbedh = bedh.sort()

			except:

				logging.error('Incorrect .bed format for .bed file ' + os.path.abspath(bed))
				sys.exit(1)


			logging.info('Organizing SVs for ' + os.path.abspath(bed))

			
			d["h{0}".format(i+1)]=dict()


			for entries in srtbedh: 

				if str(entries[0]) not in classic_chrs: #exclude errors in the first field

					logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file ' + os.path.abspath(bed))
					sys.exit(1)

				try: #exclude errors in the second field

					int(entries[1])

				except:

					logging.error('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file  ' + os.path.abspath(bed) + '. Start must be an integer')
					sys.exit(1)


				try: #exclude errors in the third field

					int(entries[2])

				except:

					logging.error('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file  ' + os.path.abspath(bed) + '. End must be an integer')
					sys.exit(1)


				if (int(entries[2]) - int(entries[1]) == 0):

					logging.error('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed ' + os.path.abspath(bed) + '.')
					sys.exit(1)



				if str(entries[3]) not in possible_variants: #exclude errors in the third field

					logging.error(str(entries[3]) + ' is not a valid variant in .bed ' + os.path.abspath(bed))
					sys.exit(1)

		
				else: #everything fine for now


					if str(entries[3]) == 'SNP':

						if str(entries[4]) not in valid_dna: #information is just a valid nucleotide

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a valid DNA base included in A,C,T,G')
							sys.exit(1)


						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))

		
					elif str(entries[3]) == 'inversion': #information must be None

						if str(entries[4]) != 'None':

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be None')
							sys.exit(1)

						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


					elif str(entries[3]) == 'deletion': #information must be None

						if str(entries[4]) != 'None':

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be None')
							sys.exit(1)

						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))



					elif str(entries[3]) == 'insertion': #information must be a valid DNA sequence

						if not (all(i in valid_dna for i in entries[4].upper())): #validate sequence

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a valid DNA string with A,C,T,G characters')
							sys.exit(1)
								
						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]).upper())]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]).upper()))



					elif str(entries[3]) == 'tandem duplication': #information must be an integer

						try:

							int(entries[4])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be an integer')
							sys.exit(1)


						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4])))



					elif str(entries[3]) == 'inverted tandem duplication': #information must be an integer

						try:

							int(entries[4])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be an integer')
							sys.exit(1)


						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4])))


					elif str(entries[3]) == 'perfect tandem repetition': #perfect tandem repetition

						entr_4 = re.split('[:]',entries[4]) #Information must contain 2 fields

						if len(entr_4) != 2:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number')
							sys.exit(1)


						if not (all(i in valid_dna for i in entr_4[0].upper())): #check if motif is vaild DNA sequence

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
							sys.exit(1)

						try:

							int(entr_4[1])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
							sys.exit(1)

						
						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))



					elif str(entries[3]) == 'approximate tandem repetition': #approximative tandem repetition

						entr_4 = re.split('[:]',entries[4]) #Information must contain 2 fields

						if len(entr_4) != 3:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number:alterations')
							sys.exit(1)

						if not (all(i in valid_dna for i in entr_4[0].upper())):

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number:alterations. Motif must be a valid DNA motif')
							sys.exit(1)

						try:

							int(entr_4[1])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number:alterations. Number must be an integer')
							sys.exit(1)


						try:

							int(entr_4[2])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number:alterations. Alterations must be an integer')
							sys.exit(1)


						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


					
					elif str(entries[3]) == 'tandem repeat expansion' or str(entries[3]) == 'tandem repeat contraction': #Information must contain 2 fields

						entr_4 = re.split('[:]',entries[4])

						if len(entr_4) != 2:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number')
							sys.exit(1)

						if not (all(i in valid_dna for i in entr_4[0].upper())):

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
							sys.exit(1)

						try:

							int(entr_4[1])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
							sys.exit(1)


						if str(entries[0]) not in d["h{0}".format(i+1)].keys():

							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

						else:


							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))
					


					elif str(entries[3]) == 'reciprocal translocation':
								

						entr_4 = re.split('[:]',str(entries[4]))

						if len(entr_4) != 5:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2')
							sys.exit(1)


						if not haplopattern.match(str(entr_4[0])):

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Haplotype must be hN, with N being any integer.')
							sys.exit(1)

						if str(entr_4[1]) not in classic_chrs:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Chromosome must be a valid chromosome')
							sys.exit(1)

						try:

							int(entr_4[2])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Breakpoint must be an integer')
							sys.exit(1)



						if str(entr_4[3]) not in ['forward', 'reverse'] or str(entr_4[4]) not in ['forward', 'reverse']:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation1:orientation2. Orientation1 and orientation2 must be forward or reverse')
							sys.exit(1)



						if str(entries[0]) not in d["h{0}".format(i+1)].keys() and str(entr_4[4]) == 'forward':


							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), 'del-ins', immutable_ref[str(entr_4[1])][(int(entr_4[2])-1)+1:(int(entr_4[2])-1)+1+(int(entries[2])-int(entries[1]))].seq)]


						elif str(entries[0]) not in d["h{0}".format(i+1)].keys() and str(entr_4[4]) == 'reverse':

							
							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), 'del-invins', immutable_ref[str(entr_4[1])][(int(entr_4[2])-1)+1:(int(entr_4[2])-1)+1+(int(entries[2])-int(entries[1]))].seq)]


						elif str(entries[0]) in d["h{0}".format(i+1)].keys() and str(entr_4[4]) == 'forward':

							
							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), 'del-ins', immutable_ref[str(entr_4[1])][(int(entr_4[2])-1)+1:(int(entr_4[2])-1)+1+(int(entries[2])-int(entries[1]))].seq))


						elif str(entries[0]) in d["h{0}".format(i+1)].keys() and str(entr_4[4]) == 'reverse':

							
							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), 'del-invins', immutable_ref[str(entr_4[1])][(int(entr_4[2])-1)+1:(int(entr_4[2])-1)+1+(int(entries[2])-int(entries[1]))].seq))


						
						if str(entr_4[0]) not in d.keys():

							d[entr_4[0]]=dict()


						if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) =='forward':


							d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2])+1, int(entr_4[2])+1+(int(entries[2])-int(entries[1])), 'del-ins', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]



						elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':


							d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2])+1, int(entr_4[2])+1+(int(entries[2])-int(entries[1])), 'del-invins', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


						elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':


							d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2])+1, int(entr_4[2])+1+(int(entries[2])-int(entries[1])), 'del-ins', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


						elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

							d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2])+1, int(entr_4[2])+1+(int(entries[2])-int(entries[1])), 'del-invins', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


					elif str(entries[3]) == 'translocation cut-paste':
								

						entr_4 = re.split('[:]',str(entries[4]))

						if len(entr_4) != 4:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation.')
							sys.exit(1)


						if not haplopattern.match(str(entr_4[0])):

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Haplotype must be hN, with N being any integer.')
							sys.exit(1)

						if str(entr_4[1]) not in classic_chrs:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Chromosome must be a valid chromosome.')
							sys.exit(1)

						try:

							int(entr_4[2])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Breakpoint must be an integer')
							sys.exit(1)


						if str(entr_4[3]) not in ['forward', 'reverse']:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Orientation must be forward or reverse')
							sys.exit(1)




						if str(entries[0]) not in d["h{0}".format(i+1)].keys():


							d["h{0}".format(i+1)][str(entries[0])] = [(int(entries[1]), int(entries[2]), 'deletion', 'None')]

						else:

							d["h{0}".format(i+1)][str(entries[0])].append((int(entries[1]), int(entries[2]), 'deletion', 'None'))

					

						if str(entr_4[0]) not in d.keys():

							d[entr_4[0]]=dict()


						if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) =='forward':


							d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2])-1, int(entr_4[2]), 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


						elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':


							d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2])-1, int(entr_4[2]), 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


						elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':


							d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2])-1, int(entr_4[2]), 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


						elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

							d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2])-1, int(entr_4[2]), 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))



			
					else: #is a translocation copy paste or interspersed duplication, they are the same


						entr_4 = re.split('[:]',str(entries[4]))

						if len(entr_4) != 4:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation.')
							sys.exit(1)


						if not haplopattern.match(str(entr_4[0])):

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Haplotype must be hN, with N being any integer.')
							sys.exit(1)

						if str(entr_4[1]) not in classic_chrs:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Chromosome must be a valid chromosome.')
							sys.exit(1)

						try:

							int(entr_4[2])

						except:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Breakpoint must be an integer')
							sys.exit(1)


						if str(entr_4[3]) not in ['forward', 'reverse']:

							logging.error('Incorrect info ' + str(entries[4]) + ' in .bed ' + os.path.abspath(bed) + ' for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation. Orientation must be forward or reverse')
							sys.exit(1)

					

						if str(entr_4[0]) not in d.keys():

							d[entr_4[0]]=dict()


						if str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) =='forward':


							d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2])-1, int(entr_4[2]), 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


						elif str(entr_4[1]) not in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':


							d[str(entr_4[0])][str(entr_4[1])] = [(int(entr_4[2])-1, int(entr_4[2]), 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


						elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'forward':


							d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2])-1, int(entr_4[2]), 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


						elif str(entr_4[1]) in d[str(entr_4[0])].keys() and str(entr_4[3]) == 'reverse':

							d[str(entr_4[0])][str(entr_4[1])].append((int(entr_4[2])-1, int(entr_4[2]), 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


	logging.info('SVs organized')						
	logging.info('Generating .fasta haplotypes with SVs')

	for dicts in d.keys():
		
		logging.info('Generating SVs for ' + str(dicts))
		ParseDict(classic_chrs, immutable_ref, d[dicts], os.path.abspath(args.output + '/' + str(dicts) + '.fa'))

	logging.info('Haplotypes with SVs generated')
	logging.info('Done')




def write_unmodified_chromosome(chromosome, seq, output_fasta):

	with open (os.path.abspath(output_fasta), 'a') as faout:

		faout.write('>' + chromosome + '\n' + seq + '\n')


def write_start_sequence(chromosome, seq, output_fasta):

	with open (os.path.abspath(output_fasta), 'a') as faout:

		faout.write('>' + chromosome + '\n' + seq )



def write_sequence_between(seq, output_fasta):

	with open (os.path.abspath(output_fasta), 'a') as faout:

		faout.write(seq)



def write_end_sequence(seq, output_fasta):

	with open (os.path.abspath(output_fasta), 'a') as faout:

		faout.write(seq + '\n')



def Change_Random_Char(word):

	length = len(word)
	word = list(word)
	k = random.sample(range(0,length),1)
	k.sort()
	nuc_list=['A','T','C','G']
   
	for index in k:

		add_rem=word[index]
		nuc_list.remove(add_rem)
		word[index] = ''.join(random.sample(nuc_list,k=1))
		nuc_list.append(add_rem)
   
	return('' . join(word))


def Delete_Random_Char(word):

	index = random.randint(0, len(word)-1)
	word = word[:index] + word[index+1:]

	return word


def Insert_Random_Char(word):

	nucs=['A','T','C','G']

	index = random.randint(0, len(word)-1)
	word = word[:index] + random.choice(nucs) + word[index:]

	return word



def Reverse(sequence, start, end):

	new_seq = sequence[start-1:end][::-1]

	return new_seq


def PTR(infofield, sequence, start, end): #new ptr

	info=re.split('[:]', infofield)
	motif,length = str(info[0]), int(info[1])
	new_seq= sequence[start-1:end] + motif*length

	return new_seq


def ATR(infofield, sequence, start, end): #new atr

	info=re.split('[:]', infofield)
	motif,length,altnum = str(info[0]), int(info[1]), int(info[2])
	new_seq= motif*length

	alterations=['insertion', 'deletion', 'substitution']

	counter=0

	while counter < altnum:

		alt_type=random.choice(alterations)

		if alt_type == 'substitution':

			new_seq=Change_Random_Char(new_seq)

		elif alt_type == 'deletion':

			new_seq=Delete_Random_Char(new_seq)

		else: #insertion

			new_seq= Insert_Random_Char(new_seq)

		counter +=1

	return sequence[start-1:end] + new_seq


def EXPTR(infofield, sequence, start, end): #expand tr

	info=re.split('[:]', infofield)
	motif,num=str(info[0]), int(info[1])
	trseq=sequence[start-1:end]

	firstbase=trseq[0]
	exprep=trseq[1:] + motif*num

	new_seq=firstbase + exprep

	return new_seq


def CTRTR(infofield, sequence, start, end): #contract tr

	info=re.split('[:]', infofield)
	motif,num=str(info[0]), int(info[1])

	trseq=sequence[start-1:end]

	firstbase=trseq[0]
	rep=trseq[1:]
	newind=len(motif)*num
	delrep=rep[newind:]
	new_seq=firstbase + delrep

	return new_seq




def ParseDict(chromosomes, fasta, dictionary, output_fasta):

	trans = str.maketrans('ATGC', 'TACG')
	skipped=0

	for chrs in chromosomes:

		chrom=fasta[chrs]
		seq=chrom[:len(chrom)].seq


		if chrs not in dictionary.keys(): #chromosome not there, write unchanged
						
			write_unmodified_chromosome(chrs, seq, output_fasta)

		else:

			alterations_list=dictionary[chrs]

			if not len(alterations_list) == 1: #else is already sorted, as it has length 1

				alterations_list=sorted(alterations_list, key=itemgetter(0,1))
				
				new_alterations_list=[]

				l=0
				
				while l < len(alterations_list):

					start,end,typ,info=alterations_list[l]

					if new_alterations_list==[]:

						new_alterations_list.append((start,end,typ,info))

					else:

						if (new_alterations_list[-1][0]<= start <= new_alterations_list[-1][1]) or (new_alterations_list[-1][0] <= end <= new_alterations_list[-1][0]):

							skipped+=1

						else:

							new_alterations_list.append((start,end,typ,info))

					l+=1

			else:

				new_alterations_list=alterations_list

			i=0

			while i < len(new_alterations_list):

				start,end,typ,info=new_alterations_list[i]

				if start==0:

					start+=1

				if i == 0: #first entry for the cromosome, write until the first variant start

					seq_until_start=seq[:start-1]

					write_start_sequence(chrs, seq_until_start, output_fasta)

				if typ == 'inversion': #inverte sequence

					alt_seq=Reverse(seq, start, end).translate(trans)

					write_sequence_between(alt_seq, output_fasta)


				elif typ == 'deletion': #write nothing; deletions and insertions are also valid for translocation, are they are translated before intro insertions and deletions

					alt_seq=''

					write_sequence_between(alt_seq, output_fasta)


				elif typ == 'insertion': #write specified insertion; deletions and insertions are also valid for translocation, are they are translated before intro insertions and deletions

					alt_seq=info

					write_sequence_between(seq[start-1:end]+alt_seq, output_fasta)

				elif typ == 'invinsertion':

					alt_seq=info[::-1].translate(trans)

					write_sequence_between(seq[start-1:end]+alt_seq, output_fasta)


				elif typ == 'del-ins':


					write_sequence_between(info, output_fasta)


				elif typ == 'del-invins':

					alt_seq=info[::-1].translate(trans)


					write_sequence_between(alt_seq, output_fasta)


				elif typ == 'tandem duplication':

					write_sequence_between(seq[start-1:end]*info, output_fasta)


				elif typ == 'inverted tandem duplication':

					write_sequence_between(seq[start-1:end] + (Reverse(seq, start, end).translate(trans) * (info-1)), output_fasta) #first part is not inverted, duplicated part it is

				elif typ == 'SNP':

					until_end=seq[start-1:end-1]

					write_sequence_between(until_end+info, output_fasta)


				elif typ == 'perfect tandem repetition': #perfect tandem repetition

					alt_seq= PTR(info, seq, start, end)

					write_sequence_between(alt_seq, output_fasta)


				elif typ == 'approximate tandem repetition': # approximate tandem repetition


					alt_seq=ATR(info,seq,start,end)

					write_sequence_between(alt_seq, output_fasta)


				elif typ == 'tandem repeat expansion': #expand a tandem repetition that is already present. Start-end are supposed to be as the one in repetitions .bed from ucsc


					alt_seq=EXPTR(info,seq,start,end)

					write_sequence_between(alt_seq, output_fasta)

				elif typ == 'tandem repeat contraction': #contract a tandem repetition that is already present. Start-end are supposed to be as the one in repetitions .bed from ucsc

					alt_seq=CTRTR(info,seq,start,end)
				
					write_sequence_between(alt_seq, output_fasta)


				if i == len(new_alterations_list) -1:

					write_end_sequence(seq[end:], output_fasta) #end not included, as it was included in the variant


				elif i < len(new_alterations_list) -1:


					nextstart=new_alterations_list[i+1][0]
					thisend=end

					write_sequence_between(seq[thisend:nextstart-1], output_fasta)

				i+=1

	if skipped > 0 :

		logging.warning('Skipped ' + str(skipped) + ' SVs for the current haplotype as they overlapped others')
