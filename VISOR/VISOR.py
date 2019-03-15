#!/usr/bin/python env

#Python 3 standard library

import os
import sys
import random
import logging
import re
from operator import itemgetter
import argparse
from argparse import HelpFormatter
import timeit
import string


#additional libraries

import pybedtools #useful to sort inside python
import pyfaidx # fastest way to deal with .fasta in python


def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 


	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference file', metavar='.fa', required=True)
	required.add_argument('-bedh1', '--bedfile_haplotype_1', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 1', metavar='.bed', required=True)
	required.add_argument('-O', '--output', help='name of the directory where the 2 .fa haplotype will be saved', metavar='folder', required=True)

	optional = parser.add_argument_group('Additional input')
	optional.add_argument('-bedh2', '--bedfile_haplotype_2', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 2', metavar='', default=None)

	args = parser.parse_args()


	if not os.path.exists(os.path.abspath(args.output)):

		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			#print('It was not possible to create the results folder. Specify a path for which you have write permissions')
			logging.error('It was not possible to create the results folder. Specify a path for which you have write permissions')
			sys.exit(1)


	logging.basicConfig(filename=os.path.abspath(args.output + '/VISOR.log'), filemode='w', level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

	start=timeit.default_timer()



	classic_chrs = ['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y', 'M']] #allowed chromosomes
	possible_variants = ['deletion', 'insertion', 'inversion', 'duplication', 'tr expansion', 'tr contraction', 'ptr', 'atr', 'translocation cut-paste', 'translocation copy-paste'] #allowed variants
	valid_dna = 'NACGT' #allowed nucleotides
	

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>')) #genome .file starts with '>'

	except:

		#print('Reference file does not exist, is not readable or is not a valid .fasta file')
		logging.error('Reference file does not exist, is not readable or is not a valid .fasta file')
		sys.exit(1)



	immutable_ref=pyfaidx.Fasta(os.path.abspath(args.genome)) #load referene, that will be used to modify real .fasta

	bedh1 = pybedtools.BedTool(os.path.abspath(args.bedfile_haplotype_1)) #this one is required
	varh1=dict()

	
	try:
	
		srtbedh1 = bedh1.sort() #sorting here is not necessary, but allows to look for general formatting errors

	except:

		#print('Incorrect .bed format for haplotype 1')
		logging.error('Incorrect .bed format for haplotype 1')
		sys.exit(1)



	logging.info('Organizing variants')

	#initialize dictionaries

	hap1dict=dict() # will be filled
	hap2dict=dict() # can be empty at the end


	for entries in srtbedh1: #read, validate and organize

		if str(entries[0]) not in classic_chrs: #exclude errors in the first field

			#print(str(entries[0]) + ' is not a valid chromosome in .bed file for haplotype 1. Allowed chromosomes are chr1-22, chrX, chrY and chrM')
			logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file for haplotype 1. Allowed chromosomes are chr1-22, chrX, chrY and chrM')
			sys.exit(1)

		try:

			int(entries[1])

		except:

			#print('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file for haplotype 1. Start must be an integer')
			logging.error('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file for haplotype 1. Start must be an integer')
			sys.exit(1)


		try:

			int(entries[2])

		except:

			#print('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file for haplotype 1. End must be an integer')
			logging.error('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file for haplotype 1. End must be an integer')
			sys.exit(1)


		if (int(entries[2]) - int(entries[1]) == 0):

			#print('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed for haplotype 1')
			logging.error('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed for haplotype 1')
			sys.exit(1)


		if int(entries[1]) == 0:

			#print('Start ' + str(entries[1]) + ' cannot be 0')
			logging.error('Start ' + str(entries[1]) + ' cannot be 0')
			sys.exit(1)


		if str(entries[3]) not in possible_variants: #stop on variant different than expected

			#print(str(entries[3]) + ' is not a valid variant in .bed for haplotype 1.')
			logging.error(str(entries[3]) + ' is not a valid variant in .bed for haplotype 1.')
			sys.exit(1)

		
		else: #everything fine
		
			if str(entries[3]) == 'inversion':

				if str(entries[4]) != 'None':

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be None')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be None')
					sys.exit(1)

				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

				else:

					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


			elif str(entries[3]) == 'deletion':

				if str(entries[4]) != 'None':

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be None')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be None')
					sys.exit(1)

				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

				else:

					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


			elif str(entries[3]) == 'insertion':

				if not (all(i in valid_dna for i in entries[4].upper())): #validate sequence

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a DNA sequence with only A,C,T,G,N chars')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a DNA sequence with A,C,T,G,N chars')
					sys.exit(1)
						
				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]).upper())]

				else:

					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]).upper()))


			elif str(entries[3]) == 'duplication':

				try:

					int(entries[4])

				except:

					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with number of duplication. Number must be an integer')
					sys.exit(1)

				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4]))]

				else:

					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4])))


			elif str(entries[3]) == 'ptr': #perfect tandem repetition

				entr_4 = re.split('[:]',entries[4])

				if len(entr_4) != 2:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
					sys.exit(1)


				elif not (all(i in valid_dna for i in entr_4[0].upper())): #check if motif is vaild DNA sequence

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
					sys.exit(1)

				try:

					int(entr_4[1])

				except:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
					sys.exit(1)

				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

				else:

					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


			elif str(entries[3]) == 'atr': #approximative tandem repetition

				entr_4 = re.split('[:]',entries[4])

				if len(entr_4) != 3:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum')
					sys.exit(1)

				elif not (all(i in valid_dna for i in entr_4[0].upper())):

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Motif must be a valid DNA motif')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Motif must be a valid DNA motif')
					sys.exit(1)

				try:

					int(entr_4[1])

				except:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Number must be an integer')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Number must be an integer')
					sys.exit(1)


				try:

					int(entr_4[2])

				except:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Altnum must be an integer')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Altnum must be an integer')
					sys.exit(1)


				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

				else:

					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))



			elif str(entries[3]) == 'tr expansion' or str(entries[3]) == 'tr contraction':

				entr_4 = re.split('[:]',entries[4])

				if len(entr_4) != 2:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
					sys.exit(1)

				elif not (all(i in valid_dna for i in entr_4[0].upper())):

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number.  Motif must be a valid DNA motif')
					sys.exit(1)

				try:

					int(entr_4[1])

				except:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
					sys.exit(1)


				if str(entries[0]) not in hap1dict:

					hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

				else:


					hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))
					

			elif str(entries[3]) == 'translocation cut-paste':
						

				entr_4 = re.split('[:]',str(entries[4]))

				if len(entr_4) != 4:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
					sys.exit(1)


				if str(entr_4[0]) not in ['h1','h2']:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
					sys.exit(1)

				if str(entr_4[1]) not in classic_chrs:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
					sys.exit(1)

				try:

					int(entr_4[2])

				except:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
					sys.exit(1)


				if str(entr_4[3]) not in ['forward', 'reverse']:

					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Orientation can be forward or reverse')
					sys.exit(1)


				if entr_4[0] == 'h1':

					if str(entries[0]) not in hap1dict:

						hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), 'deletion', 'None')]

					else:

						hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), 'deletion', 'None'))

					if str(entr_4[1]) not in hap1dict:

						if str(entr_4[3]) == 'forward':

							hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else: #orientation is reverse

							hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


					else:

						if str(entr_4[3]) == 'forward':


							hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

						else: #orientation is reverse

							hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


				elif entr_4[0] == 'h2':

					if str(entries[0]) not in hap1dict:

						hap1dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), 'deletion', 'None')]

					else:

						hap1dict[str(entries[0])].append((int(entries[1]), int(entries[2]), 'deletion', 'None'))

					if str(entr_4[1]) not in hap2dict:

						if str(entr_4[3]) == 'forward':

							hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else: #orientation is reverse

							hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

					else:

						if str(entr_4[3]) == 'forward':


							hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

						else: #orientation is reverse

							hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

			
			else: #is a translocation copy paste


				entr_4 = re.split('[:]',str(entries[4]))

				if len(entr_4) != 4:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
					sys.exit(1)


				if str(entr_4[0]) not in ['h1','h2']:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
					sys.exit(1)

				if str(entr_4[1]) not in classic_chrs:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
					sys.exit(1)

				try:

					int(entr_4[2])

				except:

					#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
					sys.exit(1)

				if str(entr_4[3]) not in ['forward', 'reverse']:

					logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 1 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Orientation can be forward or reverse')
					sys.exit(1)

				if entr_4[0] == 'h1':

					if str(entr_4[1]) not in hap1dict:

						if str(entr_4[3]) == 'forward':

							hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else: #orientation is reverse

							hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

					else:


						if str(entr_4[3]) == 'forward':


							hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

						else: #orientation is reverse

							hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


				elif entr_4[0] == 'h2':

					if str(entr_4[1]) not in hap2dict:

						if str(entr_4[3]) == 'forward':

							hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else: #orientation is reverse

							hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

					else:

						if str(entr_4[3]) == 'forward':


							hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

						else: #orientation is reverse

							hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))



	if not args.bedfile_haplotype_2 is None: #also second .bed provided

		bedh2 = pybedtools.BedTool(os.path.abspath(args.bedfile_haplotype_2))
		varh2=dict()
	
		try:

			srtbedh2 = bedh2.sort()

		except:

			#print('Incorrect .bed format for haplotype 2')
			logging.error('Incorrect .bed format for haplotype 2')
			sys.exit(1)

		for entries in srtbedh2: #read, validate and organize

			if str(entries[0]) not in classic_chrs: #exclude errors in the first field

				#print(str(entries[0]) + ' is not a valid chromosome in .bed file for haplotype 2. Allowed chromosomes are chr1-22, chrX, chrY and chrM')
				logging.error(str(entries[0]) + ' is not a valid chromosome in .bed file for haplotype 2. Allowed chromosomes are chr1-22, chrX, chrY and chrM')
				sys.exit(1)

			try:

				int(entries[1])

			except:

				#print('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file for haplotype 2. Start must be an integer')
				logging.error('Cannot convert ' + str(entries[1]) + ' to integer number in .bed file for haplotype 2. Start must be an integer')
				sys.exit(1)


			try:

				int(entries[2])

			except:

				#print('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file for haplotype 2. End must be an integer')
				logging.error('Cannot convert ' + str(entries[2]) + ' to integer number in .bed file for haplotype 2. End must be an integer')
				sys.exit(1)


			if (int(entries[2]) - int(entries[1]) == 0):

				#print('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed for haplotype 2')
				logging.error('Start ' + str(entries[1]) + ' and end ' + str(entries[2]) + ' cannot have the same value in .bed for haplotype 2')
				sys.exit(1)


			if int(entries[1]) == 0:

				#print('Start ' + str(entries[1]) + ' cannot be 0')
				logging.error('Start ' + str(entries[1]) + ' cannot be 0')
				sys.exit(1)


			if str(entries[3]) not in possible_variants: #stop on variant different than expected

				#print(str(entries[3]) + ' is not a valid variant in .bed for haplotype 2.')
				logging.error(str(entries[3]) + ' is not a valid variant in .bed for haplotype 2.')
				sys.exit(1)

			
			else: #everything fine
			
				if str(entries[3]) == 'inversion':

					if str(entries[4]) != 'None':

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be None')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be None')
						sys.exit(1)

					if str(entries[0]) not in hap2dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

					else:

						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


				elif str(entries[3]) == 'deletion':

					if str(entries[4]) != 'None':

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be None')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be None')
						sys.exit(1)

					if str(entries[0]) not in hap2dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

					else:

						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


				elif str(entries[3]) == 'insertion':

					if not (all(i in valid_dna for i in entries[4].upper())): #validate sequence

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a DNA sequence with only A,C,T,G,N chars')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a DNA sequence with A,C,T,G,N chars')
						sys.exit(1)
							
					if str(entries[0]) not in hap2dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]).upper())]

					else:

						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]).upper()))

				elif str(entries[3]) == 'duplication':

					try:

						int(entries[4])

					except:

						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with number of duplication. Number must be an integer')
						sys.exit(1)

					if str(entries[0]) not in hap1dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4]))]

					else:

						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), int(entries[4])))


				elif str(entries[3]) == 'ptr': #perfect tandem repetition

					entr_4 = re.split('[:]',entries[4])

					if len(entr_4) != 2:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
						sys.exit(1)


					elif not (all(i in valid_dna for i in entr_4[0].upper())): #check if motif is vaild DNA sequence

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
						sys.exit(1)

					try:

						int(entr_4[1])

					except:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
						sys.exit(1)

					if str(entries[0]) not in hap2dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

					else:

						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))


				elif str(entries[3]) == 'atr': #approximative tandem repetition

					entr_4 = re.split('[:]',entries[4])

					if len(entr_4) != 3:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum')
						sys.exit(1)

					elif not (all(i in valid_dna for i in entr_4[0].upper())):

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Motif must be a valid DNA motif')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Motif must be a valid DNA motif')
						sys.exit(1)

					try:

						int(entr_4[1])

					except:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Number must be an integer')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Number must be an integer')
						sys.exit(1)


					try:

						int(entr_4[2])

					except:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Altnum must be an integer')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number:altnum. Altnum must be an integer')
						sys.exit(1)


					if str(entries[0]) not in hap2dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

					else:

						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))



				elif str(entries[3]) == 'tr expansion' or str(entries[3]) == 'tr contraction':

					entr_4 = re.split('[:]',entries[4])

					if len(entr_4) != 2:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number')
						sys.exit(1)

					elif not (all(i in valid_dna for i in entr_4[0].upper())):

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Motif must be a valid DNA motif')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number.  Motif must be a valid DNA motif')
						sys.exit(1)

					try:

						int(entr_4[1])

					except:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with motif:number. Number must be an integer')
						sys.exit(1)


					if str(entries[0]) not in hap2dict:

						hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4]))]

					else:


						hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), str(entries[3]), str(entries[4])))
						

				elif str(entries[3]) == 'translocation cut-paste':
							

					entr_4 = re.split('[:]',str(entries[4]))

					if len(entr_4) != 4:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
						sys.exit(1)


					if str(entr_4[0]) not in ['h1','h2']:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
						sys.exit(1)

					if str(entr_4[1]) not in classic_chrs:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
						sys.exit(1)

					try:

						int(entr_4[2])

					except:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
						sys.exit(1)


					if str(entr_4[3]) not in ['forward', 'reverse']:

						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Orientation can be forward or reverse')
						sys.exit(1)


					if entr_4[0] == 'h2':

						if str(entries[0]) not in hap2dict:

							hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), 'deletion', 'None')]

						else:

							hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), 'deletion', 'None'))

						if str(entr_4[1]) not in hap2dict:

							if str(entr_4[3]) == 'forward':

								hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

							else: #orientation is reverse

								hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]


						else:

							if str(entr_4[3]) == 'forward':


								hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

							else: #orientation is reverse

								hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


					elif entr_4[0] == 'h1':

						if str(entries[0]) not in hap2dict:

							hap2dict[str(entries[0])] = [(int(entries[1]), int(entries[2]), 'deletion', 'None')]

						else:

							hap2dict[str(entries[0])].append((int(entries[1]), int(entries[2]), 'deletion', 'None'))

						if str(entr_4[1]) not in hap2dict:

							if str(entr_4[3]) == 'forward':

								hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

							else: #orientation is reverse

								hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else:

							if str(entr_4[3]) == 'forward':


								hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

							else: #orientation is reverse

								hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

				
				else: #is a translocation copy paste


					entr_4 = re.split('[:]',str(entries[4]))

					if len(entr_4) != 4:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation')
						sys.exit(1)


					if str(entr_4[0]) not in ['h1','h2']:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Haplotype must be h1 or h2')
						sys.exit(1)

					if str(entr_4[1]) not in classic_chrs:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Chromosomes are ch1-chr22, chrX, chrY and chrM')
						sys.exit(1)

					try:

						int(entr_4[2])

					except:

						#print('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Breakpoint must be an integer')
						sys.exit(1)

					if str(entr_4[3]) not in ['forward', 'reverse']:

						logging.error('Incorrect info ' + str(entries[4]) + ' in .bed for haplotype 2 for variant ' + str(entries[3]) + '. Must be a string with haplotype:chromosome:breakpoint:orientation Orientation can be forward or reverse')
						sys.exit(1)

					if entr_4[0] == 'h2':

						if str(entr_4[1]) not in hap2dict:

							if str(entr_4[3]) == 'forward':

								hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

							else: #orientation is reverse

								hap2dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else:


							if str(entr_4[3]) == 'forward':


								hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

							else: #orientation is reverse

								hap2dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


					elif entr_4[0] == 'h1':

						if str(entr_4[1]) not in hap2dict:

							if str(entr_4[3]) == 'forward':

								hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

							else: #orientation is reverse

								hap1dict[str(entr_4[1])] = [(int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq)]

						else:

							if str(entr_4[3]) == 'forward':


								hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'insertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))

							else: #orientation is reverse

								hap1dict[str(entr_4[1])].append((int(entr_4[2]), int(entr_4[2])+1, 'invinsertion', immutable_ref[str(entries[0])][int(entries[1])-1:int(entries[2])].seq))


	logging.info('Generating .fa file with variants for haplotype 1')
	
	ParseDict(classic_chrs, immutable_ref, hap1dict, os.path.abspath(args.output + '/h1.fa'))
	
	
	if len(hap2dict) != 0:
		
		logging.info('Generating .fa file with variants for haplotype 2')
		
		ParseDict(classic_chrs, immutable_ref, hap2dict, os.path.abspath(args.output + '/h2.fa'))



	end=timeit.default_timer()
	elapsed=(end-start)/60

	#print('Haplotypes generated in ' + str(elapsed) + ' minutes')
	logging.info('Haplotypes generated in ' + str(elapsed) + ' minutes')

	#print('Done')
	logging.info('Done')



class CustomFormat(HelpFormatter):

	def _format_action_invocation(self, action):

		if not action.option_strings:

			default = self._get_default_metavar_for_positional(action)
			metavar, = self._metavar_formatter(action, default)(1)
			
			return metavar

		else:

			parts = []

			if action.nargs == 0:

				parts.extend(action.option_strings)

			else:

				default = self._get_default_metavar_for_optional(action)
				args_string = self._format_args(action, default)
				
				for option_string in action.option_strings:

					parts.append(option_string)

				return '%s %s' % (', '.join(parts), args_string)

			return ', '.join(parts)

	def _get_default_metavar_for_optional(self, action):

		return action.dest.upper()





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
	lastbase=trseq[-1]
	exprep=trseq[1:-1] + motif*num

	new_seq=firstbase + exprep + lastbase

	return new_seq


def CTRTR(infofield, sequence, start, end): #contract tr

	info=re.split('[:]', infofield)
	motif,num=str(info[0]), int(info[1])

	trseq=sequence[start-1:end]

	firstbase=trseq[0]
	lastbase=trseq[-1]
	rep=trseq[1:-1]
	newind=len(motif)*num
	delrep=rep[newind:]
	new_seq=firstbase + delrep + lastbase

	return new_seq




def ParseDict(chromosomes, fasta, dictionary, output_fasta):

	trans = str.maketrans('ATGC', 'TACG')

	for chrs in chromosomes:

		chrom=fasta[chrs]
		seq=chrom[:len(chrom)].seq

		#skip regions for which a variant has already been inserted

		regions_seen=[]


		if chrs not in dictionary.keys(): #chromosome not there, write unchanged
						
			write_unmodified_chromosome(chrs, seq, output_fasta)

		else:


			alterations_list=dictionary[chrs]

			if not len(alterations_list) == 1: #else is already sorted, as it has length 1

				alterations_list=sorted(alterations_list, key=itemgetter(1,2))

			i=0

			while i < len(alterations_list):


				start,end,typ,info=alterations_list[i]

				if any(el[0] <=start <= el[1] or el[0] <= end <= el[1] for el in regions_seen): #skip region if start or end overlap another variant

					i+=1
					continue

				else: #region not seen, insert variant

					regions_seen.append((start,end))

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

						write_sequence_between(alt_seq, output_fasta)

					elif typ == 'duplication':

						write_sequence_between(seq[start-1:end]*info, output_fasta)


					elif typ == 'ptr': #perfect tandem repetition

						alt_seq= PTR(info, seq, start, end)

						write_sequence_between(alt_seq, output_fasta)


					elif typ == 'atr': # approximate tandem repetition


						alt_seq=ATR(info,seq,start,end)

						write_sequence_between(alt_seq, output_fasta)


					elif typ == 'tr expansion': #expand a tandem repetition that is already present. Start-end are supposed to be as the one in repetitions .bed from ucsc


						alt_seq=EXPTR(info,seq,start,end)

						write_sequence_between(alt_seq, output_fasta)

					elif typ == 'tr contraction': #contract a tandem repetition that is already present. Start-end are supposed to be as the one in repetitions .bed from ucsc

						alt_seq=CTRTR(info,seq,start,end)
				
						write_sequence_between(alt_seq, output_fasta)


					if i == len(alterations_list) -1:

						write_end_sequence(seq[end:], output_fasta) #end not included, as it was included in the variant


					elif i < len(alterations_list) -1:


						nextstart=alterations_list[i+1][0]
						thisend=end

						write_sequence_between(seq[thisend:nextstart-1], output_fasta)

					i+=1




if __name__ == '__main__':
   
	main()
