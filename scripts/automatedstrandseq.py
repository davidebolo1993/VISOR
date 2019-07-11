#!/usr/bin/python env

import argparse
from argparse import HelpFormatter
import subprocess, shlex
import os
import sys
import random
from shutil import which


def main():

	parser = argparse.ArgumentParser(prog='VISOR script', description='''Automated simulations of single cell, strand-seq data using VISOR SHORtS''', epilog='''This script has been developed by Davide Bolognini at EMBL/EBI and is provided together with VISOR''', formatter_class=CustomFormat) 

	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-as', '--affectedsample', help='folder containing one or more .fasta haplotypes with variants generated with VISOR HACk to use as templates for affected cells', metavar='folder')
	required.add_argument('-ns', '--normalsample', help='folder containing one or more .fasta haplotypes without variants to use as templates for normal cells, if normal cells have to be simulated', metavar='folder')
	required.add_argument('-g', '--genome', help='reference genome', metavar='fasta', required=True)
	required.add_argument('-sb', '--shortsbed', help='.bed file for VISOR SHORtS with regions to simulate', metavar='.bed', required=True)
	required.add_argument('-mp', '--mergerpath', help='path to ssmerger.py script provided into VISOR script folder', metavar='pyscript', required=True)
	required.add_argument('-o', '--output', help='output folder', metavar='folder', required=True)

	optional = parser.add_argument_group('Optional parameters')

	optional.add_argument('-nc', '--numberofcells', help='number of cells to simulate [100]',type=int, default=100)
	optional.add_argument('-af', '--affected', help='ratio of cells having the alterations [80.0]', type=float, default=80.0)
	optional.add_argument('-rl', '--readslength', help='short reads length [150]', type=int, default=150)
	optional.add_argument('-e', '--error', help='short reads error [0.005]', type=float, default=0.005)
	optional.add_argument('-c', '--coverage', help='coverage for simulated data [0.04]', type=float, default=0.04)
	optional.add_argument('--noise', help='percentage of noise (reads with incorrect orientation) to add to srand-seq simulations [5.0]', type=float, default=5.0)
	optional.add_argument('--scebed', help='.bed file for VISOR SHORtS describing sister chromatid exchange events in affected cells [None]', default=None)	
	optional.add_argument('-th', '--threads', help='number of cores to use for alignments [1]', type=int, default=1)

	args = parser.parse_args()

	if not os.path.exists(os.path.abspath(args.output)):
		try:

			os.makedirs(os.path.abspath(args.output))

		except:

			print('Cannot create the output folder') 	
			sys.exit(1)

	else: #path already exists

		if not os.access(os.path.abspath(args.output),os.W_OK):

			print('Missing write permissions on the output folder')			
			sys.exit(1)
			
		elif os.listdir(os.path.abspath(args.output)):

			print('The output folder is not empty. Specify another output folder or clean the previsouly chosen')
			sys.exit(1)


	external_tools=['VISOR']

	for tools in external_tools: 

		if which(tools) is None:

			print(tools + ' cannot be executed. Install ' + tools + ' and re-run this script')
			sys.exit(1)

	try:

		with open(os.path.abspath(args.genome),'r') as file:

			assert(file.readline().startswith('>'))

	except:

		print('Reference file does not exist, is not readable or is not a valid FASTA')
		sys.exit(1)


	if not os.path.exists(os.path.abspath(args.shortsbed)):

		print('The .bed file for VISOR SHORtS ' + os.path.abspath(args.shortsbed) + ' does not exist')
		sys.exit(1)


	if not os.path.exists(os.path.abspath(args.mergerpath)):

		print('Specified path ' + os.path.abspath(args.mergerpath) + ' for ssmerger.py script does not exists')
		sys.exit(1)

	if args.scebed is not None:

		if not os.path.exists(os.path.abspath(args.scebed)):

			print('Specified path ' + os.path.abspath(args.scebed) + ' for --scebed .bed file does not exist')
			sys.exit(1)


	total=args.numberofcells
	numberaffected=round((total/100)*args.affected)
	numbernotaffected=total-numberaffected

	possiblestrands=['W','C']

	if numberaffected == 0:

		print('No cells are affected, simulate only normal cells')

	else:

		print(str(numberaffected) + ' affected cells to simulate. Simulate these first...')

		done=False
		i=1

		while not done:

			print('Simulating normal cell ' + str(i))

			if args.scebed is None:

				subprocess.call(['VISOR', 'SHORtS', '-g', os.path.abspath(args.genome), '-l', str(args.readslength), '-e', str(args.error), '-bed', os.path.abspath(args.shortsbed), '-t', 'single-strand', '-s', os.path.abspath(args.normalsample), '-o', os.path.abspath(args.output + '/cell' + str(i) + '.affected'),'-th', str(args.threads), '-n', str(args.noise), '-c', str(args.coverage)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			else:

				subprocess.call(['VISOR', 'SHORtS', '-g', os.path.abspath(args.genome), '-l', str(args.length), '-e', str(args.error), '-bed', os.path.abspath(args.shortsbed), '-t', 'single-strand', '-s', os.path.abspath(args.normalsample), '-o', os.path.abspath(args.output + '/cell' + str(i) + '.affected'),'-th', str(args.threads), '-n', str(args.noise), '-c', str(args.coverage), '-scebed', os.path.abspath(args.scebed)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			subdirs=[f.path for f in os.scandir(os.path.abspath(args.output + '/cell' + str(i) + '.affected')) if f.is_dir()]
			eachhaplostrand=[]
			
			for subs in subdirs:

				eachhaplostrand.append(random.choice(possiblestrands))

			folders=' '.join(os.path.abspath(x) for x in subdirs)
			haplos=' '.join(x for x in eachhaplostrand)

			command='python ' + os.path.abspath(args.mergerpath) + ' -f ' + folders + ' -s ' + haplos + ' -o ' + os.path.abspath(args.output + '/cell' + str(i) + '.affected/BAM')
			print(shlex.split(command))
			subprocess.call(shlex.split(command), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			i+=1

			if i>numberaffected:

				done=True

	if numbernotaffected==0:

		print('No normal cells to simulate. Done')

	else:

		print(str(numbernotaffected) + ' not affected cells to simulate. Simulate these then...')

		done=False

		while not done:

			print('Simulating normal cell ' + str(i))
			
			subprocess.call(['VISOR', 'SHORtS', '-g', os.path.abspath(args.genome), '-l', str(args.readslength), '-e', str(args.error), '-bed', os.path.abspath(args.shortsbed), '-t', 'single-strand', '-s', os.path.abspath(args.normalsample), '-o', os.path.abspath(args.output + '/cell' + str(i) + '.normal'),'-th', str(args.threads), '-n', str(args.noise), '-c', str(args.coverage)],stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			subdirs=[f.path for f in os.scandir(os.path.abspath(args.output + '/cell' + str(i) + '.normal')) if f.is_dir()]
			eachhaplostrand=[]
			
			for subs in subdirs:

				eachhaplostrand.append(random.choice(possiblestrands))

			folders=' '.join(os.path.abspath(x) for x in subdirs)
			haplos=' '.join(x for x in eachhaplostrand)

			command='python ' + os.path.abspath(args.mergerpath) + ' -f ' + folders + ' -s ' + haplos + ' -o ' + os.path.abspath(args.output + '/cell' + str(i) + '.normal/BAM')
			print(shlex.split(command))
			subprocess.call(shlex.split(command), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

			i+=1

			if i>total:

				done=True

	print('Done')




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




if __name__ == '__main__':

	main()
