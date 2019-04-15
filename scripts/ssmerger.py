import argparse
from argparse import HelpFormatter
import sys
import os
import glob
import subprocess


def main():


	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This script is included in VISOR and was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	required = parser.add_argument_group('Required I/O arguments')

	required.add_argument('-f', '--folder', help='One or more folders containing single-strand .bam files generated with VISOR SHORtS', metavar='folder', nargs='+', action='append', required=True)
	required.add_argument('-s', '--strand', help='One ore more ordered acronyms for strands to merge. For example, assuming 2 input folders given (one for each haplotype in a diploid sample), using "-s C W" generates a "CW.srt.bam" file merging crick strand from first and watson from second. Choices are W and C', metavar='id', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', help='Output folder', metavar='folder', required=True)

	args = parser.parse_args()

	folders=args.folder[0]
	strands = args.strand[0]

	if len(folders) != len(strands):

		print ('Different number of inputs for -f/--folder and -s/--sample')
		sys.exit(1)

	bams =[]

	for stra,fol in zip(strands, folders):

		if len(stra) != 1 or stra not in ['W','C']:

			print('Wrong acronym ' + str(stra) + '.')
			sys.exit(1)

		if stra == 'W': #look for a watson.srt.bam

			W = glo.glob(os.path.abspath(fol) + '/*.watson.srt.bam')

			if W == []:

				print('No watson .bam found in ' + os.path.abspath(fol))
				sys.exit(1)

			else:

				bams.extend(W)

		else:

			C = glo.glob(os.path.abspath(fol) + '/*.crick.srt.bam')

			if C == []:

				print('No crick .bam found in ' + os.path.abspath(fol))
				sys.exit(1)

			else:

				bams.extend(C)
			
	with open(os.path.abspath(args.output + '/bamstomerge.txt'), 'w') as bamstomerge:

		for bam in bams:

			bamstomerge.write(bam + '\n')

	identifier = ''.join(x for x in strands)

	subprocess.call(['samtools', 'merge', '-b', os.path.abspath(args.output + '/bamstomerge.txt'), os.path.abspath(args.output + '/' + identifier + '.srt.bam')], stderr=open(os.devnull, 'wb'))
	subprocess.call(['samtools', 'index', os.path.abspath(args.output + '/' + args.identifier + '.srt.bam')], stderr=open(os.devnull, 'wb'))

	os.remove(os.path.abspath(args.output + '/bamstomerge.txt'))

if __name__ == '__main__':

	main()
