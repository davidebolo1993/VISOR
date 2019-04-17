#!/usr/bin/python env

import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='HACk,SHORtS,LASeR') #two submodules

	## HACk ##

	parser_hack = subparsers.add_parser('HACk', help='HAplotype Creator. Generates one or more haplotypes in .fasta format containing SVs specified in .bed file/s')


	required = parser_hack.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='Template reference genome', metavar='.fa', required=True)
	required.add_argument('-bed', '--bedfile', help='One or more .bed files (one for each haplotype) containing "CHROM, START, END, ALT, INFO" entries for each SV', metavar='.bed', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', help='Output folder', metavar='folder', required=True)

	parser_hack.set_defaults(func=run_subtool)


	## SHORtS ##

	parser_shorts = subparsers.add_parser('SHORtS', help='SHOrt Reads Simulator. Simulate short reads .bam files from .fasta files using regions specified in .bed file.')

	required = parser_shorts.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='Template reference genome', metavar='.fa', required=True)
	required.add_argument('-s','--sample', help='One or more folders containing .fasta haplotypes with SVs generated with VISOR HACk. If multiple folders are given, each sample is considered a subclone', metavar='folder/s',  nargs='+', action='append', required=True)
	required.add_argument('-bed','--bedfile', help='.bed file containing one or more "CHROM, START, END, CAPTURE BIAS, SAMPLE FRACTION" for regions to simulate. CAPTURE BIAS and SAMPLE FRACTION must be float pecentages', metavar='.bed', required=True)
	required.add_argument('-o','--output', help='Output folder', metavar='folder', required=True)

	
	simtype = parser_shorts.add_argument_group('Type of simulations')

	simtype.add_argument('-t','--type', help='Whether to simulate double-strand or single-strand (strand-seq) short-reads. [double-strand]', metavar='', default='double-strand', choices=['single-strand', 'double-strand'])

	wgi = parser_shorts.add_argument_group('Wgsim parameters for .fastq simulations')

	wgi.add_argument('-e', '--error', help='Base error rate [0.010]', metavar='', default=0.010, type=float)
	wgi.add_argument('-c', '--coverage', help='Mean coverage for the simulated region [30.0]', metavar='', default=30.0, type=float)
	wgi.add_argument('-l', '--length', help='Length of reads [150]', metavar='', default=150, type=int)
	wgi.add_argument('-i', '--indels', help='Fractions of indels [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-p', '--probability', help='Probability an indel is extended [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-is', '--insertsize', help='0uter distance between the two ends [500]',metavar='', default=500, type=int)
	wgi.add_argument('-sd', '--standardev', help='Standard deviation for insert size [50]',metavar='', default=50, type=int)


	optional = parser_shorts.add_argument_group('Single-strand parameters')

	optional.add_argument('-scebed', '--scebedfile', help='.bed file containing "CHROM, START, END, HAPLOTYPE" in which sister chromatid exchange will be performed. If a .bed is given, HAPLOTYPE must be in format "hN" where N is the number of the haplotype [None]', metavar='', default=None)
	optional.add_argument('-n', '--noise', help='Percentage of noise to add to single-strand .bam files [0.00]', type=float, metavar='', default=0.00)
	
	optional1 = parser_shorts.add_argument_group('Subclones simulations')

	optional1.add_argument('-cf', '--clonefraction', help='Ordered percentages for each clone specified in -s/--sample [None]', metavar='', nargs='+', action='append', default=None)
	
	optional2 = parser_shorts.add_argument_group('Additional general parameters')
		
	optional2.add_argument('-th', '--threads', help='Number of cores to use for alignments [7]', metavar='', type=int, default=7)
	optional2.add_argument('-id', '--identifier', help='Identifier to label the output [sim]', metavar='', default='sim')
		
	parser_shorts.set_defaults(func=run_subtool)


	## LASeR ##


	parser_long = subparsers.add_parser('LASeR', help='Long reAds SimulatoR. Simulate long reads .bam files from .fasta files using regions specified in .bed files. Simulations are run using pbsim.')


	required = parser_long.add_argument_group('Required I/O arguments')
	
	required.add_argument('-g','--genome', help='Template reference genome', metavar='.fa', required=True)
	required.add_argument('-s','--sample', help='One or more folders containing .fasta haplotypes with SVs generated with VISOR HACk. If multiple folders are given, each sample is considered a subclone', metavar='folder/s',  nargs='+', action='append', required=True)
	required.add_argument('-bed','--bedfile', help='.bed file containing one or more "CHROM, START, END, CAPTURE BIAS, SAMPLE FRACTION" for regions to simulate. CAPTURE BIAS and SAMPLE FRACTION must be float pecentages', metavar='.bed', required=True)
	required.add_argument('-o','--output', help='Output folder', metavar='folder', required=True)

	pbs= parser_long.add_argument_group('Pbsim parameters for .fastq simulations')

	pbs.add_argument('-a', '--accuracy', help='Mean accuracy for simulated reads [0.90]', metavar='', default=0.90, type=float)
	pbs.add_argument('-l', '--length', help='Mean length for simulated reads [8000]', metavar='', default=8000, type=int)
	pbs.add_argument('-c', '--coverage', help='Mean coverage for the simulated region [20]', metavar='', default=20.0, type=float)
	pbs.add_argument('-r', '--ratio', help='substitution:insertion:deletion ratio [30:30:40]', metavar='', default='30:30:40', type=str)

	optional = parser_long.add_argument_group('Subclones simulations')

	optional.add_argument('-cf', '--clonefraction', help='Ordered percentages for each clone specified in -s/--sample [None]', metavar='', nargs='+', action='append', default=None)

	optional1 = parser_long.add_argument_group('Additional general parameters')

	optional1.add_argument('-th', '--threads', help='Number of cores to use for alignments [7]', metavar='', type=int, default=7)
	optional1.add_argument('-id', '--identifier', help='Identifier to label the output [sim]', metavar='', default='sim')

	parser_long.set_defaults(func=run_subtool)

	args = parser.parse_args()
	args.func(parser, args)



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



def run_subtool(parser, args):

	if args.command == 'HACk': 

		from .HACk import HACk as submodule
	
	elif args.command == 'SHORtS': 

		from .SHORtS import SHORtS as submodule

	elif args.command == 'LASeR':

		from .LASeR import LASeR as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':

	main()
