#!/usr/bin/python env

import argparse
from argparse import HelpFormatter


def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='HACk, SCoRE') #two submodules

	## HACk ##

	parser_gen = subparsers.add_parser('HACk', help='HAplotypes Creator. Generate 2 haplotypes with user-defined variants; outputs are in fasta (.fa) format')

	required = parser_gen.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa', required=True)
	required.add_argument('-h1b', '--haplotype1bed', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 1', metavar='.bed', required=True)
	required.add_argument('-O', '--output', help='name of the directory where the 2 .fa haplotypes will be saved', metavar='folder', required=True)

	optional = parser_gen.add_argument_group('Additional input')
	optional.add_argument('-h2b', '--haplotype2bed', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 2', metavar='', default=None)

	parser_gen.set_defaults(func=run_subtool)

	## SCoRE ##

	parser_sim = subparsers.add_parser('SCoRE', help='Simulations CREator. Simulate .bam files with reads containing alterations from .fa files')

	required = parser_sim.add_argument_group('Required I/O arguments')
	required.add_argument('-h1fa','--haplotype1fa', help='haplotype 1 .fa file contaning alterations', metavar='.fa', required=True)
	required.add_argument('-h1b','--haplotype1bed', help='haplotype 1 .bed file containing regions to simulate with chromosome, start, end, label', metavar='.fa',required=True)
	required.add_argument('-O', '--output', help='name of the directory where the simulated .bam files will be saved', metavar='folder', required=True)



	mod = parser_sim.add_argument_group('Simulations parameters')

	mod.add_argument('-m','--mode', help='Type of simulation [classic short-read]', default='classic short-read', choices=['classic short-read', 'strand-sequencing short-read', 'classic long-read'],metavar='')

	short = parser_sim.add_argument_group('Short reads main parameters for wgsim')

	short.add_argument('-er','--errorrate', help='Base error rate [0.010]', default=0.010, type=float,metavar='')
	short.add_argument('-c','--coverage', help='Mean desired coverage [30]', default=30, type=int,metavar='')
	short.add_argument('-l1','--length1', help='Length (#bps) of first read in pair [150]', default=150, type=int,metavar='')
	short.add_argument('-l2','--length2', help='Length (#bps) of second read in pair [150]', default=150, type=int,metavar='')
	short.add_argument('-ind','--indels', help='fraction of indels [0.000000001]', default=0.000000001, type=float,metavar='')
	short.add_argument('-indext','--indelsextension', help='probability an indel is extended [0.000000001]', default=0.000000001, type=float,metavar='')


	lon = parser_sim.add_argument_group('Long reads main parameters for pbsim')

	lon.add_argument('-ml','--meanlength', help='Mean length (#bps) for reads [8000]', default=8000, type=int,metavar='')
	lon.add_argument('-ma','--meanaccuracy', help='Mean accuracy (#bps) for reads [0.90]', default=0.90, type=float,metavar='')
	lon.add_argument('-dr', '--differenceratio', help='Ratio of substitutions:insertions:deletions[30:30:40]', default='30:30:40', type=str, metavar='')
	lon.add_argument('-mc','--meancoverage', help='Mean coverage [20]', default=0.90, type=float,metavar='')


	optional = parser_sim.add_argument_group('Additional inputs')

	optional.add_argument('-h2fa','--haplotype2fa', help='haplotype 2 .fa file containing (or not) alterations', metavar='')
	optional.add_argument('-h2b','--haplotype2bed', help='haplotype 2 .bed file containing regions to simulate with chromosome, start, end, label', metavar='')

	parser_sim.set_defaults(func=run_subtool)


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
	
	elif args.command == 'SCoRE': 

		from .SCoRE import SCoRE as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':

	main()
