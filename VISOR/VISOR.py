#!/usr/bin/python env

import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI). Extensive documentation is available at: https://davidebolo1993.github.io/visordoc/''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='HACk,SHORtS,LASeR,LIKER') #two submodules

	## HACk ##

	parser_hack = subparsers.add_parser('HACk', help='HAplotype Creator. Generates one or more haplotypes in FASTA format containing SVs specified in BED file/s.')


	required = parser_hack.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='Template reference genome', metavar='FASTA', required=True)
	required.add_argument('-bed', '--bedfile', help='One or more BED files (one for each haplotype) containing "CHROM, START, END, ALT, INFO, BREAKSEQLEN" entries for each SV', metavar='BED', nargs='+', action='append', required=True)
	required.add_argument('-o', '--output', help='Output folder', metavar='FOLDER', required=True)

	parser_hack.set_defaults(func=run_subtool)


	## SHORtS ##

	parser_shorts = subparsers.add_parser('SHORtS', help='SHOrt Reads Simulator. Simulate short reads BAM files from FASTA haplotypes using regions specified in BED file.')

	required = parser_shorts.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='Template reference genome', metavar='FASTA', required=True)
	required.add_argument('-s','--sample', help='One or more folders containing FASTA haplotypes with SVs generated with VISOR HACk. If multiple folders are given, each sample is considered a subclone', metavar='FOLDER',  nargs='+', action='append', required=True)
	required.add_argument('-bed','--bedfile', help='BED file containing one or more "CHROM, START, END, CAPTURE BIAS, SAMPLE FRACTION" for regions to simulate. CAPTURE BIAS and SAMPLE FRACTION must be float pecentages', metavar='BED', required=True)
	required.add_argument('-o','--output', help='Output folder', metavar='FOLDER', required=True)
	
	simtype = parser_shorts.add_argument_group('Type of simulation')

	simtype.add_argument('-t','--type', help='Whether to simulate bulk [bulk] or strand-seq [strand-seq] data [bulk].', metavar='', default='bulk', choices=['bulk', 'strand-seq'])

	wgi = parser_shorts.add_argument_group('Wgsim parameters for FASTQ simulations')

	wgi.add_argument('-c', '--coverage', help='Mean coverage for the simulated region [30.0]', metavar='', default=30.0, type=float)
	wgi.add_argument('-e', '--error', help='Base error rate [0.010]', metavar='', default=0.010, type=float)
	wgi.add_argument('-l', '--length', help='Length of reads [150]', metavar='', default=150, type=int)
	wgi.add_argument('-i', '--indels', help='Fractions of indels [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-p', '--probability', help='Probability an indel is extended [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-is', '--insertsize', help='0uter distance between the two ends [500]',metavar='', default=500, type=int)
	wgi.add_argument('-sd', '--standardev', help='Standard deviation for insert size [50]',metavar='', default=50, type=int)

	bulk = parser_shorts.add_argument_group('Subclones parameters for bulk data')

	bulk.add_argument('--clonefraction', help='Ordered percentages for each clone specified in -s/--sample [None]', metavar='', nargs='+', action='append', default=None)

	strandseq = parser_shorts.add_argument_group('Strand-seq parameters')

	strandseq.add_argument('--scebedfile', help='BED file containing "CHROM, START, END, HAPLOTYPE" in which sister chromatid exchange will be performed. If a BED is given, HAPLOTYPE must be in format "hN" where N is the number of the haplotype [None]', metavar='', default=None)
	strandseq.add_argument('--noise', help='Percentage of noise to add to strand-seq BAM files [0.00]', type=float, metavar='', default=0.00)
	
	optional = parser_shorts.add_argument_group('Additional general parameters')
		
	optional.add_argument('--threads', help='Number of cores to use for alignments [1]', metavar='', type=int, default=1)
	optional.add_argument('--identifier', help='Identifier to label the output [sim]', metavar='', default='sim')
	optional.add_argument('--noaddtag', help='Do not tag reads in BAM by haplotype and clone number. Reads in strand-seq and 10X linked reads data are not tagged by default', action='store_false')	

	parser_shorts.set_defaults(func=run_subtool)


	## LASeR ##


	parser_long = subparsers.add_parser('LASeR', help='Long reAds SimulatoR. Simulate long reads BAM files from FASTA files using regions specified in BED file.')


	required = parser_long.add_argument_group('Required I/O arguments')
	
	required.add_argument('-g','--genome', help='Template reference genome', metavar='FASTA', required=True)
	required.add_argument('-s','--sample', help='One or more folders containing FASTA haplotypes with SVs generated with VISOR HACk. If multiple folders are given, each sample is considered a subclone', metavar='FOLDER',  nargs='+', action='append', required=True)
	required.add_argument('-bed','--bedfile', help='BED file containing one or more "CHROM, START, END, CAPTURE BIAS, SAMPLE FRACTION" for regions to simulate. CAPTURE BIAS and SAMPLE FRACTION must be float pecentages', metavar='BED', required=True)
	required.add_argument('-o','--output', help='Output folder', metavar='FOLDER', required=True)

	pbs= parser_long.add_argument_group('Pbsim parameters for FASTQ simulations')

	pbs.add_argument('-c', '--coverage', help='Mean coverage for the simulated region [20]', metavar='', default=20.0, type=float)
	pbs.add_argument('-a', '--accuracy', help='Mean accuracy for simulated reads [0.90]', metavar='', default=0.90, type=float)
	pbs.add_argument('-l', '--length', help='Mean length for simulated reads [8000]', metavar='', default=8000, type=int)
	pbs.add_argument('-r', '--ratio', help='substitution:insertion:deletion ratio [30:30:40]', metavar='', default='30:30:40', type=str)

	bulk = parser_long.add_argument_group('Subclones parameters')

	bulk.add_argument('--clonefraction', help='Ordered percentages for each clone specified in -s/--sample [None]', metavar='', nargs='+', action='append', default=None)

	optional = parser_long.add_argument_group('Additional general parameters')

	optional.add_argument('--threads', help='Number of cores to use for alignments [1]', metavar='', type=int, default=1)
	optional.add_argument('--identifier', help='Identifier to label the output [sim]', metavar='', default='sim')
	optional.add_argument('--noaddtag', help='Do not tag reads in BAM by haplotype and clone number', action='store_false')	

	parser_long.set_defaults(func=run_subtool)


	## XENIA ## [Beta version]

	parser_tenx = subparsers.add_parser('XENIA', help='10X gENomics sImulAtor. Simulate 10X Genomics FASTQ files (linked-reads FASTQ or FASTQ from single-cells for CNV detection)  from FASTA files using regions specified in BED file. Please note that this module is released in BETA version.')


	required = parser_tenx.add_argument_group('Required I/O arguments')

	required.add_argument('-s','--sample', help='A folder containing FASTA haplotypes with SVs generated with VISOR HACk', metavar='FOLDER', required=True)
	required.add_argument('-bed','--bedfile', help='BED file containing one or more "CHROM, START, END" for regions to simulate. BED for VISOR SHORtS and LASeR are accepted, but CAPTURE BIAS and SAMPLE FRACTION are ignored', metavar='BED', required=True)
	required.add_argument('-o','--output', help='Output folder', metavar='FOLDER', required=True)

	wgi = parser_tenx.add_argument_group('Wgsim parameters for FASTQ simulations')

	wgi.add_argument('-c', '--coverage', help='Mean coverage for the simulated region [30.0]', metavar='', default=30.0, type=float)
	wgi.add_argument('-e', '--error', help='Base error rate [0.010]', metavar='', default=0.010, type=float)
	wgi.add_argument('-l', '--length', help='Length of reads [150]', metavar='', default=150, type=int)
	wgi.add_argument('-i', '--indels', help='Fractions of indels [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-p', '--probability', help='Probability an indel is extended [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-is', '--insertsize', help='0uter distance between the two ends [500]',metavar='', default=500, type=int)
	wgi.add_argument('-sd', '--standardev', help='Standard deviation for insert size [50]',metavar='', default=50, type=int)

	simtype = parser_tenx.add_argument_group('Type of simulation')

	simtype.add_argument('-t','--type', help='Whether to simulate bulk [bulk] or single-cell [single-cell] data [bulk].', metavar='', default='bulk', choices=['bulk', 'single-cell'])

	molecules=parser_tenx.add_argument_group('10X linked reads')

	molecules.add_argument('--molecules_length', help='Mean molecules length [80000]', default=80000, type=int, metavar='')
	molecules.add_argument('--molecules_number', help='Mean number of molecules per GEM [10]', default=10, type=int, metavar='')
	molecules.add_argument('--molecules_coverage', help='Mean numbercoverage per molecule [0.2]', default=0.2, type=float, metavar='')

	singlecell=parser_tenx.add_argument_group('10X barcoded reads for single cell CNV detection applications')

	singlecell.add_argument('--cells_number', help='Number of cells to simulate', default=100, type=int, metavar='')

	optional = parser_tenx.add_argument_group('Additional general parameters')

	optional.add_argument('--threads', help='Number of cores to use for bulk FASTQ parallel simulations [1]', metavar='', type=int, default=1)
	optional.add_argument('--identifier', help='Identifier to label the output [sim]', metavar='', default='sim')

	parser_tenx.set_defaults(func=run_subtool)

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

	elif args.command == 'XENIA':

		from .XENIA import XENIA as submodule

	else:

		parser.print_help()

	submodule.run(parser,args)


if __name__ =='__main__':

	main()
