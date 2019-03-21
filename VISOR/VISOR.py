import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='HACk,SHORtS,LASeR') #two submodules


	## HACk ##

	parser_hack = subparsers.add_parser('HACk', help='HAplotype Creator. Generates 2 haplotypes in .fasta format containing variants specified in .bed files')


	required = parser_hack.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa', required=True)
	required.add_argument('-h1b', '--hap1bed', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 1', metavar='.bed', required=True)
	required.add_argument('-O', '--output', help='where the 2 .fa haplotypes will be saved', metavar='folder', required=True)

	optional = parser_hack.add_argument_group('Additional input')
	optional.add_argument('-h2b', '--hap2bed', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 2', metavar='', default=None)

	parser_hack.set_defaults(func=run_subtool)


	## SHORtS ##

	parser_shorts = subparsers.add_parser('SHORtS', help='SHOrt Reads Simulator. Simulate short reads .bam files from .fasta files using regions specified in .bed files. Simulations are run using wgsim')


	required = parser_shorts.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa', required=True)
	required.add_argument('-h1f','--hap1fa', help='.fasta file containing variants for haplotype 1', metavar='.fa', required=True)
	required.add_argument('-h2f','--hap2fa', help='.fasta file containing (or not) variants for haplotype 2', metavar='.fa', required=True)
	required.add_argument('-h1b','--hap1bed', help='.bed file containing regions to simulate for haplotype 1. To simulate an entire chromosome START must be 0 and END must be chromosome length.', metavar='.bed', required=True)
	required.add_argument('-h2b','--hap2bed', help='.bed file containing regions to simulate for haplotype 2. To simulate an entire chromosome START must be 0 and END must be chromosome length.', metavar='.bed', required=True)
	required.add_argument('-O','--output', help='where the simulated .bam files will be saved', metavar='folder', required=True)

	
	simtype = parser_shorts.add_argument_group('Type of simulation')

	simtype.add_argument('-t','--type', help='Whether to simulate double-strand or single-strand (strand-seq) short-reads .bam files. If simulating strand-seq short-read .bam files, for each haplotype, 2 will be created, named watson (R1 forward, R2 reverse) and crick (R1 reverse, R2 forward) [double-strand]', metavar='', default='double-strand', choices=['single-strand', 'double-strand'])

	wgi = parser_shorts.add_argument_group('Wgsim parameters for simulation')

	wgi.add_argument('-e', '--error', help='Base error rate [0.010]', metavar='', default=0.010, type=float)
	wgi.add_argument('-c', '--coverage', help='Desired coverage for the simulated region [30]', metavar='', default=30, type=int)
	wgi.add_argument('-l', '--length', help='Length of reads [150]', metavar='', default=150, type=int)
	wgi.add_argument('-i', '--indels', help='Fractions of indels [0.000000001]', metavar='', default=0.000000001, type=float)
	wgi.add_argument('-p', '--probability', help='Probability an indel is extended [0.000000001]', metavar='', default=0.000000001, type=float)


	optional = parser_shorts.add_argument_group('Additional parameters')

	optional.add_argument('-th', '--threads', help='Number of cores to use for alignments [6]', metavar='', type=int, default=6)
	optional.add_argument('-n', '--noise', help='percentage of noise to add to single-strand .bam files [0.00]', type=float, metavar='', default=0.00)

	parser_shorts.set_defaults(func=run_subtool)


	## LASeR ##


	parser_long = subparsers.add_parser('LASeR', help='Long reAds SimulatoR. Simulate long reads .bam files from .fasta files using regions specified in .bed files. Simulations are run using pbsim')


	required = parser_long.add_argument_group('Required I/O arguments')

	required.add_argument('-h1f','--hap1fa', help='.fasta file containing variants for haplotype 1', metavar='.fa', required=True)
	required.add_argument('-h2f','--hap2fa', help='.fasta file containing (or not) variants for haplotype 2', metavar='.fa', required=True)
	required.add_argument('-h1b','--hap1bed', help='.bed file containing regions to simulate for haplotype 1. To simulate an entire chromosome START must be 0 and END must be chromosome length.', metavar='.bed', required=True)
	required.add_argument('-h2b','--hap2bed', help='.bed file containing regions to simulate for haplotype 2. To simulate an entire chromosome START must be 0 and END must be chromosome length.', metavar='.bed', required=True)
	required.add_argument('-O','--output', help='where the simulated .bam files will be saved', metavar='folder', required=True)

	pbs= parser_long.add_argument_group('Pbsim parameters for simulation')


	pbs.add_argument('-a', '--accuracy', help='mean accuracy for simulated reads [0.90]', metavar='', default=0.90, type=float)
	pbs.add_argument('-l', '--length', help='mean length for simulated reads [8000]', metavar='', default=8000, type=int)
	pbs.add_argument('-c', '--coverage', help='mean coverage for the simulated region [20]', metavar='', default=20, type=int)
	pbs.add_argument('-r', '--ratio', help='substitution:insertion:deletion ratio [30:30:40]', metavar='', default='30:30:40', type=str)

	optional = parser_long.add_argument_group('Additional parameter')

	optional.add_argument('-th', '--threads', help='Number of cores to use for alignments [6]', metavar='', type=int, default=6)



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
