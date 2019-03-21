import argparse
from argparse import HelpFormatter



def main():

	parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This program was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

	subparsers = parser.add_subparsers(title='modules', dest='command', metavar='HACk, SHORtS,LASeR') #two submodules


	## HACk ##

	parser_hack = subparsers.add_parser('HACk', help='HAplotype Creator. Generates 2 haplotypes in .fasta format containing variants specified in .bed files')


	required = parser_hack.add_argument_group('Required I/O arguments')

	required.add_argument('-g','--genome', help='reference genome', metavar='.fa', required=True)
	required.add_argument('-h1b', '--haplotype1bed', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 1', metavar='.bed', required=True)
	required.add_argument('-O', '--output', help='name of the directory where the 2 .fa haplotypes will be saved', metavar='folder', required=True)

	optional = parser_hack.add_argument_group('Additional input')
	optional.add_argument('-h2b', '--haplotype2bed', help='.bed file containing "CHROM, START, END, ALT, INFO" for haplotype 2', metavar='', default=None)

	parser_hack.set_defaults(func=run_subtool)


	## SHORtS ##

	#parser_shorts = subparsers.add_parser('SHORtS', help='SHOrt Reads Simulator. Simulate short reads .bam files from .fasta files using regions specified in .bed files')


	#required = parser_shorts.add_argument_group('Required I/O arguments')


	#parser_shorts.set_defaults(func=run_subtool)




	## LASeR ##


	#parser_long = subparsers.add_parser('LASeR', help='Long reAds SimulatoR. Simulate long reads .bam files from .fasta files using regions specified in .bed files')


	#required = parser_shorts.add_argument_group('Required I/O arguments')


	#parser_shorts.set_defaults(func=run_subtool)




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
