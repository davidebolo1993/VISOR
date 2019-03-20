#!/usr/bin/python env

from collections import defaultdict
import argparse
from argparse import HelpFormatter
import subprocess
from shutil import which

#additional library

import pysam



def main():


    parser = argparse.ArgumentParser(prog='VISOR', description='''VarIants SimulatOR''', epilog='''This script from VISOR was developed by Davide Bolognini at the European Molecular Biology Laboratory/European Bioinformatic Institute (EMBL/EBI)''', formatter_class=CustomFormat) 

    required = parser.add_argument_group('Required I/O arguments')

    required.add_argument('-bam','--bamfile', help='.bam file generated using simulated reads to convert into strand-seq .bam', metavar='.bam', required=True)
    required.add_argument('-O', '--output', help='name of the directory where the .bam file will be saved', metavar='folder', required=True)


    optional = parser.add_argument_group('Strand-seq type')

    optional.add_argument('-t','--type', help='define orientations of reads 1 (R1) and reads 2 (R2) in .bam file. If "watson", all R1 are forward and all R2 are reverse, if "crick", all R1 are reverse and all R2 are forward [watson]', metavar='', default='watson', choices=['crick', 'watson'])

    args = parser.parse_args()

    tools = 'samtools'

    if which(tools) is None:

        logging.error(tools + ' was not found as an executable command. Install ' + tools + ' and re-run the script')
        sys.exit(1)



    try:

        subprocess.check_call(['samtools','quickcheck',os.path.abspath(args.bamfile)],stderr=open(os.devnull, 'wb'))

    except:

        print('.bam file does not exist, is not readable or is not a valid .bam file')
        sys.exit(1)


    bam = pysam.AlignmentFile(os.path.abspath(args.bamfile), "rb")

    if args.type == 'watson':

        watsonbam = pysam.AlignmentFile(os.path.abspath(args.output + '/watson.bam'), "wb", template=bam)

        for read1, read2 in watson_orientation(bam):
            
            watsonbam.write(read1)
            watsonbam.write(read2)

    else:

        crickbam = pysam.AlignmentFile(os.path.abspath(args.output + '/crick.bam'), "wb", template=bam)

        for read1, read2 in crick_orientation(bam):
            
            crickbam.write(read1)
            crickbam.write(read2)


def watson_orientation(bamfilein):

    bam=pysam.AlignmentFile(bamfilein, "rb")
    read_dict = defaultdict(lambda: [None, None])

    for read in bam.fetch():

        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:

            continue

        elif read.is_read1 and read.is_reverse: #if read1 is reverse skip

            continue

        elif not read.is_read1 and not read.is_reverse: #if read2 is not reverse skip

            continue

        else:

            qname = read.query_name

            if qname not in read_dict:

                if read.is_read1:

                    read_dict[qname][0] = read

                else:

                    read_dict[qname][1] = read
            else:

                if read.is_read1:

                    yield read, read_dict[qname][1]
                
                else:

                    yield read_dict[qname][0], read

                del read_dict[qname]




def crick_orientation(bamfilein):

    bam=pysam.AlignmentFile(bamfilein, "rb")
    read_dict = defaultdict(lambda: [None, None])

    for read in bam.fetch():

        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:

            continue

        elif read.is_read1 and not read.is_reverse: #if read1 is not reverse skip

            continue

        elif not read.is_read1 and  read.is_reverse: #if read2 is reverse skip

            continue

        else:

            qname = read.query_name

            if qname not in read_dict:

                if read.is_read1:

                    read_dict[qname][0] = read

                else:

                    read_dict[qname][1] = read
            else:

                if read.is_read1:

                    yield read, read_dict[qname][1]
                
                else:

                    yield read_dict[qname][0], read

                del read_dict[qname]


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
