#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import MetadataObject, SetupLogging
from spadespipeline.primer_finder_bbduk import PrimerFinder
from spadespipeline.legacy_vtyper import Custom, Filer, Vtyper
from argparse import ArgumentParser, RawTextHelpFormatter
from time import time
import logging
import os
__author__ = 'adamkoziol'


def legacy(args):
    # Prep the args object to be used in the legacy script
    SetupLogging(debug=args.debug)
    args.starttime = time()
    args.reportpath = os.path.join(args.sequencepath, 'reports')
    args.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    args.runmetadata.samples = Filer.filer(args)
    if args.analysistype == 'vtyper':
        # Perform vtx typing
        vtyper = Vtyper(inputobject=args,
                        analysistype='vtyper_legacy',
                        mismatches=args.mismatches)
        vtyper.vtyper()
    else:
        epcr = Custom(inputobject=args,
                      analysistype='custom_epcr',
                      primerfile=args.primerfile,
                      ampliconsize=args.maxampliconsize,
                      mismatches=args.mismatches,
                      export_amplicons=args.export_amplicons)
        epcr.main()


def supremacy(args):
    SetupLogging(debug=args.debug)
    # Create supremacy object
    finder = PrimerFinder(path=args.path,
                          sequence_path=args.sequencepath,
                          primer_file=args.primerfile,
                          mismatches=args.mismatches,
                          kmer_length=args.kmerlength,
                          cpus=args.cpus,
                          analysistype='ePCR')
    # Run the script
    finder.main()
    logging.info('ePCR analyses complete')


def cli():
    parser = ArgumentParser(description='Perform in silico PCR analyses')
    subparsers = parser.add_subparsers(title='Available analyses')
    # Create a parental parser from which the legacy and supremacy subparsers can inherit arguments
    parent_parser = ArgumentParser(add_help=False)
    parent_parser.add_argument('-s', '--sequencepath',
                               required=True,
                               help='Path of folder containing .fasta/.fastq(.gz) files to process.')
    parent_parser.add_argument('-m', '--mismatches',
                               default=1,
                               type=int,
                               choices=[0, 1, 2, 3],
                               help='Number of mismatches allowed [0-3]. Default is 1')
    parent_parser.add_argument('-d', '--debug',
                               action='store_true',
                               help='Allow debug-level logging to be printed to the terminal')
    # Create a subparser to run primer finder legacy
    legacy_subparser = subparsers.add_parser(parents=[parent_parser],
                                             name='legacy',
                                             description='Perform in silico PCR on FASTA-formatted files using the '
                                                         'retired ePCR tool from NCBI',
                                             formatter_class=RawTextHelpFormatter,
                                             help='Perform in silico PCR on FASTA-formatted files using the retired '
                                                  'ePCR tool from NCBI')
    legacy_subparser.add_argument('-a', '--analysistype',
                                  default='custom',
                                  choices=['vtyper', 'custom'],
                                  help='Either perform the standard vtyper analysis using the included primer files, \n'
                                       'or supply your own FASTA-formatted (multi-)primer file')
    legacy_subparser.add_argument('-pf', '--primerfile',
                                  help='Absolute path and name of the primer file (in FASTA format) to test. The file\n'
                                       'must have every primer on a separate line AND -F/-R following the name e.g.\n'
                                       '>primer1-F\n'
                                       'ATCGACTGACAC....\n'
                                       '>primer1-R\n'
                                       'ATCGATCGATCGATG....\n'
                                       '>primer2-F\n'
                                       '.......\n')
    legacy_subparser.add_argument('-e', '--export_amplicons',
                                  action='store_true',
                                  help='Export the sequence of the calculated amplicons. Default is False')
    legacy_subparser.add_argument('-mas', '--maxampliconsize',
                                  default=10000,
                                  type=int,
                                  help='Maximum size of amplicons. Default is 10000')
    legacy_subparser.set_defaults(func=legacy)
    # Create a subparser to run primer finder supremacy
    supremacy_subparser = subparsers.add_parser(parents=[parent_parser],
                                                description='Perform in silico PCR on FASTA- or FASTQ-formatted files '
                                                            'using bbduk and SPAdes',
                                                name='supremacy',
                                                formatter_class=RawTextHelpFormatter,
                                                help='Perform in silico PCR on FASTA- or FASTQ-formatted files '
                                                     'using bbduk and SPAdes')
    supremacy_subparser.add_argument('-p', '--path',
                                     required=True,
                                     help='Specify directory in which reports are to be created')
    supremacy_subparser.add_argument('-pf', '--primerfile',
                                     required=True,
                                     help='Absolute path and name of the primer file (in FASTA format) to test. \n'
                                          'The file must have every primer on a separate line AND -F/-R following \n '
                                          'the name e.g.\n'
                                          '>primer1-F\n'
                                          'ATCGACTGACAC....\n'
                                          '>primer1-R\n'
                                          'ATCGATCGATCGATG....\n'
                                          '.......\n')
    supremacy_subparser.add_argument('-n', '--cpus',
                                     default=0,
                                     help='Number of threads. Default is the number of cores in the system')
    supremacy_subparser.add_argument('-k', '--kmerlength',
                                     default='55,77,99,127',
                                     help='The range of kmers used in SPAdes assembly. Default is 55,77,99,127, but \n'
                                          'you can provide a comma-separated list of kmers e.g. 21,33,55,77,99,127, \n'
                                          'or a single kmer e.g. 33')
    supremacy_subparser.set_defaults(func=supremacy)
    # Get the arguments into an object
    arguments = parser.parse_args()
    # Run the appropriate function for each sub-parser.
    if hasattr(arguments, 'func'):
        arguments.func(arguments)
    # If the 'func' attribute doesn't exist, display the basic help
    else:
        parser.parse_args(['-h'])


if __name__ == '__main__':
    cli()
