#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import MetadataObject, SetupLogging
from genemethods.assemblypipeline.primer_finder_ipcress import CustomIP, VtyperIP
from genemethods.assemblypipeline.primer_finder_ultimatum import Ultimatum
from genemethods.assemblypipeline.primer_finder_bbduk import PrimerFinder
from genemethods.assemblypipeline.legacy_vtyper import Custom, Filer, Vtyper
from argparse import ArgumentParser, RawTextHelpFormatter
import logging
import os

__author__ = 'adamkoziol'


def legacy(args):
    # Prep the args object to be used in the legacy script
    SetupLogging(debug=args.debug)
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
                      primer_format=args.primer_format,
                      export_amplicons=args.export_amplicons)
        epcr.main()


def identity(args):
    SetupLogging(debug=args.debug)
    # Create metadata objects for the samples
    args.runmetadata = MetadataObject()
    args.runmetadata.samples = Filer.filer(args)
    if args.analysistype == 'vtyper':
        epcr = VtyperIP(metadataobject=args.runmetadata.samples,
                        analysistype=args.analysistype,
                        reportpath=os.path.join(args.sequencepath, 'reports'))
        epcr.vtyper()
    else:
        epcr = CustomIP(metadataobject=args.runmetadata.samples,
                        sequencepath=args.sequencepath,
                        reportpath=os.path.join(args.sequencepath, 'reports'),
                        primerfile=args.primerfile,
                        min_amplicon_size=args.minampliconsize,
                        max_amplicon_size=args.maxampliconsize,
                        primer_format=args.primer_format,
                        mismatches=args.mismatches,
                        export_amplicons=args.export_amplicons,
                        contigbreaks=args.contigbreaks)
        epcr.main()


def supremacy(args):
    SetupLogging(debug=args.debug)
    # Create supremacy object
    finder = PrimerFinder(sequence_path=args.sequencepath,
                          primer_file=args.primerfile,
                          mismatches=args.mismatches,
                          kmer_length=args.kmerlength,
                          cpus=args.cpus,
                          analysistype='ePCR')
    # Run the script
    finder.main()


def ultimatum(args):
    SetupLogging(debug=args.debug)
    # Create metadata objects for the samples
    args.runmetadata = MetadataObject()
    args.runmetadata.samples = Filer.filer(args)
    finder = Ultimatum(metadataobject=args.runmetadata.samples,
                       sequencepath=args.sequencepath,
                       reportpath=os.path.join(args.sequencepath, 'reports'),
                       primerfile=args.primerfile,
                       primer_format=args.primer_format,
                       mismatches=args.mismatches,
                       export_amplicons=args.export_amplicons)
    finder.main()


def cli():
    parser = ArgumentParser(description='Perform in silico PCR analyses')
    subparsers = parser.add_subparsers(title='Available analyses')
    # Create a parental parser from which the legacy, identity, and supremacy subparsers can inherit arguments
    parent_parser = ArgumentParser(add_help=False)
    parent_parser.add_argument('-s', '--sequencepath',
                               required=True,
                               help='Path of folder containing .fasta/.fastq(.gz) files to process.')
    parent_parser.add_argument('-m', '--mismatches',
                               default=1,
                               type=int,
                               choices=[0, 1, 2, 3],
                               help='Number of mismatches allowed [0-3]. Default is 1')
    parent_parser.add_argument('-pf', '--primerfile',
                               help='Absolute path and name of the primer file to test')
    parent_parser.add_argument('-f', '--primer_format',
                               choices=['epcr', 'fasta'],
                               default='fasta',
                               type=str,
                               help='Format of the supplied primer file. Choices are "epcr" and "fasta". \nDefault is '
                                    'fasta. epcr format describes one white-space delimited PCR reaction per line. \n'
                                    'The following fields be included: name of primer pair, sequence of forward '
                                    'primer, sequence of reverse primer, \nmin amplicon size, max amplicon size e.g.\n'
                                    'vtx1a CCTTTCCAGGTACAACAGCGGTT GGAAACTCATCAGATGCCATTCTGG 0 1500\n'
                                    'vtx1c CCTTTCCTGGTACAACTGCGGTT CAAGTGTTGTACGAAATCCCCTCTGA 0 1500\n'
                                    '.......\n'
                                    'fasta format must have every primer on a separate line AND -F/-R following the '
                                    'name e.g.\n'
                                    '>vtx1a-F\n'
                                    'CCTTTCCAGGTACAACAGCGGTT\n'
                                    '>vtx1a-R\n'
                                    'GGAAACTCATCAGATGCCATTCTGG\n'
                                    '>vtx1c-F\n'
                                    '.......\n')
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
    legacy_subparser.add_argument('-e', '--export_amplicons',
                                  action='store_true',
                                  help='Export the sequence of the calculated amplicons. Default is False')
    legacy_subparser.add_argument('-min', '--minampliconsize',
                                  default=0,
                                  type=int,
                                  help='Minimum size of amplicons. Default is 0')
    legacy_subparser.add_argument('-max', '--maxampliconsize',
                                  default=1500,
                                  type=int,
                                  help='Maximum size of amplicons. Default is 1500')
    legacy_subparser.set_defaults(func=legacy)
    # Create a subparser to run primer finder identity
    identity_subparser = subparsers.add_parser(parents=[parent_parser],
                                               name='identity',
                                               description='Perform in silico PCR on FASTA-formatted files using '
                                                           'ipcress from the exonerate toolbox',
                                               formatter_class=RawTextHelpFormatter,
                                               help='Perform in silico PCR on FASTA-formatted files using '
                                                    'ipcress from the exonerate toolbox')
    identity_subparser.add_argument('-a', '--analysistype',
                                    default='custom',
                                    choices=['vtyper', 'custom'],
                                    help='Either perform vtyper analysis using the included primer files,\n'
                                         'or supply your own (multi-)primer file')
    identity_subparser.add_argument('-e', '--export_amplicons',
                                    action='store_true',
                                    help='Export the sequence of the calculated amplicons. Default is False')
    identity_subparser.add_argument('-min', '--minampliconsize',
                                    default=0,
                                    type=int,
                                    help='Minimum size of amplicons. Default is 0')
    identity_subparser.add_argument('-max', '--maxampliconsize',
                                    default=1500,
                                    type=int,
                                    help='Maximum size of amplicons. Default is 1500')
    identity_subparser.add_argument('-cb', '--contigbreaks',
                                    action='store_true',
                                    help='If ipcress cannot find an amplicon, search the genome for the primers '
                                         'and return a positive result if they are found on separate contigs ')
    identity_subparser.add_argument('-rb', '--range_buffer',
                                    type=int,
                                    default=0,
                                    help='Increase the buffer size between amplicons in a sequence. Useful if you have '
                                         'overlapping primer sets (e.g. vtyper), and you are looking for the best one. '
                                         'Default is 0 "custom" analyses, and 150 for "vtyper" (vtx1c and vtx1d sets '
                                         'don\'t quite overlap, but the script requires an overlap for finding the '
                                         'best primer set)')
    identity_subparser.set_defaults(func=identity)
    # Create a subparser to run primer finder supremacy
    supremacy_subparser = subparsers.add_parser(parents=[parent_parser],
                                                description='Perform in silico PCR on FASTA- or FASTQ-formatted files '
                                                            'using bbduk and SPAdes',
                                                name='supremacy',
                                                formatter_class=RawTextHelpFormatter,
                                                help='Perform in silico PCR on FASTA- or FASTQ-formatted files '
                                                     'using bbduk and SPAdes')
    supremacy_subparser.add_argument('-n', '--cpus',
                                     default=0,
                                     help='Number of threads. Default is the number of cores in the system')
    supremacy_subparser.add_argument('-k', '--kmerlength',
                                     default='55,77,99,127',
                                     help='The range of kmers used in SPAdes assembly. Default is 55,77,99,127, but \n'
                                          'you can provide a comma-separated list of kmers e.g. 21,33,55,77,99,127, \n'
                                          'or a single kmer e.g. 33')
    supremacy_subparser.set_defaults(func=supremacy)
    ultimatum_subparser = subparsers.add_parser(parents=[parent_parser],
                                                description='Search for primer singletons using BLAST',
                                                name='ultimatum',
                                                formatter_class=RawTextHelpFormatter,
                                                help='Search for primer singletons in a sequence using BLAST')
    ultimatum_subparser.add_argument('-e', '--export_amplicons',
                                     action='store_true',
                                     help='Export the sequence of the calculated amplicons. Default is False')
    ultimatum_subparser.set_defaults(func=ultimatum)
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
    logging.info('ePCR analyses complete')
