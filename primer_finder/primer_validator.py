#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import \
    make_path, \
    MetadataObject, \
    SetupLogging
from genemethods.assemblypipeline.primer_finder_ipcress import CustomIP
from genemethods.assemblypipeline.legacy_vtyper import Filer
from argparse import ArgumentParser
import logging
import json
import os

__author__ = 'adamkoziol'


class PrimerValidator:

    def main(self):
        # Inclusion
        self.results = self.parse_outputs(
            metadata=self.inclusion_metadata,
            analysistype=self.analysistype,
            group='inclusivity',
            primer_dict=self.primer_dict,
            results_dict=self.results
        )
        # Exclusion
        self.results = self.parse_outputs(
            metadata=self.exclusion_metadata,
            analysistype=self.analysistype,
            group='exclusivity',
            primer_dict=self.primer_dict,
            results_dict=self.results
        )
        self.results = self.calculate_percentages(results_dict=self.results)
        # Print the populated results dictionary to file using json.dump()
        with open(self.json_report, 'w') as json_report:
            json.dump(self.results, fp=json_report, sort_keys=True, indent=4, separators=(',', ': '))

    @staticmethod
    def parse_outputs(metadata, analysistype, group, primer_dict, results_dict):
        """
        Extract primer finding results from the metadata objects and create a easily parseable dictionary from
        the outputs
        :param metadata: List of metadata objects for samples in the inclusivity or exclusivity panel
        :param analysistype: String of current analysis type
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param primer_dict: Dictionary of primer name: primer sequence
        :param results_dict: Dictionary of primer name: group: seqid: primer hit details
        :return: Populated results_dict
        """
        logging.info(f'Parsing {group}')
        # Extract the names of all the primer sets in the analyses
        for primer in primer_dict:
            if primer not in results_dict:
                results_dict[primer] = dict()
            # Add the inclusivity/exclusivity group to the dictionary
            if group not in results_dict[primer]:
                results_dict[primer][group] = dict()
        # Iterate through all the samples in the supplied panel (inclusivity/exclusivity)
        for sample in metadata:
            for primer in primer_dict:
                # Iterate through all the primer sets in the analyses
                for experiment in sample[analysistype].results.datastore:
                    # Remove the trailing _$PRIMER_NUMBER from the experiment to match the supplied primer name
                    # e.g. LinA_0 becomes LinA
                    reduced_experiment = '_'.join(experiment.split('_')[:-1])
                    # Ensure that the primer matches the modified primer name
                    if reduced_experiment != primer:
                        continue
                    # Extract the primer hit details from the metadata using the unmodified primer name
                    for contig, contig_results in sample[analysistype].results[experiment].datastore.items():
                        # The amplicon range attribute is unnecessary for subsequent analyses, and takes up a lot
                        # of room in the output dictionary (it is a list of every position of the hit sequence)
                        contig_results.amplicon_range = list()
                        # Use the built-in dump method of the GenObject to get the primer details into results_dict
                        results_dict[primer][group][sample.name] = contig_results.dump()
                # If there were no results for the SEQID, add an empty dictionary in place of the primer details
                if sample.name not in results_dict[primer][group]:
                    results_dict[primer][group][sample.name] = dict()
                logging.debug(f'Sample {sample.name} in panel {group}, primer set {primer} has the following outputs: '
                              f'{results_dict[primer][group][sample.name]}')
        return results_dict

    @staticmethod
    def calculate_percentages(results_dict):
        """
        For each primer set in the analyses, calculate the number of samples with amplicons created by the primer set,
        the total number of samples in the panel, and the percentage of positive/total
        :param results_dict: Dictionary of primer name: group: seqid: primer hit details
        :return: Updated results_dict
        """
        for primer, group_dict in results_dict.items():
            for group, sample_dict in group_dict.items():
                # Initialise variables to track the number of samples with primer details (the primers match), and the
                # total number of samples
                positive = 0
                total = 0
                # Extract the primer SEQID, and primer details dictionary
                for sample, contig_results in sample_dict.items():
                    # Increment the total for every sample
                    total += 1
                    # Only increment the positive if there are results
                    if contig_results:
                        positive += 1
                # Update the dictionary with the number of positive and total number of samples
                results_dict[primer][group]['positive'] = positive
                results_dict[primer][group]['total'] = total
                # Determine the percentage of positives
                try:
                    results_dict[primer][group]['percent'] = float(positive) / float(total) * 100
                except ZeroDivisionError:
                    results_dict[primer][group]['percent'] = 0.0
                logging.debug(
                    f'Primer set {primer} in group {group} had {results_dict[primer][group]["positive"]}, out of '
                    f'{results_dict[primer][group]["total"]}, for a percentage of '
                    f'{results_dict[primer][group]["percent"]}')
        return results_dict

    def __init__(self, inclusion_metadata, exclusion_metadata, primers, report_path):
        self.inclusion_metadata = inclusion_metadata
        self.exclusion_metadata = exclusion_metadata
        self.primer_dict = primers
        if report_path.startswith('~'):
            self.report_path = os.path.abspath(os.path.expanduser(os.path.join(report_path)))
        else:
            self.report_path = os.path.abspath(os.path.join(report_path))
        make_path(self.report_path)
        self.json_report = os.path.join(self.report_path, 'inclusivity_exclusivity_report.json')
        self.analysistype = 'custom_epcr'
        self.results = dict()


def cli():
    parser = ArgumentParser(description='Perform in silico PCR analyses to validate primer sets')
    parser.add_argument(
        '-i', '--inclusion_sequencepath',
        required=True,
        help='Path of folder containing inclusion .fasta files to process.')
    parser.add_argument(
        '-e', '--exclusion_sequencepath',
        required=True,
        help='Path of folder containing exclusion .fasta files to process.')
    parser.add_argument(
        '-r', '--report_path',
        required=True,
        help='Name and path of folder in which inclusivity/exclusivity reports are to be created')
    parser.add_argument(
        '-m', '--mismatches',
        default=2,
        type=int,
        choices=[0, 1, 2, 3],
        help='Number of mismatches allowed [0-3]. Default is 1')
    parser.add_argument(
        '-pf', '--primerfile',
        help='Absolute path and name of the primer file to test')
    parser.add_argument(
        '-min', '--minampliconsize',
        default=0,
        type=int,
        help='Minimum size of amplicons. Default is 0')
    parser.add_argument(
        '-max', '--maxampliconsize',
        default=1500,
        type=int,
        help='Maximum size of amplicons. Default is 1500')
    parser.add_argument(
        '-cb', '--contigbreaks',
        action='store_true',
        help='If ipcress cannot find an amplicon, search the genome for the primers '
             'and return a positive result if they are found on separate contigs ')
    parser.add_argument(
        '-rb', '--range_buffer',
        type=int,
        default=0,
        help='Increase the buffer size between amplicons in a sequence. Useful if you have '
             'overlapping primer sets (e.g. vtyper), and you are looking for the best one. '
             'Default is 0')
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Allow debug-level logging to be printed to the terminal')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=arguments.debug)
    # Create attributes to work with the CustomIP class
    arguments.primer_format = 'fasta'
    arguments.export_amplicons = False
    # Inclusion
    arguments.sequencepath = arguments.inclusion_sequencepath
    arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    arguments.runmetadata.samples = Filer.filer(arguments)
    logging.info('Processing inclusivity dataset')
    inclusion = CustomIP(
        metadataobject=arguments.runmetadata.samples,
        sequencepath=arguments.sequencepath,
        reportpath=os.path.join(arguments.sequencepath, 'reports'),
        primerfile=arguments.primerfile,
        min_amplicon_size=arguments.minampliconsize,
        max_amplicon_size=arguments.maxampliconsize,
        primer_format=arguments.primer_format,
        mismatches=arguments.mismatches,
        export_amplicons=arguments.export_amplicons,
        contigbreaks=arguments.contigbreaks)
    inclusion.main()

    # Exclusion
    arguments.sequencepath = arguments.exclusion_sequencepath
    arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    arguments.runmetadata.samples = Filer.filer(arguments)
    logging.info('Processing exclusivity dataset')
    exclusion = CustomIP(
        metadataobject=arguments.runmetadata.samples,
        sequencepath=arguments.sequencepath,
        reportpath=os.path.join(arguments.sequencepath, 'reports'),
        primerfile=arguments.primerfile,
        min_amplicon_size=arguments.minampliconsize,
        max_amplicon_size=arguments.maxampliconsize,
        primer_format=arguments.primer_format,
        mismatches=arguments.mismatches,
        export_amplicons=arguments.export_amplicons,
        contigbreaks=arguments.contigbreaks)
    exclusion.main()
    validator = PrimerValidator(inclusion_metadata=inclusion.metadata,
                                exclusion_metadata=exclusion.metadata,
                                primers=inclusion.forward_dict,
                                report_path=arguments.report_path)
    validator.main()


if __name__ == '__main__':
    cli()
    logging.info('Primer validation complete')
