#!/usr/bin/env python
# Standard library imports
from importlib import metadata
from csv import DictReader
import multiprocessing
import difflib
import logging
import shutil
import json
import os

# Third party imports
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
import pandas as pd
import xlsxwriter

# Local imports
from olctools.accessoryFunctions.accessoryFunctions import \
    GenObject, \
    make_path, \
    MetadataObject, \
    SetupLogging
from genemethods.assemblypipeline.primer_finder_ipcress import \
    CustomIP, \
    make_blastdb
from genemethods.assemblypipeline.legacy_vtyper import Filer


__version__ = metadata.version('in-silico-pcr')
__author__ = 'adamkoziol'


class PrimerValidator:

    def main(self):
        # Inclusion
        self.results, self.inclusivity_miss_dict = self.parse_outputs(
            metadata=self.inclusion_metadata,
            analysistype=self.analysistype,
            group='inclusivity',
            primer_dict=self.primer_dict,
            results_dict=self.results,
            probe_file=self.probe_file,
            amplicon_path=self.inclusion_amplicon_path,
            validator_report_path=self.report_path
        )
        # Exclusion
        self.results, self.exclusivity_miss_dict = self.parse_outputs(
            metadata=self.exclusion_metadata,
            analysistype=self.analysistype,
            group='exclusivity',
            primer_dict=self.primer_dict,
            results_dict=self.results,
            probe_file=self.probe_file,
            amplicon_path=self.exclusion_amplicon_path,
            validator_report_path=self.report_path
        )
        # Only if a probe file was supplied
        if self.probe_file:
            self.probe_dict = self.populate_probe_dict(
                probe_dict=self.probe_dict,
                probe_file=self.probe_file
            )
            self.inclusion_metadata = self.probe_blast(
                metadata=self.inclusion_metadata,
                analysistype='probe',
                group='inclusivity',
                probe_file=self.probe_file,
                primer_dict=self.primer_dict,
                headers=self.fieldnames,
                threads=self.threads
            )
            self.inclusion_metadata = self.parse_blast(
                metadata=self.inclusion_metadata,
                analysistype='probe',
                primer_dict=self.primer_dict,
                fieldnames=self.fieldnames,
                group='inclusivity',
                probe_dict=self.probe_dict,
                cutoff=self.cutoff,
                results_dict=self.results
            )
            self.create_probe_report(
                metadata=self.inclusion_metadata,
                primer_dict=self.primer_dict,
                probe_dict=self.probe_dict,
                group='inclusivity',
                report_path=self.inclusion_report_path,
                validator_report_path=self.report_path,
                miss_dict=self.inclusivity_miss_dict
            )
            self.exclusion_metadata = self.probe_blast(
                metadata=self.exclusion_metadata,
                analysistype='probe',
                group='exclusivity',
                probe_file=self.probe_file,
                primer_dict=self.primer_dict,
                headers=self.fieldnames,
                threads=self.threads,
            )
            self.exclusion_metadata = self.parse_blast(
                metadata=self.exclusion_metadata,
                analysistype='probe',
                primer_dict=self.primer_dict,
                fieldnames=self.fieldnames,
                group='exclusivity',
                probe_dict=self.probe_dict,
                cutoff=self.cutoff,
                results_dict=self.results
            )
            self.create_probe_report(
                metadata=self.exclusion_metadata,
                primer_dict=self.primer_dict,
                probe_dict=self.probe_dict,
                group='exclusivity',
                report_path=self.exclusion_report_path,
                validator_report_path=self.report_path,
                miss_dict=self.exclusivity_miss_dict
            )
        self.create_excel_report(
            report=os.path.join(self.inclusion_report_path, 'custom_epcr_report.csv'),
            group='inclusivity',
            report_path=self.report_path
        )
        self.create_excel_report(
            report=os.path.join(self.exclusion_report_path, 'custom_epcr_report.csv'),
            group='exclusivity',
            report_path=self.report_path
        )
        self.missing_report(
            miss_dict=self.inclusivity_miss_dict,
            group='inclusivity',
            analysis='primer',
            report_path=self.report_path
        )
        self.missing_report(
            miss_dict=self.exclusivity_miss_dict,
            group='exclusivity',
            analysis='primer',
            report_path=self.report_path
        )
        if self.probe_file:
            self.missing_report(
                miss_dict=self.inclusivity_miss_dict,
                group='inclusivity',
                analysis='probe',
                report_path=self.report_path
            )
            self.missing_report(
                miss_dict=self.exclusivity_miss_dict,
                group='exclusivity',
                analysis='probe',
                report_path=self.report_path
            )
        self.results = self.calculate_percentages(results_dict=self.results)
        # Print the populated results dictionary to file using json.dump()
        with open(self.json_report, 'w') as json_report:
            json.dump(self.results, fp=json_report, sort_keys=True, indent=4, separators=(',', ': '))

    @staticmethod
    def parse_outputs(metadata, analysistype, group, primer_dict, results_dict, probe_file, amplicon_path,
                      validator_report_path):
        """
        Extract primer finding results from the metadata objects and create a easily parseable dictionary from
        the outputs
        :param metadata: List of metadata objects for samples in the inclusivity or exclusivity panel
        :param analysistype: String of current analysis type
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param primer_dict: Dictionary of primer name: primer sequence
        :param results_dict: Dictionary of primer name: group: seqid: primer hit details
        :param probe_file: String of the name and path of the provided probe file
        :param amplicon_path: String of the name and path of the output amplicon file
        :param validator_report_path: String of the name and path of the folder into which the validator reports are
        to be written
        :return: Populated results_dict
        """
        logging.info(f'Parsing {group}')
        # Initialise a dictionary to store the names of samples that do not have hits
        empty_results = dict()
        if 'primer' not in empty_results:
            empty_results['primer'] = dict()
        # Extract the names of all the primer sets in the analyses
        for primer in primer_dict:
            if primer not in results_dict:
                results_dict[primer] = dict()
            if primer not in empty_results['primer']:
                empty_results['primer'][primer] = list()
            # Add the inclusivity/exclusivity group to the dictionary
            if group not in results_dict[primer]:
                results_dict[primer][group] = dict()
        # Iterate through all the samples in the supplied panel (inclusivity/exclusivity)
        for sample in metadata:
            # Initialise the .best_hits attribute as required
            if not GenObject.isattr(sample, 'probes'):
                setattr(sample, 'probe', GenObject())
            for primer in primer_dict:
                # Iterate through all the primer sets in the analyses
                for experiment in sample[analysistype].results.datastore:
                    # Remove the trailing _$PRIMER_NUMBER from the experiment to match the supplied primer name
                    # e.g. LinA_0 becomes LinA
                    reduced_experiment = '_'.join(experiment.split('_')[:-1])
                    # Ensure that the primer matches the modified primer name
                    if reduced_experiment != primer:
                        continue
                    # Create the probe_file attribute using the name and path of the file
                    sample.probe.probe_file = probe_file
                    # Extract the primer hit details from the metadata using the unmodified primer name
                    for contig, contig_results in sample[analysistype].results[experiment].datastore.items():
                        # The amplicon range attribute is unnecessary for subsequent analyses, and takes up a lot
                        # of room in the output dictionary (it is a list of every position of the hit sequence)
                        contig_results.amplicon_range = list()
                        # Use the built-in dump method of the GenObject to get the primer details into results_dict
                        results_dict[primer][group][sample.name] = contig_results.dump()
                        if probe_file:
                            PrimerValidator.write_amplicon_sequence(
                                sample=sample,
                                contig_results=contig_results,
                                primer=primer,
                                amplicon_path=amplicon_path,
                                validator_report_path=validator_report_path
                            )
                # If there were no results for the SEQID, add an empty dictionary in place of the primer details
                if sample.name not in results_dict[primer][group]:
                    results_dict[primer][group][sample.name] = dict()
                    # Add the sample name to the list of missing samples in the dictionary
                    if sample.name not in empty_results['primer'][primer]:
                        empty_results['primer'][primer].append(sample.name)
                logging.debug(f'Sample {sample.name} in panel {group}, primer set {primer} has the following outputs: '
                              f'{results_dict[primer][group][sample.name]}')
        return results_dict, empty_results

    @staticmethod
    def write_amplicon_sequence(sample, contig_results, primer, amplicon_path, validator_report_path):
        """
        Extract the amplicon sequence from the GenObject, and write it to file with SeqIO
        :param sample: Metadata object of current sample
        :param contig_results: GenObject storing ipcress results of current sample/primer combination
        :param primer: String of the name of the current primer
        :param amplicon_path: String of the name and path of the output amplicon file
        :param validator_report_path: String of the name and path of the folder into which the validator reports are
        to be written
        """
        # Initialise the header string as the sample name_primer name
        header = f'{sample.name}_{primer}'
        # Set the name of the file in which the amplicons will be saved
        amplicon_file_name = f'{header}_amplicon.fasta'
        amplicon_file = os.path.join(
            amplicon_path,
            amplicon_file_name
        )
        # Initialise GenObjects as required
        if not GenObject.isattr(sample.probe, 'amplicons'):
            sample.probe.amplicons = GenObject()
        if not GenObject.isattr(sample.probe, primer):
            sample.probe[primer] = GenObject()
        # Create the amplicon_file attribute
        sample.probe[primer].amplicon_file = amplicon_file
        # Create a sequence record with SeqRecord
        record = SeqRecord(Seq(contig_results.sequence),  # Use Seq to convert the string of the sequence
                           id=header,
                           description='')
        # Write the SeqRecord to the amplicon file
        with open(amplicon_file, 'w') as fasta:
            SeqIO.write(record, fasta, 'fasta')
        # Copy the amplicon file to the validator report path
        shutil.copyfile(amplicon_file, os.path.join(validator_report_path, amplicon_file_name))

    @staticmethod
    def probe_blast(metadata, analysistype, group, primer_dict, probe_file, headers, threads):
        """
        Run BLAST of the probe sequence against extracted amplicons
        :param metadata: List of metadata objects for samples in the inclusivity or exclusivity panel
        :param analysistype: String of current analysis type
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param primer_dict: Dictionary of primer name: primer sequence
        :param probe_file: String of the name and path of the provided probe file
        :param headers: List of all column headers used in the BLAST analyses
        :param threads: Integer of the number of threads for BLAST to use
        :return: updated metadata
        """
        logging.info(f'Running {group} probe BLAST analyses')
        for sample in metadata:
            for primer in primer_dict:
                try:
                    sample[analysistype][primer].blastresults = '{of}_{primer}_probe_blast_results.tsv' \
                        .format(
                        of=os.path.join(sample.general.outputdirectory, sample.name),
                        primer=primer
                    )
                    # Try to remove empty reports, and run the blast analyses
                    try:
                        size = os.path.getsize(sample[analysistype][primer].blastresults)
                        if size == 0:
                            os.remove(sample[analysistype][primer].blastresults)
                    except FileNotFoundError:
                        pass
                    # Check to see if the results attribute is empty
                    if not os.path.isfile(sample[analysistype][primer].blastresults):
                        db = os.path.splitext(probe_file)[0]
                        if not os.path.isfile(f'{db}.nhr'):
                            make_blastdb(formattedprimers=probe_file)
                        # BLAST command line call.
                        blastn = NcbiblastnCommandline(
                            query=sample.probe[primer].amplicon_file,
                            db=db,
                            evalue=1,
                            task='blastn-short',
                            word_size=4,
                            dust='no',
                            penalty=-2,
                            reward=3,
                            gapopen=5,
                            gapextend=5,
                            num_alignments=1000000,
                            num_threads=threads,
                            outfmt='6 qseqid sseqid positive mismatch gaps evalue bitscore slen '
                                   'length qstart qend qseq sstart send sseq',
                            out=sample[analysistype][primer].blastresults
                        )
                        # Save the blast command in the metadata
                        sample[analysistype][primer].blastcommand = str(blastn)
                        if sample.general.bestassemblyfile != 'NA':
                            # Run the blastn command
                            blastn()
                            # Add a header to the report
                            with open(sample[analysistype][primer].blastresults, 'r+') as f:
                                # Read in the information from the blastresults file
                                content = f.read()
                                # Go back to the start of the file
                                f.seek(0, 0)
                                # Write the formatted header (\n) followed by the content to the file
                                f.write('\t'.join(headers) + '\n' + content)
                except AttributeError:
                    pass
        return metadata

    @staticmethod
    def populate_probe_dict(probe_dict, probe_file):
        """
        Populate probe_dict with probe name: string of probe sequence
        :param probe_dict: Dictionary of probe name: string of probe sequence
        :param probe_file: String of the name and path of the provided probe file
        :return: updated probe_dict
        """
        # Open the probe file, and populate the dictionary with the record id: string of the record sequence
        with open(probe_file, 'r') as probes:
            for record in SeqIO.parse(probes, 'fasta'):
                probe_dict[record.id] = str(record.seq)
        return probe_dict

    @staticmethod
    def parse_blast(metadata, analysistype, primer_dict, fieldnames, group, probe_dict, cutoff, results_dict):
        """

        :param metadata: List of metadata objects for samples in the inclusivity or exclusivity panel
        :param analysistype: String of current analysis type
        :param primer_dict: Dictionary of primer name: primer sequence
        :param fieldnames: List of BLAST headers used in the analysis
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param probe_dict: Dictionary of probe name: string of probe sequence
        :param cutoff: Integer of the percent cutoff value to determine whether a probe hits an amplicon
        :param results_dict: Dictionary of primer name: group: seqid: primer hit details
        :return: updated metadata
        """
        logging.info(f'Parsing {group} BLAST outputs')
        for sample in metadata:
            for primer in primer_dict:
                try:
                    # Only process the samples that were analysed with BLAST
                    if os.path.isfile(sample[analysistype][primer].blastresults):
                        # Open blast output csv file
                        tsvfile = open(sample[analysistype][primer].blastresults)
                        # Skip header
                        tsvfile.readline()
                        # Open the sequence profile file as a dictionary
                        blastdict = DictReader(tsvfile, fieldnames=fieldnames, dialect='excel-tab')
                        hit_dict = dict()
                        # Go through each BLAST result
                        for row in blastdict:
                            # Extract the probe name from the subject_id
                            probe = row['subject_id']
                            positives = int(row['positives'])
                            # sequence_name_primer = row['query_id']
                            # Initialise the primer name in the dictionary as required
                            if primer not in hit_dict:
                                hit_dict[primer] = dict()
                            if probe not in hit_dict[primer]:
                                hit_dict[primer][probe] = dict()
                            # Populate the hit_dict
                            # Find hits that exceed the cutoff
                            percent_id = (positives / len(probe_dict[probe]) * 100)
                            if percent_id > cutoff:
                                row['percent_id'] = percent_id
                                # If the dictionary is still empty, add the number of positives: blast details
                                if not hit_dict[primer][probe]:
                                    hit_dict[primer][probe][percent_id] = row
                                # Extract the number of positives from the previous best result in the dictionary
                                previous_best = list(hit_dict[primer][probe].keys())[0]
                                # Check if the current number of positives is better than the previous best
                                if percent_id > previous_best:
                                    # Add the current positives: blast details to the dictionary
                                    hit_dict[primer][probe][percent_id] = row
                                    # Delete the previous best positives key from the dictionary
                                    del (hit_dict[primer][probe][previous_best])
                        # Iterate through hit_dict
                        for probe, detail_dict in hit_dict[primer].items():
                            # Initialise the probe attribute as required
                            if not GenObject.isattr(sample[analysistype][primer], probe):
                                sample[analysistype][primer][probe] = GenObject()
                            # Iterate through the percent identity and the details in detail_dict
                            for percent_id, details in detail_dict.items():
                                # Initialise the details attribute as required
                                if not GenObject.isattr(sample[analysistype][primer][probe], 'details'):
                                    sample[analysistype][primer][probe].details = details
                                # Add the mismatch_details key to the details dictionary. The mismatch details are
                                # a string of all the mismatches e.g. 32A>T;43A>C
                                details['mismatch_details'] = PrimerValidator.return_mismatches(
                                    sample=sample,
                                    primer=primer,
                                    probe=probe
                                )
                                # Create the details attribute as the details dictionary
                                sample[analysistype][primer][probe].details = details
                                # Populate the details dictionary with the probe information
                                if sample.name not in results_dict[primer][group]:
                                    results_dict[primer][group][sample.name] = dict()
                                if 'probe' not in results_dict[primer][group][sample.name]:
                                    results_dict[primer][group][sample.name]['probe'] = dict()
                                if probe not in results_dict[primer][group][sample.name]['probe']:
                                    results_dict[primer][group][sample.name]['probe'][probe] = dict()
                                # Update results_dict with the details dictionary
                                results_dict[primer][group][sample.name]['probe'][probe] = details
                except AttributeError:
                    pass
        return metadata

    @staticmethod
    def create_probe_report(metadata, primer_dict, probe_dict, group, report_path, validator_report_path, miss_dict):
        """
        Create an Excel report of the probe outputs
        :param metadata: List of metadata objects for samples in the inclusivity or exclusivity panel
        :param primer_dict: Dictionary of primer name: primer sequence
        :param probe_dict: Dictionary of probe name: string of probe sequence
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param report_path: String of the name and path of the folder into which the group-specific reports are to be
        written
        :param validator_report_path: String of the name and path of the folder into which the validator reports are
        to be written
        :param miss_dict: Dictionary of 'primer/probe':primer/probe name: list of samples without hits
        """
        logging.info('Creating probe report')
        # Prepare the Excel file, and set necessary variables
        report_name, report, workbook, worksheet, bold, courier, row, col = PrimerValidator.prep_worksheet(
            group=group,
            report_path=report_path,
            analysis='probe'
        )
        # Initialise a list of all the headers with 'Strain'
        headers = [
            'Sample',
            'Primer',
            'Probe',
            'PercentIdentity',
            'NumMismatches',
            'MismatchDetails',
            'QuerySequence',
            'ProbeSequence'
        ]
        # Create a dictionary to store the desired column widths
        columnwidth = {
            0: 15,  # Sample
            1: 10,  # Primer
            2: 10,  # Probe
            3: 17,  # PercentIdentity
            4: 15,  # NumMismatches
            5: 17,  # MismatchDetails
            6: 60,  # QuerySequence
            7: 60,  # ProbeSequence
        }
        # Set width of each column appropriately
        for header in headers:
            worksheet.write(row, col, header, bold)
            worksheet.set_column(col, col, columnwidth[col])
            col += 1
        missing = list()
        for sample in metadata:
            # Initialise a list to store all the data for each strain
            for primer in primer_dict:
                for probe in probe_dict:
                    # Increment the row and reset the column to zero in preparation of writing results
                    col = 0
                    try:
                        mismatch_details = str()
                        # If there are mismatches, format a string with the details e.g. 32A>T;43A>C
                        if sample.probe[primer][probe].details['percent_id'] < 100:
                            mismatch_details = PrimerValidator.return_mismatches(
                                sample=sample,
                                primer=primer,
                                probe=probe
                            )
                        # Add the appropriate information to the list
                        data = [
                            sample.name,
                            primer,
                            probe,
                            sample.probe[primer][probe].details['percent_id'],
                            sample.probe[primer][probe].details['mismatches'],
                            mismatch_details,
                            sample.probe[primer][probe].details['query_sequence'],
                            sample.probe[primer][probe].details['subject_sequence']
                        ]
                    # If there were no hits, the necessary attributes will not be present
                    except AttributeError:
                        # Initialise the necessary keys as required in the dictionary that stores samples that
                        # do not have hits for the current primer/probe combination
                        if 'probe' not in miss_dict:
                            miss_dict['probe'] = dict()
                        if probe not in miss_dict['probe']:
                            miss_dict['probe'][probe] = list()
                        # Add the sample name to the list in the dictionary
                        if sample.name not in miss_dict['probe'][probe]:
                            miss_dict['probe'][probe].append(sample.name)
                        # Only samples that had no hits at all will not have the probe_file attribute
                        if not GenObject.isattr(sample.probe, 'probe_file'):
                            # Add the sample name only once
                            if sample.name not in missing:
                                data = [sample.name]
                                missing.append(sample.name)
                            # If the sample name is already present, add an empty list
                            else:
                                data = list()
                        # Add an empty list for the primer/probe combinations that do not have hits (at least
                        # one other combination had a hit)
                        else:
                            data = list()
                    # Write out the data to the spreadsheet
                    if data:
                        row += 1
                        for results in data:
                            worksheet.write(row, col, results, courier)
                            col += 1
        # Close the workbook
        workbook.close()
        # Copy the report to the combined report path
        shutil.copyfile(report, os.path.join(validator_report_path, report_name))
        return miss_dict

    @staticmethod
    def prep_worksheet(group, report_path, analysis):
        """
        Create a workbook with a worksheet. Also set font styles for use in writing the report
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param report_path: String of the name and path of the folder into which the group-specific reports are to be
        written
        :param analysis: String of the current analysis type: primer or probe
        :return: report_name: String of the name of the report
        :return: report: String of the name and path of the report
        :return: workbook: xlsxwriter.Workbook object
        :return: worksheet: xlsxwriter.Workbook.add_worksheet object
        :return: bold: xlsxwriter.Workbook.add_format object setting font to bold, Courier New 10
        :return: courier: xlsxwriter.Workbook.add_format object setting font to regular, Courier New 10
        :return: row: Integer for the initial row to use (set to 0)
        :return: col: Integer for the intial column to use (set to 0)
        """
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        report_name = f'{group}_{analysis}_report.xlsx'
        report = os.path.join(
            report_path, report_name
        )
        workbook = xlsxwriter.Workbook(report)
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 10})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 10})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        return report_name, report, workbook, worksheet, bold, courier, row, col

    @staticmethod
    def return_mismatches(sample, primer, probe):
        """
        Create a string of the mismatch details for each primer/probe result
        :param sample: Metadata object of current sample
        :param primer: String of the name of the current primer
        :param probe: String of the name of the current probe
        :return: String of the mismatches for each probe
        """
        # Initialise a string to store the mismatch details
        mismatches = str()
        # Useful: https://towardsdatascience.com/find-the-difference-in-python-68bbd000e513
        seq_matcher = difflib.SequenceMatcher(
            None,
            sample.probe[primer][probe].details['query_sequence'],
            sample.probe[primer][probe].details['subject_sequence'])
        for tag, query_start, query_stop, subject_start, subject_stop in seq_matcher.get_opcodes():
            # Using the SequenceMatcher functionality, the only tag we're interested in is 'replace'.
            if tag == 'replace':
                # If there was a previous mismatch, append a semicolon to the string to delimit the mismatches
                if mismatches:
                    mismatches += ';'
                # Add the query stop position, the subject base at that position as well as the query
                # base at that position to the string e.g. 32A>T
                mismatches += \
                    f'{query_stop}{sample.probe[primer][probe].details["subject_sequence"][subject_start]}>' \
                    f'{sample.probe[primer][probe].details["query_sequence"][query_start]}'
        return mismatches

    @staticmethod
    def create_excel_report(report, group, report_path):
        """
        Use pandas to read in a supplied CSV report, and convert it to a properly-formatted Excel file
        :param report: String of the name and path of the CSV report to process
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param report_path: String of the name and path of the folder into which the group-specific reports are to be
        written
        """
        # Read in the CSV report
        df = pd.read_csv(report)
        # https://stackoverflow.com/a/61617835
        # Create an pandas.ExcelWriter object
        writer = pd.ExcelWriter(os.path.join(report_path, f'{group}_primer_report.xlsx'))
        # Convert the dataframe imported from the CSV report to Excel format. Set the name of the first sheet to
        # 'Sheet1', do not index (ignore the primary key column), and use 'ND' for NaN cells
        df.to_excel(writer, sheet_name='Sheet1', index=False, na_rep='ND')
        for column in df:
            # Find the length of the longest cell contents
            column_width = max(df[column].astype(str).map(len).max(), len(column))
            # If the width is less than 25, set it to 25 in order for the headers to fit
            column_width = column_width if column_width >= 25 else 25
            # Find the index of the column
            col_idx = df.columns.get_loc(column)
            # Set the width of the column
            writer.sheets['Sheet1'].set_column(col_idx, col_idx, column_width)
        # Save the changes to the writer object
        writer.save()

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
        # Add the version of the pipeline to the outputs
        results_dict['validator_version'] = __version__
        return results_dict

    @staticmethod
    def missing_report(miss_dict, group, analysis, report_path):
        """
        Create an Excel report containing the sample names that do not have results for the current group
        (inclusivity/exclusivity) and analysis type (primer/probe)
        :param miss_dict: Dictionary of 'primer/probe':primer/probe name: list of samples without hits
        :param group: String of whether the group is in the inclusivity or exclusivity panel
        :param analysis: String of the current analysis type: primer or probe
        :param report_path: String of the name and path of the folder into which the group-specific reports are to be
        written
        """
        logging.info(f'Creating report for samples with no hits for {group} {analysis}')
        # Create the workbook
        report_name, report, workbook, worksheet, bold, courier, row, col = PrimerValidator.prep_worksheet(
            group=group,
            report_path=report_path,
            analysis=f'{analysis}_missing'
        )
        # Initialise a list to store the primers/probes in the dictionary
        headers = list()
        for primer, sample_list in miss_dict[analysis].items():
            headers.append(primer)
        # Sort the headers
        headers = sorted(headers)
        # Write the primer/probe names to the header
        for header in headers:
            worksheet.write(row, col, header, bold)
            # Set the width of the column to 15 (fits OLC SEQIDs well e.g. 2014-SEQ-0276)
            worksheet.set_column(col, col, 15)
            col += 1
        # Reset the column to 0 in order to add the data
        col = 0
        # Iterate through all the primers/probes
        for header in headers:
            # Reset the row number to 1 for each list of sample names
            row = 1
            for sample_name in miss_dict[analysis][header]:
                # Write the sample name to the worksheet
                worksheet.write(row, col, sample_name, courier)
                row += 1
            # Increase the column number for each primer/probe
            col += 1
        # Close the workbook
        workbook.close()

    def __init__(self, inclusion_metadata, exclusion_metadata, primers, report_path, inclusion_report_path,
                 exclusion_report_path, probe_file, cutoff):
        self.inclusion_metadata = inclusion_metadata
        self.exclusion_metadata = exclusion_metadata
        self.primer_dict = primers
        if report_path.startswith('~'):
            self.report_path = os.path.abspath(os.path.expanduser(os.path.join(report_path)))
        else:
            self.report_path = os.path.abspath(os.path.join(report_path))
        make_path(self.report_path)
        self.inclusion_report_path = inclusion_report_path
        self.exclusion_report_path = exclusion_report_path
        self.inclusion_amplicon_path = os.path.join(
            os.path.dirname(self.inclusion_report_path),
            'amplicons'
        )
        make_path(self.inclusion_amplicon_path)
        self.exclusion_amplicon_path = os.path.join(
            os.path.dirname(self.exclusion_report_path),
            'amplicons'
        )
        make_path(self.exclusion_amplicon_path)
        if probe_file:
            if probe_file.startswith('~'):
                self.probe_file = os.path.abspath(os.path.expanduser(os.path.join(probe_file)))
            else:
                self.probe_file = os.path.abspath(os.path.join(probe_file))
        else:
            self.probe_file = str()
        self.cutoff = cutoff
        try:
            assert 49 < self.cutoff < 101
        except AssertionError:
            logging.error(f'The supplied cutoff value of {self.cutoff} is not in the acceptable range of 50 - 100')
            raise SystemExit
        if not os.path.isfile(self.probe_file) and probe_file is not None:
            logging.warning(
                f'Could not located supplied probe file: {self.probe_file}. Please ensure that you entered the name '
                f'and path of the file correctly')
            raise SystemExit
        self.json_report = os.path.join(self.report_path, 'inclusivity_exclusivity_report.json')
        self.analysistype = 'custom_epcr'
        self.results = dict()
        self.probe_dict = dict()
        self.inclusivity_miss_dict = dict()
        self.exclusivity_miss_dict = dict()
        self.threads = multiprocessing.cpu_count() - 1
        self.headers = [
            'Sample,'
            'Gene',
            'Contig',
            'GenomeLocation',
            'AmpliconSize',
            'Orientation',
            'ForwardMismatches',
            'ForwardMismatchDetails',
            'ForwardLength',
            'ReverseMismatches',
            'ReverseMismatchDetails',
            'ReverseLength',
            'ForwardPrimer',
            'ForwardQuery',
            'ReversePrimer',
            'ReverseQuery'
        ]
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = [
            'query_id',
            'subject_id',
            'positives',
            'mismatches',
            'gaps',
            'evalue',
            'bit_score',
            'subject_length',
            'alignment_length',
            'query_start',
            'query_end',
            'query_sequence',
            'subject_start',
            'subject_end',
            'subject_sequence'
        ]


def cli():
    parser = ArgumentParser(description='Perform in silico PCR analyses to validate primer sets')
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'%(prog)s commit {__version__}'
    )
    parser.add_argument(
        '-i', '--inclusion_sequencepath',
        required=True,
        help='Path of folder containing inclusion .fasta files to process.'
    )
    parser.add_argument(
        '-e', '--exclusion_sequencepath',
        required=True,
        help='Path of folder containing exclusion .fasta files to process.'
    )
    parser.add_argument(
        '-r', '--report_path',
        required=True,
        help='Name and path of folder in which inclusivity/exclusivity reports are to be created'
    )
    parser.add_argument(
        '-m', '--mismatches',
        default=2,
        type=int,
        choices=[0, 1, 2, 3],
        help='Number of mismatches allowed [0-3]. Default is 1'
    )
    parser.add_argument(
        '-pf', '--primerfile',
        help='Absolute path and name of the primer file to test'
    )
    parser.add_argument(
        '-min', '--minampliconsize',
        default=0,
        type=int,
        help='Minimum size of amplicons. Default is 0'
    )
    parser.add_argument(
        '-max', '--maxampliconsize',
        default=1500,
        type=int,
        help='Maximum size of amplicons. Default is 1500'
    )
    parser.add_argument(
        '-cb', '--contigbreaks',
        action='store_true',
        help='If ipcress cannot find an amplicon, search the genome for the primers '
             'and return a positive result if they are found on separate contigs '
    )
    parser.add_argument(
        '-rb', '--range_buffer',
        type=int,
        default=0,
        help='Increase the buffer size between amplicons in a sequence. Useful if you have '
             'overlapping primer sets (e.g. vtyper), and you are looking for the best one. '
             'Default is 0'
    )
    parser.add_argument(
        '-p', '--probe_file',
        type=str,
        default=None,
        help='Name and path of file containing probe sequences to test against extracted amplicons'
    )
    parser.add_argument(
        '-c', '--cutoff',
        type=int,
        default=80,
        help='Percent cutoff value to determine whether a probe hits an amplicon. Minimum value is 50. Maximum is 100'
             'default is 80'
    )
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Allow debug-level logging to be printed to the terminal'
    )
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
        contigbreaks=arguments.contigbreaks
    )
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
        contigbreaks=arguments.contigbreaks
    )
    exclusion.main()
    validator = PrimerValidator(
        inclusion_metadata=inclusion.metadata,
        exclusion_metadata=exclusion.metadata,
        primers=inclusion.forward_dict,
        report_path=arguments.report_path,
        inclusion_report_path=inclusion.reportpath,
        exclusion_report_path=exclusion.reportpath,
        probe_file=arguments.probe_file,
        cutoff=arguments.cutoff
    )
    validator.main()


if __name__ == '__main__':
    cli()
    logging.info('Primer validation complete')
