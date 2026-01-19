#!/usr/bin/env python
"""
Primer validator
================

Validate primer (and optional probe) performance against inclusivity and
exclusivity panels using Custom ePCR (ipcress) outputs and BLAST.

This script orchestrates:
- Running Custom ePCR on inclusivity and exclusivity datasets
- Parsing ePCR outputs into a structured results dictionary
- Optional probe-vs-amplicon BLAST and parsing the top hits
- Generating Excel reports for primers and probes
- Emitting a JSON summary including per-primer statistics
"""

# Standard library imports
from argparse import ArgumentParser
from csv import DictReader
import difflib
from importlib import metadata as meta_data
import json
import logging
import multiprocessing
import os
import shutil
import tempfile
from typing import Any, Dict, List, Optional, Tuple

# Third party imports
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import xlsxwriter

# Local imports
from genemethods.assemblypipeline.legacy_vtyper import Filer
from genemethods.assemblypipeline.primer_finder_ipcress import (
    CustomIP,
    make_blastdb,
)
from olctools.accessoryFunctions.accessoryFunctions import (
    GenObject,
    MetadataObject,
    SetupLogging,
    make_path,
)


__version__ = meta_data.version('in-silico-pcr')
__author__ = 'adamkoziol'


class PrimerValidator:
    """Validate primer and probe performance.

    This class parses Custom ePCR results for primer performance across
    inclusivity and exclusivity panels. When a probe FASTA is provided, it also
    performs BLAST of probes against extracted amplicons and aggregates those
    results.
    """

    def main(self):
        """Run validation for inclusivity and exclusivity panels.

        Steps:
        1) Parse primer finding outputs for inclusivity and exclusivity.
        2) Optionally BLAST probes vs amplicons and parse best hits.
        3) Create Excel reports for primer and (optional) probe analyses.
        4) Emit a JSON report of the consolidated results.
        """

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
            report=os.path.join(
                self.inclusion_report_path, 'custom_epcr_report.csv'
            ),
            group='inclusivity',
            report_path=self.report_path
        )
        self.create_excel_report(
            report=os.path.join(
                self.exclusion_report_path, 'custom_epcr_report.csv'
            ),
            group='exclusivity',
            report_path=self.report_path
        )
        self.missing_report(
            miss_dict=self.inclusivity_miss_dict,
            group='inclusivity',
            analysis='primer',
            report_path=self.report_path,
            primer_dict=self.primer_dict
        )
        self.missing_report(
            miss_dict=self.exclusivity_miss_dict,
            group='exclusivity',
            analysis='primer',
            report_path=self.report_path,
            primer_dict=self.primer_dict
        )
        if self.probe_file:
            self.missing_report(
                miss_dict=self.inclusivity_miss_dict,
                group='inclusivity',
                analysis='probe',
                report_path=self.report_path,
                primer_dict=self.primer_dict
            )
            self.missing_report(
                miss_dict=self.exclusivity_miss_dict,
                group='exclusivity',
                analysis='probe',
                report_path=self.report_path,
                primer_dict=self.primer_dict
            )
        self.results = self.calculate_percentages(results_dict=self.results)
        # Create a summary Excel report (best-effort from the populated results)
        try:
            self.create_summary_report()
        except Exception:
            logging.exception('Failed to create summary report')

        # Print the populated results dictionary to file using json.dump()
        with open(self.json_report, 'w', encoding='utf-8') as json_report:
            json.dump(
                self.results,
                fp=json_report,
                sort_keys=True,
                indent=4,
                separators=(',', ': ')
            )

    @staticmethod
    def parse_outputs(
        *,
        metadata: List[Any],
        analysistype: str,
        group: str,
        primer_dict: Dict[str, str],
        results_dict: Dict[str, Any],
        probe_file: Optional[str],
        amplicon_path: str,
        validator_report_path: str,
    ) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        """Parse ipcress outputs into a results dictionary.

        Args:
            metadata: Samples for the inclusivity or exclusivity panel.
            analysistype: Name of the analysis (e.g., "custom_epcr").
            group: Either "inclusivity" or "exclusivity".
            primer_dict: Mapping of primer name to sequence.
            results_dict: Accumulator mapping for results.
            probe_file: Optional path to probe FASTA file.
            amplicon_path: Directory path to write amplicon FASTAs.
            validator_report_path: Directory for validator reports.

        Returns:
            A tuple of:
            - results_dict populated with primer results for this group.
            - miss_dict tracking samples without hits per primer.
        """
        logging.info('Parsing %s', group)

        # Initialise a dictionary to store the names of samples that do not
        # have hits
        miss_dict = {}
        if 'primer' not in miss_dict:
            miss_dict['primer'] = {}

        # Extract the names of all the primer sets in the analyses
        for primer in primer_dict:
            if primer not in results_dict:
                results_dict[primer] = {}
            if primer not in miss_dict['primer']:
                miss_dict['primer'][primer] = []

            # Add the inclusivity/exclusivity group to the dictionary
            if group not in results_dict[primer]:
                results_dict[primer][group] = {}

        # Iterate through all the samples in the supplied panel
        # (inclusivity/exclusivity)
        for sample in metadata:
            # Initialise the .best_hits attribute as required
            # isattr used as a class helper from external lib.
            is_attr = GenObject.isattr  # type: ignore
            if not is_attr(sample, 'probes'):
                setattr(sample, 'probe', GenObject())
            for primer in primer_dict:
                # Iterate through all the primer sets in the analyses
                for experiment in sample[analysistype].results.datastore:
                    # Remove the trailing _$PRIMER_NUMBER from the experiment
                    # to match the supplied primer name
                    # e.g. LinA_0 becomes LinA
                    reduced_experiment = '_'.join(experiment.split('_')[:-1])

                    # Ensure that the primer matches the modified primer name
                    if reduced_experiment != primer:
                        continue

                    # Create the probe_file attribute using the name and
                    # path of the file
                    sample.probe.probe_file = probe_file

                    # Extract the primer hit details from the metadata using
                    # the unmodified primer name
                    for _, contig_results in (
                        sample[analysistype]
                        .results[experiment]
                        .datastore.items()
                    ):
                        # The amplicon range attribute is unnecessary for
                        # subsequent analyses, and takes up a lot of room in
                        # the output dictionary (it is a list of every
                        # position of the hit sequence)
                        contig_results.amplicon_range = []
                        # Use the built-in dump method of the GenObject to get
                        # the primer details into results_dict
                        results_dict[primer][group][sample.name] = \
                            contig_results.dump()
                        if probe_file:
                            PrimerValidator.write_amplicon_sequence(
                                sample=sample,
                                contig_results=contig_results,
                                primer=primer,
                                amplicon_path=amplicon_path,
                                validator_report_path=validator_report_path
                            )
                # If there were no results for the SEQID, add an empty
                # dictionary in place of the primer details
                if sample.name not in results_dict[primer][group]:
                    results_dict[primer][group][sample.name] = {}
                    # Add the sample name to the list of missing samples in
                    # the dictionary
                    if sample.name not in miss_dict['primer'][primer]:
                        miss_dict['primer'][primer].append(sample.name)

                logging.debug(
                    'Sample %s in panel %s, primer set %s has the following '
                    'outputs: %s',
                    sample.name, group, primer,
                    results_dict[primer][group][sample.name]
                )
        return results_dict, miss_dict

    @staticmethod
    def write_amplicon_sequence(
        *,
        sample: MetadataObject,
        contig_results: GenObject,
        primer: str,
        amplicon_path: str,
        validator_report_path: str,
    ) -> None:
        """Write the amplicon sequence for a sample/primer to FASTA.

        Args:
            sample: Current sample metadata object.
            contig_results: ipcress results object for sample/primer.
            primer: Primer name.
            amplicon_path: Output directory for amplicon FASTAs.
            validator_report_path: Report directory for copied FASTAs.
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
        # Use Seq to convert the string of the sequence
        record = SeqRecord(
            Seq(contig_results.sequence),
            id=header,
            description='',
        )
        # Write the SeqRecord to the amplicon file
        with open(amplicon_file, 'w', encoding='utf-8') as fasta:
            SeqIO.write(record, fasta, 'fasta')
        # Copy the amplicon file to the validator report path
        shutil.copyfile(
            amplicon_file,
            os.path.join(validator_report_path, amplicon_file_name),
        )

    @staticmethod
    def probe_blast(
        *,
        metadata: List[MetadataObject],
        analysistype: str,
        group: str,
        primer_dict: Dict[str, str],
        probe_file: str,
        headers: List[str],
        threads: int,
    ) -> List[MetadataObject]:
        """Run BLAST of probe sequences against extracted amplicons.

        Args:
            metadata: Samples in the target panel.
            analysistype: Analysis name (e.g., "probe").
            group: Either "inclusivity" or "exclusivity".
            primer_dict: Mapping of primer names to sequences.
            probe_file: Path to probe FASTA file.
            headers: Column headers to prepend to BLAST outfmt 6 output.
            threads: Number of threads for BLAST.

        Returns:
            Updated metadata list with BLAST result paths and commands.
        """
        logging.info('Running %s probe BLAST analyses', group)
        for sample in metadata:
            for primer in primer_dict:
                try:
                    base_prefix = os.path.join(
                        sample.general.outputdirectory, sample.name
                    )
                    sample[analysistype][primer].blastresults = (
                        f"{base_prefix}_{primer}_probe_blast_results.tsv"
                    )
                    # Try to remove empty reports, and run the blast analyses
                    try:
                        size = os.path.getsize(
                            sample[analysistype][primer].blastresults
                        )
                        if size == 0:
                            os.remove(
                                sample[analysistype][primer].blastresults
                            )
                    except FileNotFoundError:
                        pass
                    # Check to see if the results attribute is empty
                    if not os.path.isfile(
                        sample[analysistype][primer].blastresults
                    ):
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
                            outfmt=(
                                '6 '
                                'qseqid sseqid positive mismatch gaps evalue '
                                'bitscore slen length qstart qend qseq '
                                'sstart send sseq'
                            ),
                            out=sample[analysistype][primer].blastresults,
                        )
                        # Save the blast command in the metadata
                        sample[analysistype][primer].blastcommand = str(blastn)
                        if sample.general.bestassemblyfile != 'NA':
                            # Run the blastn command
                            blastn()
                            # Add a header to the report
                            with open(
                                sample[analysistype][primer].blastresults,
                                'r+',
                                encoding='utf-8',
                            ) as f:
                                # Read contents from the BLAST results file
                                content = f.read()
                                # Go back to the start of the file
                                f.seek(0, 0)
                                # Write header (\n) followed by content
                                f.write('\t'.join(headers) + '\n' + content)
                except AttributeError:
                    pass
        return metadata

    @staticmethod
    def populate_probe_dict(
        *, probe_dict: Dict[str, str], probe_file: str
    ) -> Dict[str, str]:
        """Populate a probe dictionary from a FASTA file.

        Args:
            probe_dict: Mapping of probe name to sequence (will be updated).
            probe_file: Path to probe FASTA file.

        Returns:
            The updated probe_dict.
        """
        # Open the probe file; populate dict with record id -> sequence string
        with open(probe_file, 'r', encoding='utf-8') as probes:
            for record in SeqIO.parse(probes, 'fasta'):
                probe_dict[record.id] = str(record.seq)
        return probe_dict

    @staticmethod
    def parse_blast(
        *,
        metadata: List[MetadataObject],
        analysistype: str,
        primer_dict: Dict[str, str],
        fieldnames: List[str],
        group: str,
        probe_dict: Dict[str, str],
        cutoff: int,
        results_dict: Dict[str, Dict[str, Dict[str, Dict[str, Any]]]],
    ) -> List[MetadataObject]:
        """Parse probe-vs-amplicon BLAST outputs and update metadata.

        Args:
            metadata: Samples in the target panel.
            analysistype: Analysis name (e.g., "probe").
            primer_dict: Mapping of primer names to sequences.
            fieldnames: Headers used for outfmt 6 parsing.
            group: Either "inclusivity" or "exclusivity".
            probe_dict: Mapping of probe names to sequences.
            cutoff: Percent identity cutoff (50–100).
            results_dict: Accumulator for results.

        Returns:
            Updated metadata.
        """
        logging.info('Parsing %s BLAST outputs', group)
        for sample in metadata:
            for primer in primer_dict:
                try:
                    # Only process the samples that were analysed with BLAST
                    if os.path.isfile(
                        sample[analysistype][primer].blastresults
                    ):
                        # Open blast output csv file
                        tsv_file = open(
                            sample[analysistype][primer].blastresults,
                            encoding='utf-8',
                        )
                        # Skip header
                        tsv_file.readline()
                        # Open the sequence profile file as a dictionary
                        blastdict = DictReader(
                            tsv_file,
                            fieldnames=fieldnames,
                            dialect='excel-tab',
                        )
                        hit_dict = dict()
                        # Go through each BLAST result
                        for row in blastdict:
                            # Extract the probe name from the subject_id
                            probe = row['subject_id']
                            positives = int(row['positives'])
                            # sequence_name_primer = row['query_id']
                            # Initialise primer name in the dictionary
                            if primer not in hit_dict:
                                hit_dict[primer] = dict()
                            if probe not in hit_dict[primer]:
                                hit_dict[primer][probe] = dict()
                            # Populate the hit_dict
                            # Find hits that exceed the cutoff
                            percent_id = (
                                positives / len(probe_dict[probe]) * 100
                            )
                            if percent_id > cutoff:
                                row['percent_id'] = percent_id
                                # If dict is empty, add positives -> details
                                if not hit_dict[primer][probe]:
                                    hit_dict[primer][probe][percent_id] = row
                                # Extract positives from previous best in dict
                                previous_best = list(
                                    hit_dict[primer][probe].keys()
                                )[0]
                                # If current positives > previous best
                                if percent_id > previous_best:
                                    # Add current positives -> details
                                    hit_dict[primer][probe][percent_id] = row
                                    # Delete the previous best positives key from the dictionary
                                    del (hit_dict[primer][probe][previous_best])
                        # Iterate through all the probes in the set of probes
                        for probe in probe_dict:
                            try:
                                # Extract the details dictionary for the current probe
                                detail_dict = hit_dict[primer][probe]
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
                            except KeyError:
                                # Populate results_dict with default values for no probe hit
                                if primer not in results_dict:
                                    results_dict[primer] = dict()
                                if group not in results_dict[primer]:
                                    results_dict[primer][group] = dict()
                                if sample.name not in results_dict[primer][group]:
                                    results_dict[primer][group][sample.name] = dict()
                                if 'probe' not in results_dict[primer][group][sample.name]:
                                    results_dict[primer][group][sample.name]['probe'] = dict()
                                if probe not in results_dict[primer][group][sample.name]['probe']:
                                    results_dict[primer][group][sample.name]['probe'][probe] = {
                                        'percent_id': 0.0,
                                        'mismatches': None,
                                        'mismatch_details': None,
                                        'query_sequence': None,
                                        'subject_sequence': None
                                    }
                except AttributeError:
                    pass
        return metadata

    @staticmethod
    def create_probe_report(
        *,
        metadata: List[Any],
        primer_dict: Dict[str, str],
        probe_dict: Dict[str, str],
        group: str,
        report_path: str,
        validator_report_path: str,
        miss_dict: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Create an Excel report for probe BLAST results.

        Args:
            metadata: Samples in the target panel.
            primer_dict: Mapping of primer names to sequences.
            probe_dict: Mapping of probe names to sequences.
            group: Either "inclusivity" or "exclusivity".
            report_path: Output directory for group reports.
            validator_report_path: Directory for aggregated reports.
            miss_dict: Tracker for samples without probe hits.

        Returns:
            The updated miss_dict.
        """
        logging.info('Creating probe report')
        # Prepare the Excel file, and set necessary variables
        (
            report_name,
            report,
            workbook,
            worksheet,
            bold,
            courier,
            row,
            col,
        ) = PrimerValidator.prep_worksheet(
            group=group, report_path=report_path, analysis='probe'
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
        column_width = {
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
            worksheet.set_column(col, col, column_width[col])
            col += 1
        missing = list()
        for sample in metadata:
            # Initialise a list to store all the data for each strain
            for primer in primer_dict:
                for probe in probe_dict:
                    # Next row; reset column to zero for writing results
                    col = 0
                    try:
                        mismatch_details = str()
                        details = sample.probe[primer][probe].details
                        # If mismatches, format details e.g. 32A>T;43A>C
                        if details['percent_id'] < 100:
                            mismatch_details = (
                                PrimerValidator.return_mismatches(
                                    sample=sample, primer=primer, probe=probe
                                )
                            )
                        # Add the appropriate information to the list
                        data = [
                            sample.name,
                            primer,
                            probe,
                            details['percent_id'],
                            details['mismatches'],
                            mismatch_details,
                            details['query_sequence'],
                            details['subject_sequence'],
                        ]
                    # If there were no hits, attributes will not be present
                    except AttributeError:
                        # Init keys in dict storing samples without hits
                        if 'probe' not in miss_dict:
                            miss_dict['probe'] = dict()
                        if primer not in miss_dict['probe']:
                            miss_dict['probe'][primer] = dict()
                        if probe not in miss_dict['probe'][primer]:
                            miss_dict['probe'][primer][probe] = list()
                        # Add the sample name to the list in the dictionary
                        if sample.name not in (
                            miss_dict['probe'][primer][probe]
                        ):
                            miss_dict['probe'][primer][probe].append(
                                sample.name
                            )
                        # Samples with no hits at all lack probe_file attr
                        if not GenObject.isattr(sample.probe, 'probe_file'):
                            # Add the sample name only once
                            if sample.name not in missing:
                                data = [sample.name]
                                missing.append(sample.name)
                            # If sample already present, add an empty list
                            else:
                                data = list()
                        # Add empty list for primer/probe combos without hits
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
        shutil.copyfile(
            report, os.path.join(validator_report_path, report_name)
        )
        return miss_dict

    @staticmethod
    def prep_worksheet(
        *, group: str, report_path: str, analysis: str
    ) -> Tuple[str, str, Any, Any, Any, Any, int, int]:
        """Create a workbook and worksheet, and return formatting helpers.

        Args:
            group: Either "inclusivity" or "exclusivity".
            report_path: Output directory for reports.
            analysis: Either "primer" or "probe" (or variants).

        Returns:
            A tuple of:
              (report_name, report, workbook, worksheet, bold, courier,
               row, col)
        """
        # Create a workbook to store the report (xlsxwriter for sizing)
        report_name = f'{group}_{analysis}_report.xlsx'
        report = os.path.join(
            report_path, report_name
        )
        workbook = xlsxwriter.Workbook(report)
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format(
            {'bold': True, 'font_name': 'Courier New', 'font_size': 10}
        )
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format(
            {'font_name': 'Courier New', 'font_size': 10}
        )
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        return (
            report_name,
            report,
            workbook,
            worksheet,
            bold,
            courier,
            row,
            col,
        )

    @staticmethod
    def return_mismatches(*, sample: Any, primer: str, probe: str) -> str:
        """Return a compact string describing probe mismatches.

        Args:
            sample: Current sample metadata object.
            primer: Primer name.
            probe: Probe name.

        Returns:
            A semicolon-delimited string of mismatches (e.g., "32A>T;43A>C").
        """
        # Initialise a string to store the mismatch details
        mismatches = str()
    # Diff reference: https://bit.ly/3yUMn55
        seq_matcher = difflib.SequenceMatcher(
            None,
            sample.probe[primer][probe].details['query_sequence'],
            sample.probe[primer][probe].details['subject_sequence'])
        for (
            tag,
            query_start,
            query_stop,
            subject_start,
            _,
        ) in seq_matcher.get_opcodes():
            # Only 'replace' tags are interesting
            if tag == 'replace':
                # Append ';' between mismatch entries
                if mismatches:
                    mismatches += ';'
                # Add "posSUBJ>QUERY" (e.g., 32A>T)
                subj_base = sample.probe[primer][probe].details[
                    'subject_sequence'
                ][subject_start]
                query_base = sample.probe[primer][probe].details[
                    'query_sequence'
                ][query_start]
                mismatches += f'{query_stop}{subj_base}>{query_base}'
        return mismatches

    @staticmethod
    def create_excel_report(
        *, report: str, group: str, report_path: str
    ) -> None:
        """Convert a CSV report into a formatted Excel workbook.

        Args:
            report: Path to the CSV report to convert.
            group: Either "inclusivity" or "exclusivity".
            report_path: Output directory for Excel reports.
        """
        # Read in the CSV report
        df = pd.read_csv(report)
        # https://stackoverflow.com/a/61617835
        # Create an pandas.ExcelWriter object
        out_xlsx = os.path.join(report_path, f'{group}_primer_report.xlsx')
        with pd.ExcelWriter(out_xlsx) as writer:
            # Convert CSV-derived dataframe to Excel
            # Set the name of the first sheet to 'Sheet1', do not index, and
            # use 'ND' for NaN cells.
            df.to_excel(writer, sheet_name='Sheet1', index=False, na_rep='ND')
            for column in df:
                # Find the length of the longest cell contents
                column_width = max(
                    df[column].astype(str).map(len).max(), len(str(column))
                )
                # If the width is < 25, set it to 25 to fit headers
                column_width = column_width if column_width >= 25 else 25
                # Find the index of the column
                col_idx = df.columns.get_loc(column)
                # Set the width of the column
                writer.sheets['Sheet1'].set_column(
                    col_idx, col_idx, column_width
                )
        # Context manager closes the writer automatically

    @staticmethod
    def calculate_percentages(
        *, results_dict: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Add positive/total and percent summaries to the results dict.

        Args:
            results_dict: Mapping of primer -> group -> sample -> details.

        Returns:
            The updated results_dict with summary fields.
        """
        for primer, group_dict in results_dict.items():
            for group, sample_dict in group_dict.items():
                # Track positives (primers match) and total samples
                positive = 0
                total = 0
                # Extract the primer SEQID, and primer details dictionary
                for _, contig_results in sample_dict.items():
                    # Increment the total for every sample
                    total += 1
                    # Only increment the positive if there are results
                    if contig_results:
                        positive += 1
                # Update dict with positive and total counts
                results_dict[primer][group]['positive'] = positive
                results_dict[primer][group]['total'] = total
                # Determine the percentage of positives
                try:
                    results_dict[primer][group]['percent'] = (
                        float(positive) / float(total) * 100
                    )
                except ZeroDivisionError:
                    results_dict[primer][group]['percent'] = 0.0
                logging.debug(
                    'Primer set %s in group %s had %s, out of %s, for a '
                    'percentage of %s',
                    primer,
                    group,
                    results_dict[primer][group]['positive'],
                    results_dict[primer][group]['total'],
                    results_dict[primer][group]['percent'],
                )
        # Add the version of the pipeline to the outputs
        results_dict['validator_version'] = __version__
        return results_dict

    def create_summary_report(self) -> None:
        """
        Create a summary Excel workbook at self.report_path/summary_report.xlsx

        Sheets created:
        - PrimerVerifierResults: per-sample summary of amplicons and mismatches
        - Overview: per-primer/per-panel percent/positive/total rows
        - RawSummaryJSON: flattened JSON view of the summary
        """
        if xlsxwriter is None:
            return

        summary = self.results or {}

        # Build a mapping (PANEL, SEQID) -> primer_sets and per-primer seq data
        entries: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for primer_set_name, primer_data in summary.items():
            if primer_set_name == 'validator_version':
                continue
            if not isinstance(primer_data, dict):
                continue
            for panel_name, panel_dict in primer_data.items():
                if not isinstance(panel_dict, dict):
                    continue
                for seqid, seq_data in panel_dict.items():
                    if seqid in ('percent', 'positive', 'total'):
                        continue
                    key = (panel_name.upper(), seqid)
                    ent = entries.setdefault(
                        key, {"primer_sets": [], "per_primer": {}}
                    )
                    if primer_set_name not in ent["primer_sets"]:
                        ent["primer_sets"].append(primer_set_name)
                    # Store the sequence data for this primer set; use None
                    # for empty entries
                    ent["per_primer"][primer_set_name] = \
                        seq_data if seq_data else None

        # Helper to flatten summary counts for Overview sheet
        def _collect_stats_rows(node, path):
            rows = []
            if isinstance(node, dict):
                if (
                    ("total" in node) or ("positive" in node)
                    or ("percent" in node)
                ):
                    primer = path[0] if len(path) > 0 else ""
                    panel = path[1] if len(path) > 1 else ""
                    genus = path[2] if len(path) > 2 else ""
                    rows.append(
                        (
                            primer,
                            panel,
                            genus,
                            node.get("total"),
                            node.get("positive"),
                            node.get("percent")
                        )
                    )
                for k, v in node.items():
                    rows.extend(_collect_stats_rows(v, path + [str(k)]))
            elif isinstance(node, list):
                for i, v in enumerate(node):
                    rows.extend(_collect_stats_rows(v, path + [str(i)]))
            return rows

        # Helper to emit flattened JSON
        def _emit_flat(ws, node, prefix, row_idx):
            if isinstance(node, dict):
                for k in sorted(node.keys()):
                    v = node[k]
                    next_prefix = prefix + ("." if prefix else "") + str(k)
                    row_idx = _emit_flat(ws, v, next_prefix, row_idx)
            elif isinstance(node, list):
                for i, v in enumerate(node):
                    next_prefix = prefix + ("/" if prefix else "") + str(i)
                    row_idx = _emit_flat(ws, v, next_prefix, row_idx)
            else:
                try:
                    val = json.dumps(node)
                except Exception:
                    val = str(node)
                ws.write(row_idx, 0, prefix)
                ws.write(row_idx, 1, val)
                row_idx += 1
            return row_idx

        tmp = tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False)
        tmp_path = tmp.name
        tmp.close()

        try:
            workbook = xlsxwriter.Workbook(tmp_path)
            header_fmt = workbook.add_format(
                {"bold": True, "bg_color": "#F2F2F2", "border": 1}
            )
            cell_fmt = workbook.add_format({"border": 1})
            pct_fmt = workbook.add_format({"num_format": "0.00%", "border": 1})
            int_fmt = workbook.add_format({"num_format": "0", "border": 1})

            # PrimerVerifierResults (first sheet)
            ws0 = workbook.add_worksheet("PrimerVerifierResults")
            headers0 = [
                "SEQID",
                "PANEL",
                "primer set",
                "amplicon length",
                "contig",
                "location",
                "direction",
                "forward mismatch",
                "forward mismatch details",
                "reverse mismatch",
                "reverse mismatch details",
                "probe mismatches",
                "probe mismatch details",
                "Total Mismatches",
            ]
            for col, h in enumerate(headers0):
                ws0.write(0, col, h, header_fmt)

            # Write rows from entries
            panel_order = {"EXCLUSIVITY": 1, "INCLUSIVITY": 0}
            sorted_keys = sorted(
                entries.keys(),
                key=lambda k: (panel_order.get(k[0], 2), k[1])
            )
            write_row = 1

            def _get_field(sd, keys):
                if not sd or not isinstance(sd, dict):
                    return None
                for k in keys:
                    if k in sd:
                        return sd.get(k)
                    # lower-case fallback
                    lk = k.lower()
                    if lk in sd:
                        return sd.get(lk)
                return None

            for panel, seqid in sorted_keys:
                ent = entries[(panel, seqid)]
                primer_list = sorted(ent["primer_sets"])
                # Primer sets that have seq_data (non-empty)
                present = [p for p in primer_list if ent["per_primer"].get(p)]

                def _write_row_from_sd(sd, primer_cell):
                    # sd may be a dict or None
                    if not sd:
                        amplicon_length = "ND"
                        contig = "ND"
                        location = "ND"
                        direction = "ND"
                        forward_mismatch = "ND"
                        forward_mismatch_details = ""
                        reverse_mismatch = "ND"
                        reverse_mismatch_details = ""
                        probe_mismatches = "ND"
                        probe_mismatch_details = ""
                        total_mismatch = "ND"
                    else:
                        amplicon_length = _get_field(
                            sd,
                            [
                                "amplicon_size",
                                "AmpliconSize",
                                "amplicon_length"
                            ]
                        ) or "ND"
                        contig = _get_field(
                            sd,
                            [
                                "contig",
                                "Contig"
                            ]
                        ) or "ND"
                        location = _get_field(
                            sd,
                            [
                                "location",
                                "amplicon_range",
                                "GenomeLocation"
                            ]
                        ) or "ND"
                        direction = _get_field(
                            sd,
                            [
                                "orientation",
                                "direction"
                            ]
                        ) or "ND"
                        forward_mismatch = _get_field(
                            sd,
                            [
                                "forward_mismatch",
                                "ForwardMismatches"
                            ]
                        ) or "ND"
                        forward_mismatch_details = _get_field(
                            sd,
                            [
                                "forward_mismatch_details",
                                "ForwardMismatchDetails"
                            ]
                        ) or ""
                        reverse_mismatch = _get_field(
                            sd,
                            [
                                "reverse_mismatch",
                                "ReverseMismatches"
                            ]
                        ) or "ND"
                        reverse_mismatch_details = _get_field(
                            sd,
                            [
                                "reverse_mismatch_details",
                                "ReverseMismatchDetails"
                            ]
                        ) or ""

                        # Probe extraction (handle nested/flat forms)
                        probe_mismatches = "ND"
                        probe_mismatch_details = ""
                        probe_data = sd.get(
                            "probe"
                        ) if isinstance(sd, dict) else None
                        if probe_data and isinstance(probe_data, dict):
                            try:
                                if "probe" in probe_data and isinstance(
                                    probe_data["probe"], dict
                                ):
                                    pinfo = probe_data["probe"]
                                else:
                                    kprobe = next(iter(probe_data.keys()))
                                    pinfo = probe_data[kprobe] if isinstance(
                                        probe_data[kprobe], dict
                                    ) else probe_data
                                probe_mismatches = pinfo.get(
                                    "mismatches",
                                    pinfo.get("mismatch", "ND")
                                )
                                probe_mismatch_details = pinfo.get(
                                    "mismatch_details", ""
                                )
                            except Exception:
                                probe_mismatches = probe_data.get(
                                    "mismatches", "ND"
                                )

                        # Compute total mismatches
                        def _parse_mismatch_value(val):
                            if val is None:
                                return None
                            sval = str(val).strip()
                            if sval == "" or sval.upper() == "ND":
                                return None
                            try:
                                return int(sval)
                            except Exception:
                                try:
                                    return int(float(sval))
                                except Exception:
                                    return None

                        fval = _parse_mismatch_value(forward_mismatch)
                        rval = _parse_mismatch_value(reverse_mismatch)
                        pval = _parse_mismatch_value(probe_mismatches)
                        if fval is None and rval is None and pval is None:
                            total_mismatch = "ND"
                        else:
                            total_mismatch = \
                                (fval or 0) + (rval or 0) + (pval or 0)

                    ws0.write(write_row, 0, seqid, cell_fmt)
                    ws0.write(write_row, 1, panel, cell_fmt)
                    ws0.write(write_row, 2, primer_cell, cell_fmt)
                    ws0.write(write_row, 3, amplicon_length, cell_fmt)
                    ws0.write(write_row, 4, contig, cell_fmt)
                    ws0.write(write_row, 5, location, cell_fmt)
                    ws0.write(write_row, 6, direction, cell_fmt)
                    ws0.write(write_row, 7, forward_mismatch, cell_fmt)
                    ws0.write(write_row, 8, forward_mismatch_details, cell_fmt)
                    ws0.write(write_row, 9, reverse_mismatch, cell_fmt)
                    ws0.write(
                        write_row, 10, reverse_mismatch_details, cell_fmt
                    )
                    ws0.write(write_row, 11, probe_mismatches, cell_fmt)
                    ws0.write(write_row, 12, probe_mismatch_details, cell_fmt)
                    ws0.write(write_row, 13, total_mismatch, cell_fmt)

                if present:
                    for primer_name in present:
                        sd = ent["per_primer"][primer_name]
                        _write_row_from_sd(sd, primer_name)
                        write_row += 1
                else:
                    # No primer set had data: write grouped row with ND values
                    primer_sets = ", ".join(primer_list)
                    sd = None
                    _write_row_from_sd(sd, primer_sets)
                    write_row += 1

            ws0.autofilter(0, 0, max(write_row - 1, 0), len(headers0) - 1)
            ws0.freeze_panes(1, 0)
            ws0.set_column(0, 0, 20)
            ws0.set_column(1, 1, 12)
            ws0.set_column(2, 2, 30)
            ws0.set_column(3, 6, 14)
            ws0.set_column(7, 12, 22)
            ws0.set_column(13, 13, 16)

            # Overview sheet
            ws = workbook.add_worksheet("Overview")
            headers = [
                "Primer", "Panel", "Genus", "Total", "Positive", "Percent"
            ]
            for col, h in enumerate(headers):
                ws.write(0, col, h, header_fmt)
            stats_rows = _collect_stats_rows(summary, [])
            seen = set()
            write_row = 1
            for r in stats_rows:
                if r in seen:
                    continue
                seen.add(r)
                ws.write(write_row, 0, r[0], cell_fmt)
                ws.write(write_row, 1, r[1], cell_fmt)
                ws.write(write_row, 2, r[2], cell_fmt)
                if r[3] is None:
                    ws.write_blank(write_row, 3, None, cell_fmt)
                else:
                    ws.write_number(write_row, 3, r[3], int_fmt)
                if r[4] is None:
                    ws.write_blank(write_row, 4, None, cell_fmt)
                else:
                    ws.write_number(write_row, 4, r[4], int_fmt)
                percent_val = r[5]
                if isinstance(percent_val, (int, float)):
                    if percent_val > 1.0:
                        percent_val = percent_val / 100.0
                    ws.write_number(write_row, 5, percent_val, pct_fmt)
                elif percent_val is None:
                    ws.write_blank(write_row, 5, None, cell_fmt)
                else:
                    ws.write(write_row, 5, str(percent_val), cell_fmt)
                write_row += 1

            ws.autofilter(0, 0, max(write_row - 1, 0), len(headers) - 1)
            ws.freeze_panes(1, 0)
            ws.set_column(0, 0, 22)
            ws.set_column(1, 1, 22)
            ws.set_column(2, 2, 18)
            ws.set_column(3, 4, 10)
            ws.set_column(5, 5, 10)

            # RawSummaryJSON
            ws2 = workbook.add_worksheet("RawSummaryJSON")
            ws2.write(0, 0, "Path", header_fmt)
            ws2.write(0, 1, "Value", header_fmt)
            end_row = _emit_flat(ws2, summary, "", 1)
            ws2.autofilter(0, 0, max(end_row - 1, 0), 1)
            ws2.freeze_panes(1, 0)
            ws2.set_column(0, 0, 70)
            ws2.set_column(1, 1, 80)

            workbook.close()

            # Save into report path
            dest = os.path.join(self.report_path, "summary_report.xlsx")
            shutil.copyfile(tmp_path, dest)
        finally:
            try:
                os.remove(tmp_path)
            except Exception:
                pass

    @staticmethod
    def missing_report(
        *,
        miss_dict: Dict[str, Any],
        group: str,
        analysis: str,
        report_path: str,
        primer_dict: Dict[str, str],
    ) -> None:
        """Create an Excel report listing samples with no hits.

        Args:
            miss_dict: Tracker of missing results per primer/probe.
            group: Either "inclusivity" or "exclusivity".
            analysis: "primer" or "probe".
            report_path: Output directory for group reports.
            primer_dict: Mapping of primer names to sequences.
        """
        logging.info(
            'Creating report for samples with no hits for %s %s',
            group,
            analysis,
        )
        # Create the workbook
        (
            _,
            __,
            workbook,
            worksheet,
            bold,
            courier,
            row,
            col,
        ) = PrimerValidator.prep_worksheet(
            group=group,
            report_path=report_path,
            analysis=f'{analysis}_missing',
        )
        # Initialise a list to store the primers/probes in the dictionary
        headers = list()
    # Separate list for probe name only (not primer+probe in header)
        probe_headers = list()
    # If a panel had no sequences, dict may not be initialised
        try:
            if analysis == 'primer':
                for primer in miss_dict[analysis].keys():
                    headers.append(primer)
            # Find the primer/probe combination
            elif analysis == 'probe':
                for primer in sorted(primer_dict):
                    for probe, _ in sorted(
                        miss_dict[analysis][primer].items()
                    ):
                        probe_headers.append(probe)
                        headers.append(f'{primer}_{probe}')
            else:
                logging.error(
                    'Analysis type %s is unsupported. Something has gone '
                    'wrong!',
                    analysis,
                )
                raise SystemExit
        except KeyError:
            pass
        # Sort the headers
        headers = sorted(headers)
        # Write the primer/probe names to the header
        for header in headers:
            worksheet.write(row, col, header, bold)
            # Set column width to 15 (fits OLC SEQIDs)
            worksheet.set_column(col, col, 15)
            col += 1
        # Reset the column to 0 in order to add the data
        col = 0
        # Iterate through all the primers/probes
        if analysis == 'primer':
            for header in headers:
                # Reset the row number to 1 for each list of sample names
                row = 1
                for sample_name in miss_dict[analysis][header]:
                    # Write the sample name to the worksheet
                    worksheet.write(row, col, sample_name, courier)
                    row += 1
                # Increase the column number for each primer/probe
                col += 1
    # Probes have an extra nested dict for primer/probe combos
        elif analysis == 'probe':
            for primer in sorted(primer_dict):
                # Enumerate through the headers to get the index
                for i, header in enumerate(probe_headers):
                    # Ensure that the primer matches the name in the headers
                    if headers[i].startswith(primer):
                        # Reset row number to 1 for each list of sample names
                        row = 1
                        for sample_name in miss_dict[analysis][primer][header]:
                            # Write the sample name to the worksheet
                            worksheet.write(row, col, sample_name, courier)
                            row += 1
                        # Increase the column number for each primer/probe
                        col += 1
        # Close the workbook
        workbook.close()

    def __init__(
        self,
        *,
        inclusion_metadata: List[MetadataObject],
        exclusion_metadata: List[MetadataObject],
        primers: Dict[str, str],
        report_path: str,
        inclusion_report_path: str,
        exclusion_report_path: str,
        probe_file: Optional[str],
        cutoff: int,
    ) -> None:
        self.inclusion_metadata = inclusion_metadata
        self.exclusion_metadata = exclusion_metadata
        self.primer_dict = primers
        if report_path.startswith('~'):
            self.report_path = os.path.abspath(
                os.path.expanduser(os.path.join(report_path))
            )
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
                self.probe_file = os.path.abspath(
                    os.path.expanduser(os.path.join(probe_file))
                )
            else:
                self.probe_file = os.path.abspath(os.path.join(probe_file))
        else:
            self.probe_file = str()
        self.cutoff = cutoff
        try:
            assert 49 < self.cutoff < 101
        except AssertionError as exc:
            logging.error(
                'The supplied cutoff value of %s is not in the acceptable '
                'range of 50 - 100',
                self.cutoff,
            )
            raise SystemExit from exc
        if not os.path.isfile(self.probe_file) and probe_file is not None:
            logging.warning(
                'Could not locate supplied probe file: %s. Please ensure the '
                'name and path are correct.',
                self.probe_file,
            )
            raise SystemExit
        self.json_report = os.path.join(
            self.report_path, 'inclusivity_exclusivity_report.json'
        )
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


def cli() -> None:
    """Entry point for running primer validation from the command line."""
    parser = ArgumentParser(
        description='Perform in silico PCR analyses to validate primer sets'
    )
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
        help=(
            'Name and path of folder in which inclusivity/exclusivity '
            'reports are to be created'
        ),
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
        help=(
            'If ipcress cannot find an amplicon, search the genome for the '
            'primers and return a positive result if they are found on '
            'separate contigs'
        ),
    )
    parser.add_argument(
        '-rb', '--range_buffer',
        type=int,
        default=0,
        help=(
            'Increase the buffer size between amplicons in a sequence. '
            'Useful if you have overlapping primer sets (e.g. vtyper), and '
            'you are looking for the best one. Default is 0'
        ),
    )
    parser.add_argument(
        '-p', '--probe_file',
        type=str,
        default=None,
        help=(
            'Name and path of file containing probe sequences to test '
            'against extracted amplicons'
        ),
    )
    parser.add_argument(
        '-c', '--cutoff',
        type=int,
        default=80,
        help=(
            'Percent cutoff to determine whether a probe hits an amplicon. '
            'Minimum value is 50. Maximum is 100. Default is 80'
        ),
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
    # External library type mismatch: samples is a list at runtime
    arguments.runmetadata.samples = Filer.filer(arguments)  # type: ignore
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
        contigbreaks=arguments.contigbreaks,
    )
    inclusion.main()

    # Exclusion
    arguments.sequencepath = arguments.exclusion_sequencepath
    arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    # External library type mismatch: samples is a list at runtime
    arguments.runmetadata.samples = Filer.filer(arguments)  # type: ignore
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
        contigbreaks=arguments.contigbreaks,
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
        cutoff=arguments.cutoff,
    )
    validator.main()


if __name__ == '__main__':
    cli()
    logging.info('Primer validation complete')
