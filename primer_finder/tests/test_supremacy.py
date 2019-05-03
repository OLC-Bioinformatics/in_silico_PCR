#!/usr/bin/env python 3
from spadespipeline.primer_finder_bbduk import PrimerFinder
from validator_helper import validate
import pytest
import shutil
import os

# Define global variables
test_path = os.path.abspath(os.path.dirname(__file__))
test_file_root = os.path.join(test_path, 'test_data')
test_sequences_path = os.path.join(test_file_root, 'sequences', 'supremacy')
test_desired_output = os.path.join(test_file_root, 'desired_output')
test_primer_folder = os.path.join(test_file_root, 'primer_files')

__author__ = 'adamkoziol'


def create_validator(reference_csv, test_csv):
    global validator, column_list
    column_list = validate.find_all_columns(csv_file=reference_csv,
                                            columns_to_exclude=['Sample', 'Contig', 'GenomeLocation'])
    validator = validate.Validator(reference_csv=reference_csv,
                                   test_csv=test_csv,
                                   column_list=column_list,
                                   identifying_column='Sample')


def headers_same():
    # Check that the reference and test CSV files have the same headers.
    headers_are_same = validator.same_columns_in_ref_and_test()
    assert headers_are_same


def columns_present():
    # Check that all the columns in the column_list are in both the reference and test CSV files.
    columns = validator.all_test_columns_in_ref_and_test()
    assert columns


def samples_present():
    # Check that both test and reference have all the same sample names.
    samples = validator.check_samples_present()
    assert samples


def values():
    # Finally, check that Categorical columns are identical for test and reference, and range values
    # are within acceptable limits
    values_are_ok = validator.check_columns_match()
    assert values_are_ok


def clear_analyses():
    for root, dirs, files in os.walk(test_sequences_path):
        for directory in dirs:
            if os.path.isdir(os.path.join(root, directory)):
                shutil.rmtree(os.path.join(root, directory))


def remove_consolidated_report():
    shutil.rmtree(os.path.join(test_file_root, 'consolidated_report'))


def remove_detailed_reports():
    shutil.rmtree(os.path.join(test_file_root, 'detailed_reports'))


def test_empty():
    with pytest.raises(TypeError):
        pf = PrimerFinder()


def test_o157_lmhyla_stn():
    pf = PrimerFinder(path=test_file_root,
                      sequence_path=test_sequences_path,
                      primer_file=os.path.join(test_primer_folder, 'O157_LMhlyA_stn.txt'),
                      mismatches=0,
                      kmer_length='55',
                      analysistype='O157_LMhylA_stn')
    pf.main()


def test_validate_o157_lmhyla_stn():
    create_validator(reference_csv=os.path.join(test_desired_output, 'O157_LMhylA_stn_report.csv'),
                     test_csv=os.path.join(test_file_root, 'consolidated_report', 'O157_LMhylA_stn_report.csv'))
    headers_same()
    columns_present()
    samples_present()
    values()
    clear_analyses()
    remove_consolidated_report()
    remove_detailed_reports()


def test_degenerate_primers():
    pf = PrimerFinder(path=test_file_root,
                      sequence_path=test_sequences_path,
                      primer_file=os.path.join(test_primer_folder, 'degenerate.txt'),
                      mismatches=0,
                      kmer_length='55',
                      analysistype='degenerate')
    pf.main()


def test_validate_degenerate():
    create_validator(reference_csv=os.path.join(test_desired_output, 'degenerate_report.csv'),
                     test_csv=os.path.join(test_file_root, 'consolidated_report', 'degenerate_report.csv'))
    headers_same()
    columns_present()
    samples_present()
    values()
    clear_analyses()
    remove_consolidated_report()
    remove_detailed_reports()
