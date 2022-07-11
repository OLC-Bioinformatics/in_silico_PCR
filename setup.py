#!/usr/bin/env python
from setuptools import setup, find_packages
from distutils.util import convert_path
import os
__author__ = 'adamkoziol'

# Find the version
version = dict()
with open(convert_path(os.path.join('primer_finder', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)

setup(
    name="in_silico_pcr",
    version=version['__version__'],
    scripts=[
        os.path.join('primer_finder', 'primer_finder.py'),
        os.path.join('primer_finder', 'primer_validator.py')
    ],
    packages=find_packages(),
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/in_silico_PCR",
)
