#!/usr/bin/env python3
from setuptools import setup, find_packages
import os
setup(
    name="in_silico_pcr",
    version="0.0.2.2",
    scripts=[os.path.join('primer_finder', 'primer_finder.py')],
    packages=find_packages(),
    include_package_data=True,
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/adamkoziol/in_silico_PCR",
)
