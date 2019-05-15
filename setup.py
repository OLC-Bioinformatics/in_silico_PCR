#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="in_silico_pcr",
    version="0.0.2.1",
    scripts=['primer_finder/primer_finder.py'],
    packages=find_packages(),
    author="Adam Koziol",
    author_email="adam.koziol@canada.ca",
    url="https://github.com/adamkoziol/in_silico_PCR",
)
