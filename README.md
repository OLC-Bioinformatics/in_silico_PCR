## Description

This tool is an attempt to do _in silico_ PCRs from .fastq(.gz) or .fasta files. 

The script uses bbduk to bait reads containing the primer sequences from the .fastq(.gz) files. It 
performs a second round of baiting on the original .fastq(.gz) files with the newly created baited
.fastq files in order to hopefully have sequence data to span the entire amplicon. Resulting
double-baited read files are assembled into contigs using SPAdes. 
The assemblies are BLASTed against the primer file to determine if both forward and reverse 
primers can be found in a single contig, thus a valid PCR product. PCR product length is also 
reported.

The longer the reads in the FASTQ file(s), the better the assembly and the fewer false negatives.

## External Dependencies
1. conda


## Installation

`conda install -c adamkoziol in_silico_pcr`

## Inputs

1. Primer pair list (fasta). Primer names have to end with “-F” or “-R”. Note: it is possible to have an integer 
following the direction: >vtx1a-F1 or >vtx1a-F are both acceptable
2. Raw reads (FASTQ) or assemblies (FASTA)

## Required arguments

````
primer_finder.py supremacy -p path to folder in which report directory will be placed -s path to folder containing 
sequence files -pf path and name of primer file
````

### Optional arguments

`-m number of mismatches: Number of mismatches allowed [0-3]. Default is 1`

`-n number of threads: Number of threads. Default is the number of cores in the system`

`-k kmer length: The range of kmers used in SPAdes assembly. Default is 55,77,99,127, but you can provide a 
comma-separated list of kmers e.g. 21,33,55,77,99,127 or a single kmer e.g. 33`


## Testing

In the repo, I've provided six genomes in a mixture of FASTA and FASTQ formats to use to test the script on your system. 
The FASTA files are assemblies, while the FASTQ files are pre-baited files to reduce size. 
The report you create should match the one in the 'desired_outputs' folder (there may be small
differences when it comes to the order of genes).

```
git clone https://github.com/adamkoziol/in_silico_PCR.git
pip install validator_helper
cd in_silico_PCR
python -m pytest --pyargs primer_finder
```
