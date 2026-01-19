# primer_validator

Validate primer sets against inclusion and exclusion genome collections using in silico PCR and BLAST‑based probe evaluation. This tool runs ipcress on FASTA inputs to call amplicons, optionally evaluates probe sequences against extracted amplicons, and produces inclusion/exclusion summary reports and a combined JSON.

---

## Features
- Inclusion and exclusion analyses on FASTA datasets (assembly inputs)
- In silico PCR via ipcress (exonerate) with mismatch control and size bounds
- Optional probe evaluation using BLAST+ against extracted amplicons
- Per‑group CSV reports and an overall JSON summary with hit statistics

---

## Installation

`conda install -c olcbioinformatics in_silico_pcr`

you may need to add the bioconda channel

`conda config --add channels bioconda`

## Inputs

1. Primer pair list (fasta). Primer names have to end with “-F” or “-R”. Note: it is possible to have an integer 
following the direction: `>vtx1a-F1` or `>vtx1a-F` are both acceptable
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

```bash
```bash
conda create -n primer_validator -c olcbioinformatics in_silico_pcr
conda activate primer_validator
```

---

## Usage
Run the validator on inclusion and exclusion FASTA folders. Reports for each set are written under their respective `<sequence_path>/reports`. A combined JSON is written to the user‑supplied `--report_path`.

General help:
```bash
python -m primer_finder.primer_validator -h
```

Basic run:
```bash
python -m primer_finder.primer_validator \
  -i /path/to/inclusion_fastas \
  -e /path/to/exclusion_fastas \
  -r /path/to/validator_reports \
  -pf /path/to/primers.fasta \
  -m 2 \
  -min 0 \
  -max 1500
```

Optional probe evaluation:
```bash
python -m primer_finder.primer_validator \
  -i /path/to/inclusion_fastas \
  -e /path/to/exclusion_fastas \
  -r /path/to/validator_reports \
  -pf /path/to/primers.fasta \
  -p /path/to/probes.fasta \
  -c 80
```

### Arguments
- `-i, --inclusion_sequencepath` Path to inclusion FASTA files (required)
- `-e, --exclusion_sequencepath` Path to exclusion FASTA files (required)
- `-r, --report_path` Destination folder for combined outputs (required)
- `-pf, --primerfile` Primer FASTA; names must end with `-F` or `-R` (required)
- `-m, --mismatches` Allowed mismatches [0–3] (default: 2)
- `-min, --minampliconsize` Minimum amplicon size (default: 0)
- `-max, --maxampliconsize` Maximum amplicon size (default: 1500)
- `-cb, --contigbreaks` Count primers found on separate contigs as positive
- `-rb, --range_buffer` Buffer between overlapping primer sets (default: 0)
- `-p, --probe_file` Optional probe FASTA for BLAST validation
- `-c, --cutoff` Probe hit percent identity cutoff [50–100] (default: 80)
- `-d, --debug` Enable debug logging to terminal

---

## Outputs
- Summary Excel: `<report_path>/summary_report.xlsx`
- Inclusion CSV: `<inclusion_sequencepath>/reports/custom_epcr_report.csv`
- Exclusion CSV: `<exclusion_sequencepath>/reports/custom_epcr_report.csv`
- Combined JSON: `<report_path>/inclusivity_exclusivity_report.json`
- Amplicon FASTA folders created alongside each reports directory

---

## Notes
- Primer FASTA headers must end with `-F` or `-R` (e.g., `>vtx1a-F`, `>vtx1a-R`)
- Input data for this tool is FASTA (assemblies); raw FASTQ is not required
- Ensure BLAST+ and exonerate are installed and available on PATH

---

## License
MIT

---

## Contact
- Author: Adam Koziol
- Email: adam.koziol@inspection.gc.ca