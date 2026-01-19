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

### Conda (recommended)
Install external tools from Bioconda, then Python dependencies via pip:

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