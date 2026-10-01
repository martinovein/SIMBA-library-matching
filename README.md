# SIMBA NNK Variant Counting

This repository contains a Python pipeline for exact raw-count processing of
SIMBA NNK single-end FASTQ data.

## Inspect-first findings

The real inputs contain three workbook sheets: `NNK_library`, `Fixed_controls`,
and `Samples`. The manifest lists 121 unique gzip FASTQ files, and the current
`data/` directory contains the same 121 filenames with no missing, duplicated,
or unexpected files.

Representative FASTQ records show Cutadapt-retained reads in the primary
orientation with:

- retained 5prime anchor: `TATGCAGCCTCTTCATCTGGC`
- retained 3prime anchor: `GGTGGAGGATCCGGAGG`
- NNK coding length: 54 bp, expected read length 92 bp
- fixed-control coding length: 48 bp, expected read length 86 bp

The full Excel flanks are not retained as complete read flanks, so the script
uses the observed retained anchors by default. Both anchors can be overridden.

## Usage

```bash
python scripts/count_nnk_fastq.py \
  --workbook library-description.xlsx \
  --data-dir data \
  --output-dir outputs/nnk_counting \
  --threads 8 \
  --min-base-quality 20 \
  --force
```

Useful validation runs:

```bash
python scripts/count_nnk_fastq.py --dry-run --force
python scripts/count_nnk_fastq.py --max-reads-per-file 1000 --threads 4 --force
```

The script writes only to the output directory and does not modify
`library-description.xlsx` or `data/`.

## Outputs

Primary outputs are written under `outputs/nnk_counting/`:

- DNA, protein, and fixed-control raw count tables
- complete mapping annotations
- validated sample manifest
- unexpected sequence table
- Enrich2 count files/configurations
- `QC_report.xlsx`
- machine-readable run parameters, summary, checkpoints, and log

No pseudocounts, normalization, fitness estimates, or SIMBA scores are produced.
Counts are exact integer read counts after extraction and quality filtering.

## Requirements

Install dependencies with:

```bash
python -m pip install -r requirements.txt
```

Enrich2 was not installed in the inspected environment. The pipeline still
generates Enrich2-oriented files and records parser validation as unavailable.
Install Enrich2 separately in a compatible environment to validate/run configs.

## Tests

```bash
python -m pytest tests
```

