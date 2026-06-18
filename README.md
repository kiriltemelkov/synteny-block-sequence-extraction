# synteny-block-sequence-extraction
Python algorithm for extracting synteny-block sequences from genome FASTA files using GFF coordinates.

## Overview
This repository contains a Python script that recovers block sequences from genome FASTA files using coordinates stored in a GFF3 file. For each annotated block occurrence, the script
- reads the genomic coordinates and strand,
- extracts the corresponding DNA sequence from the correct genome FASTA file,
- reverse-complements the sequence if the strand is negative,
- groups sequences by block identifier, and
- writes one FASTA file per block.

This step is intended as the extraction stage of a larger pipeline that constructs elastic-degenerate strings from synteny blocks.

## Usage
python extraction.py --gff output.gff --genomes ./genomes  --outdir extracted_blocks

## Command-line arguments

* --gff :  Path to the input GFF3 file
* --genomes :  Path to the directory containing genome FASTA files
* --outdir :  Path to the output directory where block FASTA files will be written
* --feature-type :  Optional filter: process only GFF rows whose third column matches this value exactly

## Requirements
- Python 3.9 or newer
- `pyfastx`

Install the dependency with:
```bash
pip install pyfastx
