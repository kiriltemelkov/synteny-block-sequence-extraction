#!/usr/bin/env python3
"""
Extract synteny block sequences from a GFF file using genome FASTA files

What this algorithm does:
- Parses a GFF3 file.
- For each feature line, reads:
    - seqid
    - start
    - end
    - strand
    - ID
    - genome
- Finds the matching genome FASTA file in a given directory
- Extracts the sequence from the specified seqid and coordinates
- Reverse-complements the sequence if strand == '-'
- Groups sequences by block ID
- Writes one FASTA file per ID into an output directory

Usage: python (name of the algorithm).py --gff output.gff --genomes ./genomes  --outdir extracted_blocks

Needed: pip install pyfastx
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
import pyfastx


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract sequences from genome FASTA files based on GFF coordinates and group them by ID."
    )
    parser.add_argument(
        "--gff",
        required=True,
        help="Input GFF3 file"
    )

    parser.add_argument(
        "--genomes",
        required=True,
        help="Directory containing genome FASTA files (.fna, .fasta, .fa)"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory where one FASTA file per ID will be written"
    )

    parser.add_argument(
        "--feature-type",
        default=None,
        help="Optional: only process rows whose 3rd GFF column matches this value exactly"
    )
    return parser.parse_args()


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]



def parse_gff_attributes(attr_text: str) -> Dict[str, str]:
    """
    Parse GFF attribute column like: ID=5653;genome=GCF_000005845.2_ASM584v2_genomic.fna
    """
    attrs = {}
    parts = attr_text.strip().split(";")
    for part in parts:
        if not part:
            continue
        if "=" in part:
            key, value = part.split("=", 1)
            attrs[key] = value
    return attrs



def build_genome_file_map(genome_dir: Path) -> Dict[str, Path]:
    """
    Build a mapping from basename -> full path for genome FASTA files.

    Supports extensions like: .fna, .fa, .fasta, .fna.gz, .fa.gz, .fasta.gz
    """
    allowed_suffixes = {
        ".fna", ".fa", ".fasta",
        ".fna.gz", ".fa.gz", ".fasta.gz"
    }

    genome_map: Dict[str, Path] = {}

    for path in genome_dir.iterdir():
        if not path.is_file():
            continue

        name = path.name.lower()
        matched = False
        for suffix in allowed_suffixes:
            if name.endswith(suffix):
                genome_map[path.name] = path
                matched = True
                break

        if not matched:
            continue

    return genome_map


def get_pyfastx_handle(cache: Dict[Path, pyfastx.Fasta], fasta_path: Path) -> pyfastx.Fasta:
    if fasta_path not in cache:
        cache[fasta_path] = pyfastx.Fasta(str(fasta_path), build_index=True)
    return cache[fasta_path]


def wrap_fasta_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def sanitize_filename(name: str) -> str:
    """
    Make ID safe as a filename
    """
    safe = []
    for ch in name:
        if ch.isalnum() or ch in ("-", "_", "."):
            safe.append(ch)
        else:
            safe.append("_")
    return "".join(safe)




def main() -> None:
    args = parse_args()

    gff_path = Path(args.gff)
    genome_dir = Path(args.genomes)
    outdir = Path(args.outdir)

    if not gff_path.is_file():
        sys.exit(f"ERROR: GFF file not found: {gff_path}")
    if not genome_dir.is_dir():
        sys.exit(f"ERROR: Genome directory not found: {genome_dir}")

    outdir.mkdir(parents=True, exist_ok=True)

    genome_file_map = build_genome_file_map(genome_dir)
    if not genome_file_map:
        sys.exit(
            "ERROR: No genome FASTA files found in the genomes directory.\n"
            "Expected files ending in .fna, .fa, .fasta, optionally .gz."
        )

    fasta_cache: Dict[Path, pyfastx.Fasta] = {}

    # ID -> list of (header, sequence)
    grouped_sequences: Dict[str, List[Tuple[str, str]]] = defaultdict(list)

    line_number = 0
    extracted_count = 0
    skipped_count = 0



    with gff_path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line_number += 1
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                print(f"WARNING: Skipping malformed GFF line {line_number}: expected 9 columns", file=sys.stderr)
                skipped_count += 1
                continue

            seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes = cols

            if args.feature_type is not None and feature_type != args.feature_type:
                continue

            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                print(f"WARNING: Invalid start/end on line {line_number}", file=sys.stderr)
                skipped_count += 1
                continue

            if start < 1 or end < 1 or end < start:
                print(f"WARNING: Invalid coordinate range on line {line_number}: {start}-{end}", file=sys.stderr)
                skipped_count += 1
                continue

            attrs = parse_gff_attributes(attributes)
            block_id = attrs.get("ID")
            genome_name = attrs.get("genome")

            if not block_id:
                print(f"WARNING: Missing ID attribute on line {line_number}", file=sys.stderr)
                skipped_count += 1
                continue

            if not genome_name:
                print(f"WARNING: Missing genome attribute on line {line_number} for ID={block_id}", file=sys.stderr)
                skipped_count += 1
                continue

            fasta_path = genome_file_map.get(genome_name)
            if fasta_path is None:
                print(
                    f"WARNING: Genome file '{genome_name}' not found in {genome_dir} "
                    f"(line {line_number}, ID={block_id})",
                    file=sys.stderr
                )
                skipped_count += 1
                continue

            try:
                fasta = get_pyfastx_handle(fasta_cache, fasta_path)
            except Exception as e:
                print(f"WARNING: Could not open FASTA file '{fasta_path}': {e}", file=sys.stderr)
                skipped_count += 1
                continue

            try:
                # pyfastx uses 0-based slicing like Python.
                # GFF is 1-based inclusive, so converting:
                # start..end  -> [start-1 : end]
                seq = fasta[seqid][start - 1:end].seq
            except KeyError:
                print(
                    f"WARNING: seqid '{seqid}' not found in genome '{genome_name}' "
                    f"(line {line_number}, ID={block_id})",
                    file=sys.stderr
                )
                skipped_count += 1
                continue
            except Exception as e:
                print(
                    f"WARNING: Could not extract {seqid}:{start}-{end} from '{genome_name}' "
                    f"(line {line_number}, ID={block_id}): {e}",
                    file=sys.stderr
                )
                skipped_count += 1
                continue

            if strand == "-":
                seq = reverse_complement(seq)
            elif strand != "+":
                # If strand is '.' or something else, keep as extracted.
                pass

            header = f"{block_id}|{genome_name}|{seqid}:{start}-{end}({strand})"
            grouped_sequences[block_id].append((header, seq))
            extracted_count += 1



    # Write one FASTA file per ID
    for block_id, records in grouped_sequences.items():
        out_path = outdir / f"{sanitize_filename(block_id)}.fasta"
        with out_path.open("w", encoding="utf-8") as out_fh:
            for header, seq in records:
                out_fh.write(f">{header}\n")
                out_fh.write(wrap_fasta_sequence(seq))
                out_fh.write("\n")

    print(f"Done.")
    print(f"Extracted sequences: {extracted_count}")
    print(f"Skipped lines:       {skipped_count}")
    print(f"IDs written:         {len(grouped_sequences)}")
    print(f"Output directory:    {outdir}")


if __name__ == "__main__":
    main()