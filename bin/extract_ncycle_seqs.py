#!/usr/bin/env python3
"""
extract_ncycle_seqs.py

Extracts amino acid sequences for nitrogen cycling KO hits from Prokka .faa
files, grouped by gene name (e.g. narG, nifH). Supports two sources:

  --bins_tsv / --bins_faa    : KOfamscan hits + .faa files from MAG bins
  --assembly_tsv / --assembly_faa : KOfamscan hits + .faa from full assembly

Outputs one multi-FASTA per gene, ready for MAFFT alignment.

Headers:  >{source}|{locus_tag}|{ko}|{tier}
  where source = bin name (for MAGs) or "unbinned" (for assembly contigs)
"""
import argparse
import os
import glob
from collections import defaultdict

# Map KO IDs to common gene names for file naming and grouping
KO_TO_GENE = {
    # Nitrogen fixation
    'K02588': 'nifH',
    'K02586': 'nifD',
    'K02591': 'nifK',
    # Nitrification - ammonia oxidation
    'K10944': 'amoA',
    'K10945': 'amoB',
    'K10946': 'amoC',
    'K10535': 'hao',
    # Nitrification - nitrite oxidation / denitrification nitrate reductase
    'K00370': 'narG',
    'K00371': 'narH',
    'K00374': 'narI',
    # Denitrification
    'K00368': 'nirK',
    'K15864': 'nirS',
    'K04561': 'norB',
    'K02305': 'norC',
    'K00376': 'nosZ',
    # DNRA
    'K00362': 'nirB',
    'K00363': 'nirD',
    # Assimilatory nitrate reduction
    'K00367': 'narB',
    'K10534': 'NR',
    'K00372': 'nasA',
    'K00366': 'nirA',
    'K17877': 'NIT6',
    # Anammox
    'K20932': 'hzsA',
    'K20933': 'hzsB',
    'K20934': 'hzsC',
    'K20935': 'hzo',
    'K22896': 'nrfA_anammox_1',
    'K22897': 'nrfA_anammox_2',
    'K22898': 'nrfA_anammox_3',
    'K22899': 'nrfA_anammox_4',
}

KEY_GENES = {'narG', 'nirS', 'nirK', 'nosZ', 'nifH', 'norB', 'hao', 'amoA'}


def parse_faa(filepath):
    """Parse a FASTA file into a dict of {locus_tag: sequence}."""
    seqs = {}
    header = None
    parts = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    seqs[header] = ''.join(parts)
                header = line[1:].split()[0]  # locus tag only
                parts = []
            else:
                parts.append(line)
    if header is not None:
        seqs[header] = ''.join(parts)
    return seqs


def load_faa_index(faa_dir):
    """Build {name: {locus_tag: sequence}} from a Prokka output directory.

    Handles both layouts:
      - bins:     faa_dir/{bin}/{bin}.faa  (subdirectories)
      - assembly: faa_dir/{sample}.faa     (flat)
    """
    index = {}
    # Subdirectory layout (bins)
    for faa_path in sorted(glob.glob(os.path.join(faa_dir, '*', '*.faa'))):
        name = os.path.basename(os.path.dirname(faa_path))
        index[name] = parse_faa(faa_path)
    # Flat layout (assembly)
    for faa_path in sorted(glob.glob(os.path.join(faa_dir, '*.faa'))):
        name = os.path.splitext(os.path.basename(faa_path))[0]
        if name not in index:
            index[name] = parse_faa(faa_path)
    return index


def load_kofamscan_hits(tsv_path, ncycle_kos):
    """Load filtered KOfamscan TSV, keep only nitrogen cycle KOs."""
    hits = []
    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
            if not fields or len(fields) < len(header):
                continue
            row = dict(zip(header, fields))
            if row['ko'] in ncycle_kos:
                hits.append(row)
    return hits


def extract_sequences(hits, faa_index, source_label, gene_seqs):
    """Look up amino acid sequences for KO hits and append to gene_seqs dict.

    source_label: used in FASTA header. For bins this is the bin name itself,
                  for assembly this is 'unbinned'.
    """
    missing = 0
    found = 0
    for hit in hits:
        name = hit['bin']       # bin name or sample name (from .out filename)
        locus_tag = hit['gene']
        ko = hit['ko']
        tier = hit['tier']
        gene_name = KO_TO_GENE[ko]

        if name not in faa_index:
            missing += 1
            continue
        if locus_tag not in faa_index[name]:
            missing += 1
            continue

        seq = faa_index[name][locus_tag]
        # For bins, label = bin name; for assembly, label = "unbinned"
        label = name if source_label == 'bins' else 'unbinned'
        fasta_header = f">{label}|{locus_tag}|{ko}|{tier}"
        gene_seqs[gene_name].append((fasta_header, seq))
        found += 1

    return found, missing


def main():
    parser = argparse.ArgumentParser(
        description='Extract amino acid sequences for nitrogen cycling genes '
                    'from Prokka annotations, using KOfamscan hits as the index.'
    )
    parser.add_argument('--bins_tsv',
                        help='Filtered KOfamscan TSV from MAG bins')
    parser.add_argument('--bins_faa',
                        help='Prokka bins annotation directory')
    parser.add_argument('--assembly_tsv',
                        help='Filtered KOfamscan TSV from full assembly')
    parser.add_argument('--assembly_faa',
                        help='Prokka assembly annotation directory')
    parser.add_argument('-o', '--outdir', required=True,
                        help='Output directory for per-gene FASTA files')
    parser.add_argument('--key_only', action='store_true',
                        help='Only output key marker genes: '
                             'narG, nirS, nirK, nosZ, nifH, norB, hao, amoA')
    args = parser.parse_args()

    if not args.bins_tsv and not args.assembly_tsv:
        parser.error('Provide at least one of --bins_tsv or --assembly_tsv')

    os.makedirs(args.outdir, exist_ok=True)
    ncycle_kos = set(KO_TO_GENE.keys())
    gene_seqs = defaultdict(list)

    # ---- MAG bins ----
    if args.bins_tsv and args.bins_faa:
        hits = load_kofamscan_hits(args.bins_tsv, ncycle_kos)
        faa_index = load_faa_index(args.bins_faa)
        found, missing = extract_sequences(hits, faa_index, 'bins', gene_seqs)
        print(f"Bins: {found} sequences extracted, {missing} not found "
              f"({len(faa_index)} bins loaded)")

    # ---- Full assembly (unbinned + binned contigs) ----
    if args.assembly_tsv and args.assembly_faa:
        hits = load_kofamscan_hits(args.assembly_tsv, ncycle_kos)
        faa_index = load_faa_index(args.assembly_faa)
        found, missing = extract_sequences(hits, faa_index, 'assembly', gene_seqs)
        print(f"Assembly: {found} sequences extracted, {missing} not found")

    # ---- Write per-gene FASTA files ----
    genes_to_write = set(gene_seqs.keys())
    if args.key_only:
        genes_to_write = genes_to_write & KEY_GENES

    written = 0
    for gene_name in sorted(genes_to_write):
        entries = gene_seqs[gene_name]
        if len(entries) < 2:
            print(f"  {gene_name}: {len(entries)} sequence(s) -- skipping "
                  f"(need >= 2 for alignment)")
            continue

        outpath = os.path.join(args.outdir, f"{gene_name}.faa")
        with open(outpath, 'w') as out:
            for header, seq in entries:
                out.write(f"{header}\n{seq}\n")
        print(f"  {gene_name}: {len(entries)} sequences -> {outpath}")
        written += 1

    print(f"\nWrote {written} gene FASTA files to {args.outdir}")


if __name__ == '__main__':
    main()
