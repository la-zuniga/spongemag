#!/usr/bin/env python3
"""
map_kegg_api.py

Takes the filtered KOFamScan TSV and queries the KEGG REST API to:
  1. Map each KO to its KEGG pathways
  2. Fetch all KOs expected in each pathway
  3. Compute per-bin pathway completeness
  4. Output a summary TSV
"""
import argparse
import time
import requests
import pandas as pd
from collections import defaultdict


# ── KEGG REST helpers ──────────────────────────────────────────────────────────

BASE_URL = "https://rest.kegg.jp"


def kegg_get(endpoint, retries=3, delay=0.35):
    """
    Politely query the KEGG REST API.
    KEGG asks for no more than ~3 requests/second from a single IP.
    """
    url = f"{BASE_URL}/{endpoint}"
    for attempt in range(retries):
        try:
            response = requests.get(url, timeout=15)
            if response.status_code == 200:
                time.sleep(delay)
                return response.text
            elif response.status_code == 404:
                return None
        except requests.RequestException as e:
            print(f"  Warning: request failed ({e}), retrying {attempt+1}/{retries}")
            time.sleep(2 ** attempt)
    return None


def parse_kegg_table(text):
    """Parse a two-column tab-separated KEGG response into a list of tuples."""
    if not text:
        return []
    rows = []
    for line in text.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) == 2:
            rows.append((parts[0].strip(), parts[1].strip()))
    return rows


def get_pathways_for_ko(ko_id):
    """Return list of pathway IDs containing this KO."""
    text = kegg_get(f"link/pathway/{ko_id}")
    rows = parse_kegg_table(text)
    return [r[1].replace('path:', '') for r in rows if 'map' in r[1]]


def get_ko_list_for_pathway(pathway_id):
    """Return set of all KO IDs expected in a given pathway."""
    text = kegg_get(f"link/ko/{pathway_id}")
    rows = parse_kegg_table(text)
    return set(r[1].replace('ko:', '') for r in rows)


def get_pathway_name(pathway_id):
    """Return the human-readable name of a pathway."""
    text = kegg_get(f"get/{pathway_id}")
    if not text:
        return pathway_id
    for line in text.split('\n'):
        if line.startswith('NAME'):
            return line.replace('NAME', '').strip().rstrip(' - Reference pathway')
    return pathway_id


# ── Core logic ─────────────────────────────────────────────────────────────────

def load_kofamscan_results(filepath, min_tier='putative'):
    """Load filtered KOFamScan TSV, optionally filtering by tier."""
    df = pd.read_csv(filepath, sep='\t')
    if min_tier == 'confirmed':
        df = df[df['tier'] == 'confirmed']
    df = df.sort_values('score', ascending=False)
    df = df.drop_duplicates(subset=['bin', 'gene', 'ko'])
    return df


def map_pathways(df, target_pathways=None):
    """
    For each unique KO in the results, fetch its KEGG pathways.
    Returns:
        ko_to_pathways  : dict  KO -> list of pathway IDs
        pathway_names   : dict  pathway_id -> human readable name
        pathway_ko_sets : dict  pathway_id -> set of all expected KOs
    """
    unique_kos = df['ko'].unique().tolist()
    print(f"Querying KEGG for {len(unique_kos)} unique KOs...")

    ko_to_pathways = {}
    pathway_names = {}
    pathway_ko_sets = {}

    for i, ko in enumerate(unique_kos):
        print(f"  [{i+1}/{len(unique_kos)}] {ko}", end=' ')
        pathways = get_pathways_for_ko(ko)

        if target_pathways:
            pathways = [p for p in pathways if p in target_pathways]

        ko_to_pathways[ko] = pathways
        print(f"-> {pathways if pathways else 'no map pathways'}")

        for pathway_id in pathways:
            if pathway_id not in pathway_ko_sets:
                print(f"    Fetching KO list for {pathway_id}...")
                pathway_ko_sets[pathway_id] = get_ko_list_for_pathway(pathway_id)
                pathway_names[pathway_id] = get_pathway_name(pathway_id)
                print(
                    f"    {pathway_id}: {pathway_names[pathway_id]} "
                    f"({len(pathway_ko_sets[pathway_id])} KOs total)"
                )

    return ko_to_pathways, pathway_names, pathway_ko_sets


def compute_completeness(df, ko_to_pathways, pathway_names, pathway_ko_sets):
    """
    For each bin, compute how many KOs from each pathway are present.
    Returns a DataFrame with completeness stats per bin per pathway.
    """
    records = []

    for bin_name, bin_df in df.groupby('bin'):
        bin_kos = set(bin_df['ko'].unique())

        bin_pathways = set()
        for ko in bin_kos:
            bin_pathways.update(ko_to_pathways.get(ko, []))

        for pathway_id in bin_pathways:
            expected_kos = pathway_ko_sets[pathway_id]
            present_kos = bin_kos & expected_kos
            completeness = len(present_kos) / len(expected_kos) if expected_kos else 0

            pathway_bin_df = bin_df[bin_df['ko'].isin(present_kos)]
            confirmed_kos = set(pathway_bin_df[pathway_bin_df['tier'] == 'confirmed']['ko'])
            putative_kos = set(pathway_bin_df[pathway_bin_df['tier'] == 'putative']['ko'])

            records.append({
                'bin': bin_name,
                'pathway_id': pathway_id,
                'pathway_name': pathway_names.get(pathway_id, pathway_id),
                'pathway_total_kos': len(expected_kos),
                'present_kos': len(present_kos),
                'confirmed_kos': len(confirmed_kos),
                'putative_kos': len(putative_kos),
                'completeness': round(completeness, 4),
                'present_ko_ids': ','.join(sorted(present_kos)),
                'confirmed_ko_ids': ','.join(sorted(confirmed_kos)),
                'putative_ko_ids': ','.join(sorted(putative_kos)),
            })

    return pd.DataFrame(records).sort_values(
        ['bin', 'completeness'], ascending=[True, False]
    )


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Map KOFamScan results to KEGG pathways via the KEGG REST API.'
    )
    parser.add_argument(
        '-i', '--input', required=True,
        help='Filtered KOFamScan TSV from filter_kofamscan.py'
    )
    parser.add_argument(
        '-o', '--output', required=True,
        help='Output pathway completeness TSV'
    )
    parser.add_argument(
        '--min_tier', default='putative',
        choices=['confirmed', 'putative'],
        help='Minimum tier to include [default: putative]'
    )
    parser.add_argument(
        '--min_completeness', type=float, default=0.1,
        help='Minimum pathway completeness fraction to report (0-1) [default: 0.1]'
    )
    parser.add_argument(
        '--target_pathways', nargs='+', default=None,
        help='Optional pathway IDs to restrict to (e.g. map00910 map00920)'
    )
    args = parser.parse_args()

    df = load_kofamscan_results(args.input, args.min_tier)
    print(
        f"Loaded {len(df)} hits across "
        f"{df['bin'].nunique()} bins and "
        f"{df['ko'].nunique()} unique KOs"
    )

    ko_to_pathways, pathway_names, pathway_ko_sets = map_pathways(
        df,
        target_pathways=set(args.target_pathways) if args.target_pathways else None
    )

    results = compute_completeness(df, ko_to_pathways, pathway_names, pathway_ko_sets)
    results = results[results['completeness'] >= args.min_completeness]

    results.to_csv(args.output, sep='\t', index=False)
    print(f"\nWritten {len(results)} pathway-bin records to {args.output}")

    print("\n── Summary ──────────────────────────────────────────────────────")
    for _, row in results.iterrows():
        print(
            f"  {row['bin']:30s}  {row['pathway_id']}  "
            f"{row['pathway_name'][:40]:40s}  "
            f"{row['present_kos']}/{row['pathway_total_kos']} KOs  "
            f"({row['completeness']*100:.1f}%)"
        )


if __name__ == '__main__':
    main()
