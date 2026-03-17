#!/usr/bin/env python3
"""
map_kegg_api.py

Takes the filtered KOFamScan TSV and:
  1. Queries the KEGG REST API to map KOs to KEGG pathways
  2. Computes completeness of curated nitrogen cycling reactions
  3. Outputs a KEGG pathway completeness TSV and a nitrogen cycle TSV
"""
import argparse
import time
import requests
import urllib3
import pandas as pd

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

BASE_URL = "https://rest.kegg.jp"


# ── Curated nitrogen cycle reaction definitions ────────────────────────────────
#
# Each step has:
#   kos:      set of KO IDs — any one being present counts as the step being present
#   required: if True, this step must be present for the reaction to be considered
#             functionally complete
#
# Sources: Kuypers et al. 2018 (Nature Reviews Microbiology),
#          Shafer et al. 2022, KEGG map00910, MetaCyc nitrogen pathways

NITROGEN_REACTIONS = {

    'nitrogen_fixation': {
        'description': 'Biological nitrogen fixation (N2 -> NH3)',
        'steps': {
            'nifH (nitrogenase Fe protein)':        {'kos': {'K02588'}, 'required': True},
            'nifD (nitrogenase alpha subunit)':     {'kos': {'K02586'}, 'required': True},
            'nifK (nitrogenase beta subunit)':      {'kos': {'K02591'}, 'required': True},
            'nifE/anfG (alt nitrogenase)':          {'kos': {'K00531'}, 'required': False},
        }
    },

    'nitrification_ammonia_oxidation': {
        'description': 'Ammonia oxidation to nitrite (NH3 -> NO2-)',
        'steps': {
            'amoA (ammonia monooxygenase A)':       {'kos': {'K10944'}, 'required': True},
            'amoB (ammonia monooxygenase B)':       {'kos': {'K10945'}, 'required': True},
            'amoC (ammonia monooxygenase C)':       {'kos': {'K10946'}, 'required': False},
            'hao (hydroxylamine oxidoreductase)':   {'kos': {'K10535'}, 'required': True},
        }
    },

    'nitrification_nitrite_oxidation': {
        'description': 'Nitrite oxidation to nitrate (NO2- -> NO3-)',
        'steps': {
            'nxrA/narG (nitrite oxidoreductase A)': {'kos': {'K00370'}, 'required': True},
            'nxrB/narH (nitrite oxidoreductase B)': {'kos': {'K00371'}, 'required': True},
        }
    },

    'denitrification': {
        'description': 'Denitrification (NO3- -> N2)',
        'steps': {
            'narG/narZ (nitrate reductase alpha)':  {'kos': {'K00370'}, 'required': True},
            'narH/narY (nitrate reductase beta)':   {'kos': {'K00371'}, 'required': True},
            'narI (nitrate reductase gamma)':       {'kos': {'K00374'}, 'required': False},
            'nirK or nirS (nitrite -> NO)':         {'kos': {'K00368', 'K15864'}, 'required': True},
            'norB (NO reductase B)':                {'kos': {'K04561'}, 'required': True},
            'norC (NO reductase C)':                {'kos': {'K02305'}, 'required': False},
            'nosZ (N2O reductase)':                 {'kos': {'K00376'}, 'required': True},
        }
    },

    'dissimilatory_nitrate_reduction': {
        'description': 'DNRA — dissimilatory nitrate reduction to ammonium (NO3- -> NH4+)',
        'steps': {
            'narG/narZ (nitrate reductase alpha)':  {'kos': {'K00370'}, 'required': True},
            'narH/narY (nitrate reductase beta)':   {'kos': {'K00371'}, 'required': True},
            'narI (nitrate reductase gamma)':       {'kos': {'K00374'}, 'required': False},
            'nirB (nitrite reductase large)':       {'kos': {'K00362'}, 'required': True},
            'nirD (nitrite reductase small)':       {'kos': {'K00363'}, 'required': True},
        }
    },

    'assimilatory_nitrate_reduction': {
        'description': 'Assimilatory nitrate reduction (NO3- -> NH4+ for biosynthesis)',
        'steps': {
            'narB/NR (assimilatory nitrate red.)':  {'kos': {'K00367', 'K10534'}, 'required': True},
            'nasA (assimilatory nitrate red.)':     {'kos': {'K00372'},           'required': False},
            'nirA (ferredoxin-nitrite red.)':       {'kos': {'K00366'},           'required': True},
            'nirB (NADH nitrite red. large)':       {'kos': {'K00362'},           'required': False},
            'nirD (NADH nitrite red. small)':       {'kos': {'K00363'},           'required': False},
            'K17877 (NAD(P)H nitrite red.)':        {'kos': {'K17877'},           'required': False},
        }
    },

    'anammox': {
        'description': 'Anaerobic ammonia oxidation (NH4+ + NO2- -> N2)',
        'steps': {
            'hzsA (hydrazine synthase A)':          {'kos': {'K20932'}, 'required': True},
            'hzsB (hydrazine synthase B)':          {'kos': {'K20933'}, 'required': True},
            'hzsC (hydrazine synthase C)':          {'kos': {'K20934'}, 'required': False},
            'hzo (hydrazine dehydrogenase)':        {'kos': {'K20935'}, 'required': True},
            'nrfA/nirS (nitrite red. anammox)':     {'kos': {'K22896', 'K22897',
                                                             'K22898', 'K22899'}, 'required': True},
        }
    },
}


# ── KEGG REST helpers ──────────────────────────────────────────────────────────

def kegg_get(endpoint, retries=3, delay=0.35):
    url = f"{BASE_URL}/{endpoint}"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=15, verify=False)
            if r.status_code == 200:
                time.sleep(delay)
                return r.text
            elif r.status_code == 404:
                return None
        except requests.RequestException as e:
            print(f"  Warning: request failed ({e}), retrying {attempt+1}/{retries}")
            time.sleep(2 ** attempt)
    return None


def parse_kegg_table(text):
    if not text:
        return []
    rows = []
    for line in text.strip().split('\n'):
        parts = line.split('\t')
        if len(parts) == 2:
            rows.append((parts[0].strip(), parts[1].strip()))
    return rows


def get_pathways_for_ko(ko_id):
    text = kegg_get(f"link/pathway/{ko_id}")
    rows = parse_kegg_table(text)
    return [r[1].replace('path:', '') for r in rows if 'map' in r[1]]


def get_ko_list_for_pathway(pathway_id):
    text = kegg_get(f"link/ko/{pathway_id}")
    rows = parse_kegg_table(text)
    return set(r[1].replace('ko:', '') for r in rows)


def get_pathway_name(pathway_id):
    text = kegg_get(f"get/{pathway_id}")
    if not text:
        return pathway_id
    for line in text.split('\n'):
        if line.startswith('NAME'):
            return line.replace('NAME', '').strip().rstrip(' - Reference pathway')
    return pathway_id


# ── KEGG pathway completeness ──────────────────────────────────────────────────

def load_kofamscan_results(filepath, min_tier='putative'):
    df = pd.read_csv(filepath, sep='\t')
    if min_tier == 'confirmed':
        df = df[df['tier'] == 'confirmed']
    df = df.sort_values('score', ascending=False)
    df = df.drop_duplicates(subset=['bin', 'gene', 'ko'])
    return df


def map_kegg_pathways(df, target_pathways=None):
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
                print(f"    {pathway_id}: {pathway_names[pathway_id]} "
                      f"({len(pathway_ko_sets[pathway_id])} KOs total)")

    return ko_to_pathways, pathway_names, pathway_ko_sets


def compute_kegg_completeness(df, ko_to_pathways, pathway_names, pathway_ko_sets):
    records = []
    for bin_name, bin_df in df.groupby('bin'):
        bin_kos = set(bin_df['ko'].unique())
        bin_pathways = set()
        for ko in bin_kos:
            bin_pathways.update(ko_to_pathways.get(ko, []))

        for pathway_id in bin_pathways:
            expected = pathway_ko_sets[pathway_id]
            present = bin_kos & expected
            completeness = len(present) / len(expected) if expected else 0

            sub = bin_df[bin_df['ko'].isin(present)]
            confirmed = set(sub[sub['tier'] == 'confirmed']['ko'])
            putative  = set(sub[sub['tier'] == 'putative']['ko'])

            records.append({
                'bin':              bin_name,
                'pathway_id':       pathway_id,
                'pathway_name':     pathway_names.get(pathway_id, pathway_id),
                'pathway_total_kos': len(expected),
                'present_kos':      len(present),
                'confirmed_kos':    len(confirmed),
                'putative_kos':     len(putative),
                'completeness':     round(completeness, 4),
                'present_ko_ids':   ','.join(sorted(present)),
                'confirmed_ko_ids': ','.join(sorted(confirmed)),
                'putative_ko_ids':  ','.join(sorted(putative)),
            })

    if not records:
        print("  Warning: no KEGG pathway records found.")
        return pd.DataFrame(columns=[
            'bin','pathway_id','pathway_name','pathway_total_kos','present_kos',
            'confirmed_kos','putative_kos','completeness',
            'present_ko_ids','confirmed_ko_ids','putative_ko_ids'
        ])

    return pd.DataFrame(records).sort_values(['bin','completeness'], ascending=[True,False])


# ── Nitrogen cycle reaction completeness ───────────────────────────────────────

def compute_nitrogen_completeness(df):
    """
    Evaluate per-bin completeness of each curated nitrogen cycling reaction.

    A step is 'present' if any of its KOs appear in the bin.
    Completeness = present steps / total steps.
    required_steps_complete = True only if every required=True step is present.
    """
    records = []

    for bin_name, bin_df in df.groupby('bin'):
        bin_kos           = set(bin_df['ko'].unique())
        confirmed_bin_kos = set(bin_df[bin_df['tier'] == 'confirmed']['ko'])
        putative_bin_kos  = set(bin_df[bin_df['tier'] == 'putative']['ko'])

        for reaction_name, reaction in NITROGEN_REACTIONS.items():
            steps = reaction['steps']
            present_steps    = []
            missing_required = []
            step_detail      = []
            present_kos_in_reaction = set()

            for step_name, step in steps.items():
                hit = step['kos'] & bin_kos
                if hit:
                    present_steps.append(step_name)
                    present_kos_in_reaction.update(hit)
                    tier = 'confirmed' if (step['kos'] & confirmed_bin_kos) else 'putative'
                    step_detail.append(f"{step_name}:{','.join(sorted(hit))}({tier})")
                else:
                    if step['required']:
                        missing_required.append(step_name)
                    step_detail.append(f"{step_name}:absent")

            completeness       = len(present_steps) / len(steps) if steps else 0
            required_complete  = len(missing_required) == 0
            confirmed_present  = present_kos_in_reaction & confirmed_bin_kos
            putative_present   = present_kos_in_reaction & putative_bin_kos

            if not present_kos_in_reaction:
                confidence = 'none'
            elif confirmed_present and putative_present:
                confidence = 'mixed'
            elif putative_present:
                confidence = 'putative'
            else:
                confidence = 'confirmed'

            records.append({
                'bin':                    bin_name,
                'reaction':               reaction_name,
                'description':            reaction['description'],
                'total_steps':            len(steps),
                'present_steps':          len(present_steps),
                'completeness':           round(completeness, 4),
                'required_steps_complete': required_complete,
                'confidence':             confidence,
                'missing_required_steps': '; '.join(missing_required) or 'none',
                'present_ko_ids':         ','.join(sorted(present_kos_in_reaction)),
                'confirmed_ko_ids':       ','.join(sorted(confirmed_present)),
                'putative_ko_ids':        ','.join(sorted(putative_present)),
                'step_detail':            ' | '.join(step_detail),
            })

    if not records:
        print("  Warning: no nitrogen cycle records found.")
        return pd.DataFrame(columns=[
            'bin','reaction','description','total_steps','present_steps','completeness',
            'required_steps_complete','confidence','missing_required_steps',
            'present_ko_ids','confirmed_ko_ids','putative_ko_ids','step_detail'
        ])

    return pd.DataFrame(records).sort_values(['bin','completeness'], ascending=[True,False])


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Map KOFamScan results to KEGG pathways and nitrogen cycle reactions.'
    )
    parser.add_argument('-i', '--input',     required=True,
                        help='Filtered KOFamScan TSV from filter_kofamscan.py')
    parser.add_argument('-o', '--output',    required=True,
                        help='Output KEGG pathway completeness TSV')
    parser.add_argument('--nitrogen_output', required=True,
                        help='Output nitrogen cycle reaction completeness TSV')
    parser.add_argument('--min_tier', default='putative',
                        choices=['confirmed','putative'],
                        help='Minimum evidence tier to include [default: putative]')
    parser.add_argument('--min_completeness', type=float, default=0.05,
                        help='Min KEGG pathway completeness to report (0-1) [default: 0.05]')
    parser.add_argument('--target_pathways', nargs='+', default=None,
                        help='Restrict KEGG query to these pathway IDs (e.g. map00910)')
    args = parser.parse_args()

    df = load_kofamscan_results(args.input, args.min_tier)
    print(f"Loaded {len(df)} hits across {df['bin'].nunique()} bins "
          f"and {df['ko'].nunique()} unique KOs")

    # KEGG pathway completeness
    print("\n── KEGG pathway mapping ─────────────────────────────────────────")
    ko_to_pathways, pathway_names, pathway_ko_sets = map_kegg_pathways(
        df,
        target_pathways=set(args.target_pathways) if args.target_pathways else None
    )
    kegg_df = compute_kegg_completeness(df, ko_to_pathways, pathway_names, pathway_ko_sets)
    kegg_df = kegg_df[kegg_df['completeness'] >= args.min_completeness]
    kegg_df.to_csv(args.output, sep='\t', index=False)
    print(f"Written {len(kegg_df)} KEGG pathway-bin records to {args.output}")

    # Nitrogen cycle reaction completeness
    print("\n── Nitrogen cycle reaction completeness ─────────────────────────")
    n_df = compute_nitrogen_completeness(df)
    n_df.to_csv(args.nitrogen_output, sep='\t', index=False)
    print(f"Written {len(n_df)} nitrogen cycle records to {args.nitrogen_output}")

    # Print summary table to stdout
    print(f"\n  {'reaction':<40} {'steps':>6}  {'%':>6}  {'req':>4}  confidence")
    print("  " + "-" * 70)
    for _, row in n_df.iterrows():
        req = '✓' if row['required_steps_complete'] else '✗'
        print(f"  {row['reaction']:<40} "
              f"{row['present_steps']}/{row['total_steps']}  "
              f"{row['completeness']*100:>5.1f}%  "
              f"  {req}   {row['confidence']}")


if __name__ == '__main__':
    main()


































































































































































































































