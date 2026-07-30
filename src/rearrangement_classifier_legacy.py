"""
Rearrangement classifier (rebuild, ring-bond cleavage excluded -- Phase 1).

For each compound:
  1. Parse MOL V2000 block, find non-ring (bridge) bonds.
  2. For each bridge, split into two fragments; compute candidate fragment
     formulas at delta-H in {-2,-1,0,+1,+2} relative to the simple valence-cap
     baseline.
  3. Convert to [fragment+H]+ / [fragment-H]- ion masses.
  4. Match against observed MS/MS peaks (all collision energies present) at a
     given ppm tolerance.
  5. Classify each matched peak: 'direct' (best match has deltaH=0) or
     'rearranged' (only nonzero deltaH matches).

Ring-bond cleavage is deliberately excluded in this phase -- known limitation,
addressed in phase 2.
"""

import sys
import json
from collections import defaultdict
from mol_core import (
    split_sdf_records, parse_molblock, formula_and_mass,
    PROTON_MASS, MONO_MASS, VALENCE
)

PPM_TOL = 20  # matches Winnie's stated METLIN fragment accuracy (0-20 ppm)
DELTAS = [-2, -1, 0, 1, 2]


def parse_peak_fields(props):
    """Return list of (polarity, collision_energy_label, [(mz, intensity), ...])."""
    out = []
    for key, val in props.items():
        if not key.startswith('MASS SPECTRAL PEAKS'):
            continue
        polarity = 'POSITIVE' if key.endswith('POSITIVE') else (
            'NEGATIVE' if key.endswith('NEGATIVE') else None)
        if polarity is None:
            continue
        peaks = []
        for line in val.strip().split('\n'):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                mz = float(parts[0])
                inten = float(parts[1])
            except ValueError:
                continue
            peaks.append((mz, inten))
        out.append((polarity, key, peaks))
    return out


def fragment_counts_at_delta(mol, atom_set, removed_bond_idx, cleave_atom, delta):
    """Element counts for a fragment at a given deltaH relative to the valence-cap baseline."""
    from collections import defaultdict as dd
    counts = dd(int)
    b = mol.bonds[removed_bond_idx]
    removed_order = 1.5 if b.order == 4 else b.order

    for aidx in atom_set:
        atom = mol.atoms[aidx]
        el = atom.element
        if el not in VALENCE:
            return None
        counts[el] += 1
        if aidx == cleave_atom:
            bond_sum = mol.bond_order_sum(aidx, exclude_bond_idx=removed_bond_idx)
            # NOTE: no "+ removed_order" cap here. Adding it double-counts the
            # migrated hydrogen: VALENCE - bond_sum (excluding the broken bond)
            # already gives the valence-satisfied H count for this atom, treating
            # the broken bond's position as if it were simply absent. Delta then
            # sweeps outward from that true zero-rearrangement anchor.
            implicit_h = VALENCE[el] - bond_sum
        else:
            bond_sum = mol.bond_order_sum(aidx, exclude_bond_idx=None)
            implicit_h = VALENCE[el] - bond_sum
        counts['H'] += max(0, round(implicit_h))

    counts['H'] += delta
    if counts['H'] < 0:
        return None
    return dict(counts)


def classify_compound(molblock_lines, props, ppm_tol=PPM_TOL):
    mol = parse_molblock(molblock_lines)
    peak_groups = parse_peak_fields(props)
    if not peak_groups:
        return None

    # Precursor ion mass, to exclude intact-precursor peaks (esp. at 0eV) from
    # being scored as "unmatched fragments" -- they were never fragments.
    all_atoms = set(range(len(mol.atoms)))
    precursor_counts = defaultdict(int)
    ok = True
    for atom in mol.atoms:
        if atom.element not in VALENCE:
            ok = False
            break
        precursor_counts[atom.element] += 1
        bond_sum = mol.bond_order_sum(atom.idx)
        h = VALENCE[atom.element] - bond_sum
        precursor_counts['H'] += max(0, round(h))
    precursor_neutral_mass = None
    if ok:
        _, precursor_neutral_mass = formula_and_mass(dict(precursor_counts))
    precursor_mz_pos = precursor_neutral_mass + PROTON_MASS if precursor_neutral_mass else None
    precursor_mz_neg = precursor_neutral_mass - PROTON_MASS if precursor_neutral_mass else None

    bridges = [bi for bi in range(len(mol.bonds)) if mol.is_bridge(bi)]

    # Precompute all candidate ion masses: list of (mass, deltaH, bridge_idx, side)
    candidates_pos = []  # neutral fragment mass -> will add/subtract proton per polarity
    candidates_neg = []
    frag_cache = {}  # (bridge_idx, side) -> neutral mass at deltaH=0 baseline, for logging

    for bi in bridges:
        b = mol.bonds[bi]
        atoms_a = mol.component_after_removal(bi, b.a1)
        all_atoms = set(range(len(mol.atoms)))
        atoms_b = all_atoms - atoms_a

        for side_name, atom_set, cleave_atom in (
            ('A', atoms_a, b.a1), ('B', atoms_b, b.a2)
        ):
            for delta in DELTAS:
                counts = fragment_counts_at_delta(mol, atom_set, bi, cleave_atom, delta)
                if counts is None:
                    continue
                formula, neutral_mass = formula_and_mass(counts)
                if neutral_mass is None:
                    continue
                candidates_pos.append((neutral_mass + PROTON_MASS, delta, bi, side_name, formula))
                candidates_neg.append((neutral_mass - PROTON_MASS, delta, bi, side_name, formula))

    results = {
        'n_bridges': len(bridges),
        'n_bonds_total': len(mol.bonds),
        'peaks': []
    }

    for polarity, key, peaks in peak_groups:
        cand_list = candidates_pos if polarity == 'POSITIVE' else candidates_neg
        precursor_mz = precursor_mz_pos if polarity == 'POSITIVE' else precursor_mz_neg
        for mz, inten in peaks:
            tol_da = mz * ppm_tol * 1e-6
            if precursor_mz is not None and abs(mz - precursor_mz) <= tol_da:
                results['peaks'].append({
                    'energy_label': key, 'polarity': polarity, 'mz': mz,
                    'intensity': inten, 'classification': 'precursor',
                    'best_delta': None, 'best_formula': None,
                })
                continue
            matches = [c for c in cand_list if abs(c[0] - mz) <= tol_da]
            if not matches:
                classification = 'unmatched'
                best = None
            else:
                direct_matches = [c for c in matches if c[1] == 0]
                if direct_matches:
                    classification = 'direct'
                    best = min(direct_matches, key=lambda c: abs(c[0] - mz))
                else:
                    classification = 'rearranged'
                    best = min(matches, key=lambda c: abs(c[0] - mz))
            results['peaks'].append({
                'energy_label': key, 'polarity': polarity, 'mz': mz,
                'intensity': inten, 'classification': classification,
                'best_delta': best[1] if best else None,
                'best_formula': best[4] if best else None,
            })
    return results


def main(sdf_path, n_compounds, out_path):
    summary = {
        'n_compounds_attempted': 0,
        'n_compounds_ok': 0,
        'n_compounds_error': 0,
        'errors': [],
        'total_intensity_direct': 0.0,
        'total_intensity_rearranged': 0.0,
        'total_intensity_unmatched': 0.0,
        'total_intensity_precursor': 0.0,
        'per_compound': []
    }

    for i, (molblock_lines, props) in enumerate(split_sdf_records(sdf_path)):
        if i >= n_compounds:
            break
        summary['n_compounds_attempted'] += 1
        metlin_id = props.get('METLIN ID', f'idx{i}')
        name = props.get('n', '')
        if molblock_lines is None:
            summary['n_compounds_error'] += 1
            summary['errors'].append({'i': i, 'id': metlin_id, 'error': 'no molblock parsed'})
            continue
        try:
            res = classify_compound(molblock_lines, props)
        except Exception as e:
            summary['n_compounds_error'] += 1
            summary['errors'].append({'i': i, 'id': metlin_id, 'error': str(e)})
            continue

        if res is None:
            summary['n_compounds_error'] += 1
            summary['errors'].append({'i': i, 'id': metlin_id, 'error': 'no peak data'})
            continue

        summary['n_compounds_ok'] += 1
        c_direct = sum(p['intensity'] for p in res['peaks'] if p['classification'] == 'direct')
        c_rearr = sum(p['intensity'] for p in res['peaks'] if p['classification'] == 'rearranged')
        c_unmatched = sum(p['intensity'] for p in res['peaks'] if p['classification'] == 'unmatched')
        c_precursor = sum(p['intensity'] for p in res['peaks'] if p['classification'] == 'precursor')
        summary['total_intensity_direct'] += c_direct
        summary['total_intensity_rearranged'] += c_rearr
        summary['total_intensity_unmatched'] += c_unmatched
        summary['total_intensity_precursor'] += c_precursor
        summary['per_compound'].append({
            'i': i, 'id': metlin_id, 'name': name,
            'n_bridges': res['n_bridges'], 'n_bonds_total': res['n_bonds_total'],
            'n_peaks': len(res['peaks']),
            'intensity_direct': c_direct, 'intensity_rearranged': c_rearr,
            'intensity_unmatched': c_unmatched, 'intensity_precursor': c_precursor,
        })

    matched_total = summary['total_intensity_direct'] + summary['total_intensity_rearranged']
    summary['pct_matched_intensity_requiring_rearrangement'] = (
        100.0 * summary['total_intensity_rearranged'] / matched_total if matched_total > 0 else None
    )

    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"Attempted: {summary['n_compounds_attempted']}")
    print(f"OK: {summary['n_compounds_ok']}  Errors: {summary['n_compounds_error']}")
    print(f"Matched intensity - direct: {summary['total_intensity_direct']:.1f}  "
          f"rearranged: {summary['total_intensity_rearranged']:.1f}  "
          f"unmatched: {summary['total_intensity_unmatched']:.1f}  "
          f"precursor: {summary['total_intensity_precursor']:.1f}")
    if summary['pct_matched_intensity_requiring_rearrangement'] is not None:
        print(f"% of MATCHED intensity requiring rearrangement: "
              f"{summary['pct_matched_intensity_requiring_rearrangement']:.1f}%")
    if summary['errors']:
        print(f"\nFirst few errors:")
        for e in summary['errors'][:5]:
            print(f"  {e}")


if __name__ == '__main__':
    sdf_path = sys.argv[1] if len(sys.argv) > 1 else 'METLIN_Core_v12.sdf'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    out = sys.argv[3] if len(sys.argv) > 3 else 'results.json'
    main(sdf_path, n, out)
