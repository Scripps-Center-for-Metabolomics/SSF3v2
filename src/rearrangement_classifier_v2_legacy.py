"""
Rearrangement classifier v2 -- Phase 2: multi-bond (k=2,3) and ring-inclusive
cleavage, extending the Phase 1 (single-bond, non-ring) classifier.

Cleavage enumeration:
  k=1: existing bridge-bond logic (unchanged from Phase 1).
  k=2, k=3: combinations of bonds drawn from the SAME ring system (a fused-ring
    cluster identified via connectivity among non-bridge bonds), kept only if
    removing them splits the molecule into EXACTLY 2 connected components.
    Combinations that split into >2 (or 0/1) components are skipped rather
    than guessed at -- no neutral-loss or internal-fragment convention is
    assumed, since that wasn't specified anywhere recoverable.

This deliberately does NOT enumerate arbitrary whole-molecule k=2/3 bond
combinations (combinatorially infeasible for large molecules, and mostly
invalid anyway -- cutting two unrelated bridges almost never yields exactly
2 fragments). Restricting to same-ring-system combinations is both
chemically motivated (aromatic/ring losses like C2H2, common in EI/CID
spectra) and computationally tractable, since ring systems are typically
small even in fused polycyclic structures.

Known limitation: does not consider combinations mixing a ring bond with an
adjacent bridge bond. Ring systems larger than MAX_RING_SYSTEM_BONDS bonds
skip k=2/3 enumeration (documented combinatorial ceiling).
"""

import sys
import json
import itertools
from collections import defaultdict
from mol_core import (
    split_sdf_records, parse_molblock, formula_and_mass,
    PROTON_MASS, MONO_MASS, VALENCE
)
from rearrangement_classifier_legacy import parse_peak_fields, PPM_TOL, DELTAS

MAX_RING_SYSTEM_BONDS = 40  # combinatorial ceiling; documented limitation


def fragment_counts_at_delta_multi(mol, atom_set, removed_bond_idxs, delta):
    """Generalized version: removed_bond_idxs is a set (size 1, 2, or 3).
    Every atom in atom_set has all removed bonds touching it excluded when
    computing its implicit-H valence anchor -- atoms untouched by any removed
    bond are unaffected (exclusion set has no effect on them)."""
    counts = defaultdict(int)
    for aidx in atom_set:
        atom = mol.atoms[aidx]
        el = atom.element
        if el not in VALENCE:
            return None
        counts[el] += 1
        bond_sum = mol.bond_order_sum(aidx, exclude_bond_idxs=removed_bond_idxs)
        implicit_h = VALENCE[el] - bond_sum
        counts['H'] += max(0, round(implicit_h))
    counts['H'] += delta
    if counts['H'] < 0:
        return None
    return dict(counts)


def enumerate_cleavage_sets(mol, max_k=3):
    """Yields (frozenset(bond_idxs), (atoms_a, atoms_b)) for every valid
    cleavage: k=1 bridges (always exactly 2 components by definition), plus
    k=2/3 same-ring-system combinations that yield exactly 2 components."""
    # k=1: bridges
    for bi in range(len(mol.bonds)):
        if mol.is_bridge(bi):
            b = mol.bonds[bi]
            atoms_a = mol.component_after_removal(bi, b.a1)
            all_atoms = set(range(len(mol.atoms)))
            atoms_b = all_atoms - atoms_a
            yield frozenset([bi]), (atoms_a, atoms_b)

    if max_k < 2:
        return

    # k=2/3: same-ring-system combinations
    for ring_bonds in mol.ring_systems():
        ring_bonds = sorted(ring_bonds)
        if len(ring_bonds) > MAX_RING_SYSTEM_BONDS:
            continue  # combinatorial ceiling, documented limitation
        for k in range(2, max_k + 1):
            for combo in itertools.combinations(ring_bonds, k):
                comps = mol.components_after_removal(combo)
                if len(comps) == 2:
                    yield frozenset(combo), (comps[0], comps[1])


def classify_compound_v2(molblock_lines, props, ppm_tol=PPM_TOL, max_k=3):
    mol = parse_molblock(molblock_lines)
    peak_groups = parse_peak_fields(props)
    if not peak_groups:
        return None

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

    candidates_pos = []
    candidates_neg = []
    n_cleavage_sets = 0

    for bond_set, (atoms_a, atoms_b) in enumerate_cleavage_sets(mol, max_k=max_k):
        n_cleavage_sets += 1
        for side_name, atom_set in (('A', atoms_a), ('B', atoms_b)):
            for delta in DELTAS:
                counts = fragment_counts_at_delta_multi(mol, atom_set, bond_set, delta)
                if counts is None:
                    continue
                formula, neutral_mass = formula_and_mass(counts)
                if neutral_mass is None:
                    continue
                candidates_pos.append((neutral_mass + PROTON_MASS, delta, formula))
                candidates_neg.append((neutral_mass - PROTON_MASS, delta, formula))

    results = {'n_cleavage_sets': n_cleavage_sets, 'n_bonds_total': len(mol.bonds), 'peaks': []}

    for polarity, key, peaks in peak_groups:
        cand_list = candidates_pos if polarity == 'POSITIVE' else candidates_neg
        precursor_mz = precursor_mz_pos if polarity == 'POSITIVE' else precursor_mz_neg
        for mz, inten in peaks:
            tol_da = mz * ppm_tol * 1e-6
            if precursor_mz is not None and abs(mz - precursor_mz) <= tol_da:
                results['peaks'].append({
                    'energy_label': key, 'polarity': polarity, 'mz': mz,
                    'intensity': inten, 'classification': 'precursor',
                })
                continue
            matches = [c for c in cand_list if abs(c[0] - mz) <= tol_da]
            if not matches:
                classification = 'unmatched'
            else:
                direct_matches = [c for c in matches if c[1] == 0]
                classification = 'direct' if direct_matches else 'rearranged'
            results['peaks'].append({
                'energy_label': key, 'polarity': polarity, 'mz': mz,
                'intensity': inten, 'classification': classification,
            })
    return results


def main(sdf_path, n_compounds, out_path, max_k=3):
    summary = {
        'n_compounds_attempted': 0, 'n_compounds_ok': 0, 'n_compounds_error': 0,
        'errors': [],
        'total_intensity_direct': 0.0, 'total_intensity_rearranged': 0.0,
        'total_intensity_unmatched': 0.0, 'total_intensity_precursor': 0.0,
        'per_compound': []
    }

    for i, (molblock_lines, props) in enumerate(split_sdf_records(sdf_path)):
        if i >= n_compounds:
            break
        summary['n_compounds_attempted'] += 1
        metlin_id = props.get('METLIN ID', f'idx{i}')
        if molblock_lines is None:
            summary['n_compounds_error'] += 1
            summary['errors'].append({'i': i, 'id': metlin_id, 'error': 'no molblock parsed'})
            continue
        try:
            res = classify_compound_v2(molblock_lines, props, max_k=max_k)
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
            'i': i, 'id': metlin_id, 'n_cleavage_sets': res['n_cleavage_sets'],
            'n_bonds_total': res['n_bonds_total'], 'n_peaks': len(res['peaks']),
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


if __name__ == '__main__':
    sdf_path = sys.argv[1] if len(sys.argv) > 1 else 'METLIN_Core_v12.sdf'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    out = sys.argv[3] if len(sys.argv) > 3 else 'results_v2.json'
    max_k = int(sys.argv[4]) if len(sys.argv) > 4 else 3
    main(sdf_path, n, out, max_k=max_k)
