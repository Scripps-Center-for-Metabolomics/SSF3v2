"""
hydrogen_redistribution.py
===========================

The corrected bond-cleavage classifier for hydrogen redistribution during
collision-induced dissociation (CID), as described in the Methods section
"Bond-Cleavage Classifier: Baseline Correction" of the manuscript.

BACKGROUND
----------
An earlier version of this classifier (see `src/rearrangement_classifier_legacy.py`
and `src/rearrangement_classifier_v2_legacy.py`, retained here for
transparency and side-by-side comparison) defined its zero-rearrangement
(delta-H = 0) reference by valence-completing each candidate fragment
independently at the cleavage site. This does not conserve atoms when a
cleavage's two complementary fragments are considered jointly -- e.g.
ethanol's C-C cleavage gives independent delta-H=0 fragments of CH4 + CH4O
(8 H total) against ethanol's own 6 H.

THE CORRECTION
--------------
This module implements two changes over the legacy classifier:

1. **Atom-inventory baseline** (`fragment_counts_atom_inventory`): every atom,
   including the one at the cleavage site, retains exactly its
   precursor-derived implicit hydrogen count. The two complementary
   fragments of any cleavage now always sum back to the precursor's own
   formula.

2. **Dual ionization mechanism**: because the atom-inventory baseline
   fragment is, by construction, an odd-electron (radical) species at the
   cleavage site, two distinct ionization hypotheses are tested for every
   candidate rather than one uniform convention:
     (a) Direct cation/anion formation via loss/gain of a single electron
         (delta-H = 0 only -- this is a single, specific hypothesis, not a
         family to sweep across delta; see module docstring note below).
     (b) Protonation/deprotonation of a genuine closed-shell neutral, swept
         across delta in {-2, -1, 0, +1, +2}, as in the legacy classifier.

VALIDATION
----------
Mechanism (a) was checked against known chemistry: acylium formation from
methyl acetate reproduces the real acetylium cation mass (43.0178 Da)
without an added proton, which the legacy baseline could not do (it required
44.0257 Da, i.e. assumed protonation).

Both the legacy and corrected classifiers were re-validated with a
chance-matching null model (see `convergence_analysis.py`) to confirm the
corrected baseline's candidate space does not inflate matches by chance
beyond the legacy baseline's own null rate.

USAGE
-----
    from hydrogen_redistribution import classify_compound_corrected

    result = classify_compound_corrected(molblock_lines, sdf_properties)
    for peak in result['peaks']:
        print(peak['classification'], peak['intensity'])  # 'direct' | 'rearranged' | 'unexplained' | 'precursor'

Phase 1 (single-bond, non-ring cleavage) and Phase 2 (ring-inclusive,
k=1..3 bond combinations) are both provided; Phase 2 is the scope used for
the manuscript's reported corpus-wide statistics.
"""
from collections import defaultdict

from mol_core import (
    parse_molblock, formula_and_mass, PROTON_MASS, ELECTRON_MASS, VALENCE
)
from rearrangement_classifier_legacy import parse_peak_fields, PPM_TOL, DELTAS
from rearrangement_classifier_v2_legacy import enumerate_cleavage_sets


# ---------------------------------------------------------------------------
# Atom-inventory baseline
# ---------------------------------------------------------------------------
def fragment_counts_atom_inventory(mol, atom_set, delta):
    """Atom counts for a candidate fragment under the corrected, mass-
    conserving baseline. Every atom in atom_set -- including any atom at a
    cleavage site -- retains its full, precursor-derived implicit hydrogen
    count; `delta` is then added/subtracted uniformly on top of that
    baseline. This does not depend on which specific bond(s) were cut, only
    on the final atom set, so it generalizes to both Phase 1 (single-bond)
    and Phase 2 (multi-bond, ring-inclusive) cleavage enumeration without
    modification.
    """
    counts = defaultdict(int)
    for aidx in atom_set:
        atom = mol.atoms[aidx]
        el = atom.element
        if el not in VALENCE:
            return None
        counts[el] += 1
        bond_sum = mol.bond_order_sum(aidx, exclude_bond_idx=None)
        implicit_h = VALENCE[el] - bond_sum
        counts['H'] += max(0, round(implicit_h))
    counts['H'] += delta
    if counts['H'] < 0:
        return None
    return dict(counts)


def _build_candidates(mol, cleavage_iter):
    """Shared candidate-generation logic for both Phase 1 and Phase 2. Yields
    (ion_mass_positive, ion_mass_negative, delta) triples across all
    cleavages and both ionization mechanisms."""
    candidates_pos, candidates_neg = [], []
    for atom_set in cleavage_iter:
        # Mechanism (a): direct ion via electron loss/gain, delta=0 only.
        counts0 = fragment_counts_atom_inventory(mol, atom_set, 0)
        if counts0 is not None:
            _, nm0 = formula_and_mass(counts0)
            if nm0 is not None:
                candidates_pos.append((nm0 - ELECTRON_MASS, 0))
                candidates_neg.append((nm0 + ELECTRON_MASS, 0))
        # Mechanism (b): protonated/deprotonated neutral, full delta sweep.
        for delta in DELTAS:
            counts = fragment_counts_atom_inventory(mol, atom_set, delta)
            if counts is None:
                continue
            _, neutral_mass = formula_and_mass(counts)
            if neutral_mass is None:
                continue
            candidates_pos.append((neutral_mass + PROTON_MASS, delta))
            candidates_neg.append((neutral_mass - PROTON_MASS, delta))
    return candidates_pos, candidates_neg


def _classify_with_candidates(mol, props, candidates_pos, candidates_neg, ppm_tol):
    peak_groups = parse_peak_fields(props)
    if not peak_groups:
        return None

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

    results = {'peaks': []}
    for polarity, key, peaks in peak_groups:
        cand_list = candidates_pos if polarity == 'POSITIVE' else candidates_neg
        precursor_mz = precursor_mz_pos if polarity == 'POSITIVE' else precursor_mz_neg
        for mz, inten in peaks:
            if inten <= 0:
                continue
            tol_da = mz * ppm_tol * 1e-6
            if precursor_mz is not None and abs(mz - precursor_mz) <= tol_da:
                results['peaks'].append({'classification': 'precursor', 'intensity': inten})
                continue
            matches = [c for c in cand_list if abs(c[0] - mz) <= tol_da]
            if not matches:
                classification = 'unexplained'
            else:
                direct_matches = [c for c in matches if c[1] == 0]
                classification = 'direct' if direct_matches else 'rearranged'
            results['peaks'].append({'classification': classification, 'intensity': inten})
    return results


def classify_compound_corrected(molblock_lines, props, ppm_tol=PPM_TOL):
    """Phase 1 (single-bond, non-ring cleavage) classification under the
    corrected, dual-mechanism atom-inventory baseline."""
    mol = parse_molblock(molblock_lines)

    def cleavage_iter():
        for bi in range(len(mol.bonds)):
            if not mol.is_bridge(bi):
                continue
            b = mol.bonds[bi]
            atoms_a = mol.component_after_removal(bi, b.a1)
            all_atoms = set(range(len(mol.atoms)))
            atoms_b = all_atoms - atoms_a
            yield atoms_a
            yield atoms_b

    candidates_pos, candidates_neg = _build_candidates(mol, cleavage_iter())
    return _classify_with_candidates(mol, props, candidates_pos, candidates_neg, ppm_tol)


def classify_compound_corrected_phase2(molblock_lines, props, ppm_tol=PPM_TOL, max_k=3):
    """Phase 2 (ring-inclusive, k=1..3 bond combinations) classification
    under the corrected, dual-mechanism atom-inventory baseline. This is the
    scope used for the manuscript's reported corpus-wide statistics."""
    mol = parse_molblock(molblock_lines)

    def cleavage_iter():
        for bond_set, (atoms_a, atoms_b) in enumerate_cleavage_sets(mol, max_k=max_k):
            yield atoms_a
            yield atoms_b

    candidates_pos, candidates_neg = _build_candidates(mol, cleavage_iter())
    return _classify_with_candidates(mol, props, candidates_pos, candidates_neg, ppm_tol)


if __name__ == '__main__':
    import sys
    print(__doc__)
    print("This module is a library; import classify_compound_corrected(_phase2) "
          "from your own script, or see notebooks/walkthrough.ipynb for an "
          "end-to-end example.")
    sys.exit(0)
