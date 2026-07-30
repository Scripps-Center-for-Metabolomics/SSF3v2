"""
recurrent_ion_discovery.py
===========================

Defines the 56 recurrent ("universal") product-ion formulas across five
chemically distinct families, and the isotope-satellite discriminator used
to validate that a given mass match is a genuine, independently-formed ion
rather than a natural-abundance isotope satellite of a different peak.

FAMILIES
--------
Family 1 -- aromatic/allylic hydrocarbon cations (positive mode, 19 ions)
Family 2 -- iminium/ammonium cations (positive mode, 14 ions)
Family 3 -- saturated alkyl carbocations (positive mode, 7 ions)
Family 4 -- small heteroatom-containing anions (negative mode, 11 ions)
Family 5 -- acylium/oxocarbenium cations (positive mode, 5 prevalent ions;
            a 6th, CHO+, is confirmed but negligible in prevalence and is
            excluded from the headline "56 ions" count used throughout the
            manuscript)

All reference masses are monoisotopic; matching uses a 20 ppm tolerance
(`ION_PPM_TOL`), consistent with the tolerance used in the original
`universal_cation_test.py` discovery script.

ISOTOPE-SATELLITE VALIDATION
-----------------------------
A candidate ion match at mass M could, in principle, be the natural 13C (or
other) isotope satellite of an unrelated peak at M - 1.00336, rather than a
genuine, independently-formed species. `passes_isotope_discriminator` checks
for a peak at that lighter mass and, if present, requires the candidate/
parent intensity ratio to exceed 15% (i.e. larger than a plausible single-
13C natural-abundance satellite) before counting the candidate as genuine.
Applying this filter was found to change corpus-wide prevalence estimates by
well under 1 percentage point (see Methods, "Isotope-Satellite Validation"),
confirming that spurious isotope-satellite contamination is not a material
concern at the 20 ppm tolerance used.

USAGE
-----
    from recurrent_ion_discovery import match_universal_ion, FAMILIES

    tag = match_universal_ion(mz=91.0542, polarity='POSITIVE')
    # -> ('Family1_aromatic', 'C7H7+')
"""

FAMILY1_AROMATIC_ALLYLIC = {  # positive mode, 19 ions
    'C3H3+': 39.02293, 'C3H5+': 41.03858, 'C4H5+': 53.03858, 'C4H7+': 55.05423,
    'C5H5+': 65.03858, 'C5H7+': 67.05423, 'C5H9+': 69.06988, 'C6H5+': 77.03858,
    'C6H7+': 79.05423, 'C6H9+': 81.06988, 'C6H11+': 83.08553, 'C7H7+': 91.05423,
    'C7H9+': 93.06988, 'C7H11+': 95.08553, 'C8H9+': 105.06988, 'C8H11+': 107.08553,
    'C8H13+': 109.10118, 'C9H7+': 115.05423, 'C10H7+': 127.05423,
}
FAMILY2_IMINIUM_AMMONIUM = {  # positive mode, 14 ions
    'C2H6N+': 44.04948, 'C3H6N+': 56.04948, 'C4H8N+': 70.06513,
    'C5H10N+': 84.08078, 'C6H6N+': 92.04948,
    'C2H4N+': 42.0338, 'C3H8N+': 58.0651, 'C4H6N+': 68.0495, 'C5H6N+': 80.0495,
    'C5H8N+': 82.0651, 'C6H8N+': 94.0651, 'C7H6N+': 104.0495, 'C7H8N+': 106.0651,
    'C8H8N+': 118.0651,
}
FAMILY3_ALKYL = {  # positive mode, 7 ions
    'C2H5+': 29.03858, 'C3H7+': 43.05423, 'C4H9+': 57.06988,
    'C5H11+': 71.08553, 'C6H13+': 85.10118, 'C7H15+': 99.11683, 'C8H17+': 113.13248,
}
FAMILY4_ANION = {  # negative mode, 11 ions
    'CNO-': 41.99836, 'C2H3O-': 43.01784, 'C2HO-': 41.00219,
    'HCOO-': 44.99765, 'C2H3O2-': 59.01330,
    'SCN-': 57.97497, 'C7H4NO-': 118.02929, 'C6H5O-': 93.03404,
    'C6H6N-': 92.05005, 'C6H4O-': 92.02620, 'C2H2N-': 40.01870,
}
FAMILY5_ACYLIUM = {  # positive mode, 5 prevalent ions (headline count)
    'C7H5O+': 105.03349, 'C2H3O+': 43.01784, 'C2H5O+': 45.03349,
    'C4H7O+': 71.04914, 'C3H5O+': 57.03349,
}
FAMILY5_ACYLIUM_EXTENDED = dict(FAMILY5_ACYLIUM, **{'CHO+': 29.00220})  # +1 negligible ion

FAMILIES = {
    'Family1_aromatic': (FAMILY1_AROMATIC_ALLYLIC, 'POSITIVE'),
    'Family2_iminium': (FAMILY2_IMINIUM_AMMONIUM, 'POSITIVE'),
    'Family3_alkyl': (FAMILY3_ALKYL, 'POSITIVE'),
    'Family4_anion': (FAMILY4_ANION, 'NEGATIVE'),
    'Family5_acylium': (FAMILY5_ACYLIUM, 'POSITIVE'),
}
N_TOTAL_IONS = sum(len(d) for d, _ in FAMILIES.values())
assert N_TOTAL_IONS == 56, f"Expected 56 headline ions, got {N_TOTAL_IONS}"

ION_PPM_TOL = 20
ISOTOPE_SHIFT = 1.00336  # ~13C - 12C mass difference
ISOTOPE_RATIO_THRESHOLD = 0.15  # candidate/parent ratio above which a match
                                  # is NOT dismissed as a natural isotope satellite


def match_universal_ion(mz, polarity, ppm_tol=ION_PPM_TOL):
    """Return (family_name, ion_name) if mz matches one of the 56 recurrent
    ions at this polarity within tolerance, else None."""
    for fam_name, (ion_dict, fam_polarity) in FAMILIES.items():
        if polarity != fam_polarity:
            continue
        for ion_name, ref in ion_dict.items():
            tol = ref * ppm_tol * 1e-6
            if abs(mz - ref) <= tol:
                return fam_name, ion_name
    return None


def _find_peak(peaks, target_mz, tol_ppm=ION_PPM_TOL):
    """Return the intensity of the peak nearest target_mz within tolerance,
    or None if no such peak exists."""
    tol = target_mz * tol_ppm * 1e-6
    best = None
    for mzv, inten in peaks:
        if inten > 0 and abs(mzv - target_mz) <= tol:
            if best is None or inten > best:
                best = inten
    return best


def passes_isotope_discriminator(candidate_mz, candidate_intensity, same_polarity_peaks,
                                   ratio_threshold=ISOTOPE_RATIO_THRESHOLD):
    """Return True if a candidate ion match at candidate_mz should be counted
    as a genuine, independently-formed ion rather than dismissed as a likely
    natural-abundance isotope satellite of a lighter, unrelated peak.

    If no peak exists at (candidate_mz - ISOTOPE_SHIFT), there is no
    candidate "parent" peak for this to be a satellite of, so the match is
    accepted. If such a parent peak exists, the match is only accepted if
    candidate/parent intensity exceeds ratio_threshold (a simple 13C
    satellite would be a small fraction of its parent's intensity).
    """
    parent_mz = candidate_mz - ISOTOPE_SHIFT
    parent_intensity = _find_peak(same_polarity_peaks, parent_mz)
    if parent_intensity is None:
        return True
    return (candidate_intensity / parent_intensity) > ratio_threshold


if __name__ == '__main__':
    print(__doc__)
    print(f"Total headline ions across 5 families: {N_TOTAL_IONS}")
    for fam, (ions, pol) in FAMILIES.items():
        print(f"  {fam} ({pol}): {len(ions)} ions")
