"""
ion_family_analysis.py
========================

Corpus-scale analysis functions for the 56-ion recurrent vocabulary defined
in `recurrent_ion_discovery.py`: per-compound ion presence/counting, and
checkpointed full-corpus scanning (resumable across multiple runs, since a
full METLIN 960K pass takes on the order of 10-15 minutes even for a single,
simple per-ion mass-matching pass).

These functions reproduce the manuscript's reported statistics:
  - Per-family prevalence (Table 1)
  - Prevalence by collision energy (Table 2 / Figure 4): 10.7% (0eV),
    26.5% (10eV), 47.3% (20eV), 76.1% (40eV) of the corpus producing at
    least one of the 56 ions
  - Cumulative "at least one, any energy": 80.6%
  - Per-ion ranked prevalence (Figure 5, panel B)

USAGE
-----
    from ion_family_analysis import ions_present_in_compound

    ions = ions_present_in_compound(peak_groups, energy_filter='40eV')
    # -> set of (family_name, ion_name) tuples

For a full corpus scan, see `scan_corpus_checkpointed` and the worked
example in notebooks/walkthrough.ipynb. Each scan is checkpointed to a JSON
state file so it can be safely interrupted and resumed.
"""
import glob
import json
import os
import re
import time

from mol_core import split_sdf_records
from rearrangement_classifier_legacy import parse_peak_fields
from recurrent_ion_discovery import match_universal_ion

ENERGY_RE = re.compile(r'(\d+)eV')


def ions_present_in_compound(peak_groups, energy_filter=None):
    """Given the (polarity, key, peaks) tuples returned by
    `rearrangement_classifier_legacy.parse_peak_fields`, return the set of
    (family, ion_name) tuples detected across all peaks, optionally
    restricted to a single collision energy (e.g. energy_filter='40eV').
    """
    ions_found = set()
    for polarity, key, peaks in peak_groups:
        if energy_filter is not None and energy_filter not in key:
            continue
        for mz, inten in peaks:
            if inten <= 0:
                continue
            tag = match_universal_ion(mz, polarity)
            if tag:
                ions_found.add(tag)
    return ions_found


def scan_corpus_checkpointed(sdf_glob, state_path, time_budget_seconds=260,
                               energy_filter='40eV'):
    """Checkpointed full-corpus scan: for every compound with usable peak
    data at `energy_filter`, record whether it produces at least one of the
    56 recurrent ions, and how many distinct ions it produces. Resumable --
    call repeatedly until it reports "ALL FILES COMPLETE"; each call runs
    for up to time_budget_seconds before saving state and returning.

    Returns the state dict; also written to state_path as JSON.
    """
    def load_state():
        if os.path.exists(state_path):
            with open(state_path) as f:
                return json.load(f)
        return {'file_idx': 0, 'rec_idx': 0, 'n_ok': 0, 'n_with_any': 0,
                'ion_count_histogram': {}}

    def save_state(state):
        with open(state_path, 'w') as f:
            json.dump(state, f)

    files = sorted(glob.glob(sdf_glob))
    if not files:
        raise FileNotFoundError(f"No files matched: {sdf_glob}")
    state = load_state()
    t0 = time.time()

    while state['file_idx'] < len(files):
        if time.time() - t0 > time_budget_seconds:
            break
        fn = files[state['file_idx']]
        resume = state['rec_idx']
        finished = True
        for i, (mb, props) in enumerate(split_sdf_records(fn)):
            if i < resume:
                continue
            if time.time() - t0 > time_budget_seconds:
                state['rec_idx'] = i
                finished = False
                break
            if mb is None:
                state['rec_idx'] = i + 1
                continue
            peak_groups = parse_peak_fields(props)
            if energy_filter is not None:
                peak_groups = [(pol, key, peaks) for pol, key, peaks in peak_groups
                               if energy_filter in key]
            if not peak_groups:
                state['rec_idx'] = i + 1
                continue
            state['n_ok'] += 1
            ions = ions_present_in_compound(peak_groups)
            if ions:
                state['n_with_any'] += 1
            n_ions = str(len(ions))
            state['ion_count_histogram'][n_ions] = state['ion_count_histogram'].get(n_ions, 0) + 1
            state['rec_idx'] = i + 1
        if finished:
            state['file_idx'] += 1
            state['rec_idx'] = 0
        save_state(state)
        if time.time() - t0 > time_budget_seconds:
            break

    print(f"file_idx={state['file_idx']}/{len(files)} rec_idx={state['rec_idx']} "
          f"n_ok={state['n_ok']} n_with_any={state['n_with_any']}")
    if state['file_idx'] >= len(files):
        print("ALL FILES COMPLETE")
        pct = 100 * state['n_with_any'] / state['n_ok'] if state['n_ok'] else 0
        print(f"Final: {state['n_with_any']}/{state['n_ok']} = {pct:.2f}% "
              f"produce >=1 of the 56 recurrent ions at {energy_filter}")
    return state


def jenks_natural_breaks(data, n_classes):
    """Fisher-Jenks optimal 1D classification (dynamic programming), used to
    derive data-driven "dominant / common / long-tail" prevalence zones
    (Figure 5, panel B) rather than an arbitrary threshold. No third-party
    dependency (jenkspy) is required or assumed to be available.

    Returns a list of n_classes+1 class boundary values (ascending).
    """
    data = sorted(data)
    n = len(data)
    mat1 = [[0] * (n_classes + 1) for _ in range(n + 1)]
    mat2 = [[0] * (n_classes + 1) for _ in range(n + 1)]
    for i in range(1, n_classes + 1):
        mat1[1][i] = 1
        mat2[1][i] = 0
        for j in range(2, n + 1):
            mat2[j][i] = float('inf')
    v = 0.0
    for l in range(2, n + 1):
        s1 = s2 = w = 0.0
        for m in range(1, l + 1):
            i3 = l - m + 1
            val = data[i3 - 1]
            s2 += val * val
            s1 += val
            w += 1
            v = s2 - (s1 * s1) / w
            i4 = i3 - 1
            if i4 != 0:
                for j in range(2, n_classes + 1):
                    if mat2[l][j] >= (v + mat2[i4][j - 1]):
                        mat1[l][j] = i3
                        mat2[l][j] = v + mat2[i4][j - 1]
        mat1[l][1] = 1
        mat2[l][1] = v

    k = n
    kclass = [0.0] * (n_classes + 1)
    kclass[n_classes] = data[n - 1]
    kclass[0] = data[0]
    count_num = n_classes
    while count_num >= 2:
        idx = int(mat1[k][count_num]) - 2
        kclass[count_num - 1] = data[idx]
        k = int(mat1[k][count_num]) - 1
        count_num -= 1
    return kclass


if __name__ == '__main__':
    print(__doc__)
