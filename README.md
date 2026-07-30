# SSF3: Convergent Fragmentation and Hydrogen Redistribution Limit Structure–Spectrum Fidelity

Code and processed data supporting the manuscript investigating why tandem
mass spectral similarity is an imperfect proxy for structural relatedness
across the METLIN reference library.

> **Status note:** this repository was assembled to accompany a manuscript
> currently on hold pending review of a companion paper (SSF1). Several
> items below are placeholders that need to be completed before this
> repository is made public or submitted alongside a reviewer response --
> search this README and the codebase for `[CONFIRM ...]` and `PLACEHOLDER`.

---

## 1. The scientific question and principal finding

Molecular networking and spectral-library search rest on an implicit
assumption: tandem mass spectra that resemble one another arise from
structurally related molecules. A companion analysis (SSF1) found this
correspondence -- Structure–Spectrum Fidelity (SSF) -- to be moderate
rather than high. This work identifies and quantifies **two distinct,
compounding mechanisms** behind that gap:

1. **Convergent fragmentation**: chemically diverse, structurally unrelated
   compounds funnel onto a small, shared vocabulary of just **56 recurrent
   product-ion formulas** across five chemically distinct families. This
   inflates apparent spectral similarity between otherwise dissimilar
   compounds.
2. **Hydrogen-redistribution mismatch**: two structurally similar compounds
   can differ substantially in how much of their fragmentation requires net
   hydrogen redistribution (rather than direct, valence-conserving
   cleavage). This mismatch is a significant, independent predictor of
   *reduced* spectral similarity between related compounds -- a finding
   that only emerged after a real baseline error in the original classifier
   was found and corrected (see `scripts/hydrogen_redistribution.py` for
   the full technical account).

Both mechanisms are quantified at full corpus scale (958,450+ compounds,
METLIN 960K) and validated against known chemistry (e.g. acylium cation
formation) and chance-matching null models.

---

## 2. How the 56 ions were discovered and classified

The 56 recurrent ions were identified via a bottom-up discovery procedure
(see `original_universal_cation_test.py` for the original discovery script,
preserved verbatim) and organized into five families by chemical mechanism
of formation:

| Family | Ion type | Polarity | Count |
|---|---|---|---|
| 1 | Aromatic/allylic hydrocarbon cations | Positive | 19 |
| 2 | Iminium/ammonium cations | Positive | 14 |
| 3 | Saturated alkyl carbocations | Positive | 7 |
| 4 | Small heteroatom-containing anions | Negative | 11 |
| 5 | Acylium/oxocarbenium cations | Positive | 5 |

(A 6th Family 5 ion, CHO+, is confirmed but negligible in prevalence and
excluded from the headline "56" count -- see
`recurrent_ion_discovery.py`.)

Each candidate mass match is validated against an isotope-satellite
discriminator (`recurrent_ion_discovery.passes_isotope_discriminator`) to
rule out natural-abundance isotope contamination; this was found to change
prevalence estimates by well under 1 percentage point at the 20 ppm
matching tolerance used throughout.

---

## 3. Analyses corresponding to each manuscript figure and table

| Manuscript item | Script(s) | Processed data |
|---|---|---|
| Figure 1 (raw SSF discontinuity) | `generate_figures.figure1_ssf_discontinuity` | -- (regenerate from your own structural/spectral similarity pairs file) |
| Figure 2 (shared-cation regression) | `convergence_analysis.py`, `generate_figures.figure2_shared_cation_regression` | `data/processed/shared_cation_regression_FINAL.csv` |
| Redistribution-mismatch regression (Results) | `hydrogen_redistribution.py`, `convergence_analysis.ols_residual` | `data/processed/redistribution_mismatch_regression_corrected.csv` |
| Table 1 / Figure 5 (family prevalence, vocabulary) | `recurrent_ion_discovery.py`, `ion_family_analysis.py` | `data/processed/figure5_per_ion_prevalence_and_families.json`, `data/processed/figure5_jenks_natural_breaks.json` |
| Table 2 / Figure 4 (prevalence by collision energy) | `ion_family_analysis.scan_corpus_checkpointed` | `data/processed/corpus_scan__count_histogram_all_energies.json`, `corpus_scan__polarity_split.json`, `corpus_scan__polarity_counts.json` |
| Bond-environment analysis | `hydrogen_redistribution.py` (Phase 2) | `data/processed/corpus_scan__bond_env_phase2.json` |
| Non-convergent-compound characterization (Discussion) | `ion_family_analysis.py` (peak count / MW / polarity profiling) | `data/processed/corpus_scan__zero_ion_profile.json`, `corpus_scan__ceiling_check.json` |
| Universal Neutral Losses section | **`scripts/neutral_loss_analysis.py` is a placeholder -- see file docstring. This code was not developed in the session that produced the rest of this repository and needs to be sourced separately before submission.** | -- |

`figures/` contains the actual published PNGs for Figures 1, 2, and 5, plus
the graphical abstract, for direct reference alongside the code that
produced them.

---

## 4. Installation and running the code

```bash
git clone <repository-url>
cd SSF3-recurrent-fragmentation
pip install -r requirements.txt
```

All analysis scripts live in `scripts/`; shared, dependency-free parsing
utilities (MOL/SDF parsing, the legacy pre-correction classifier kept for
side-by-side comparison) live in `src/`. Scripts add `src/` to their import
path automatically when run from the repository root; if importing from
elsewhere, add both `src/` and `scripts/` to `PYTHONPATH`:

```bash
export PYTHONPATH="$PWD/src:$PWD/scripts:$PYTHONPATH"
```

**Quickest way to see everything work**: open `notebooks/walkthrough.ipynb`,
which runs the full pipeline (classification, ion detection, and the actual
shared-cation regression, using real processed data) against a tiny,
synthetic, freely-redistributable example dataset in under a minute:

```bash
jupyter notebook notebooks/walkthrough.ipynb
```

**To run a real full-corpus scan** (requires the METLIN 960K SDF corpus --
see §5 below), the scanning functions in `ion_family_analysis.py` are
checkpointed and safely resumable:

```python
from ion_family_analysis import scan_corpus_checkpointed

state = scan_corpus_checkpointed(
    sdf_glob='/path/to/metlin_960k_chunks/*.sdf',
    state_path='my_scan_state.json',
    energy_filter='40eV',
)
# Call again (same state_path) to resume; prints "ALL FILES COMPLETE" when done.
```

---

## 5. What public processed data are included

`data/processed/` contains derived, aggregate result files (JSON checkpoint
states and CSVs of computed per-pair regression variables) -- **not** raw
compound structures or spectra. These are small (a few MB total) and
freely redistributable, since they contain only aggregate statistics and
computed similarity/classification values, not the underlying proprietary
structures and spectra themselves.

`data/example/` contains a tiny, hand-written, clearly-synthetic SDF file
(two illustrative compounds, methyl acetate and toluene) for demonstrating
the pipeline without requiring access to real METLIN data.

---

## 6. METLIN data access

**The METLIN 960K reference library itself (raw structures and spectra) is
not included in this repository and cannot be freely redistributed here.**
All full-corpus statistics reported in the manuscript were computed against
this corpus; reproducing them from scratch requires separate access.

[PLACEHOLDER -- confirm and insert the actual access process, e.g.:]
> METLIN data access can be requested at [URL], subject to [terms]. Contact
> [name/email] for questions about data access for reproduction purposes.

---

## 7. Citation and contact

See `CITATION.cff` for the structured citation record.

**[PLACEHOLDER -- manuscript citation once finalized:]**
> [Authors]. [Title]. *[Journal]*, [year]. DOI: [pending]

**Contact:** [PLACEHOLDER -- corresponding author name and email]

---

## Known gaps in this repository (confirm/resolve before publishing)

- `scripts/neutral_loss_analysis.py` is an unimplemented placeholder (see file).
- `notebooks/walkthrough.ipynb` demonstrates the pipeline on synthetic data
  only; a full-corpus reproduction requires separate METLIN access (§5–6).
- `LICENSE` defaults to MIT; confirm against institutional/journal policy.
- `CITATION.cff` and the citation section above have placeholder author,
  DOI, and URL fields.
- Winnie Uritboonthai's and Linh Hoang's sign-off on the dual-mechanism
  classifier design (`hydrogen_redistribution.py`) was flagged repeatedly
  during development as the one significant validation item never formally
  closed out -- confirm this has happened before treating the classifier as
  fully validated.
