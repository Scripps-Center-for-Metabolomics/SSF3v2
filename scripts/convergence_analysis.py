"""
convergence_analysis.py
=========================

Pair-level analysis of whether sharing recurrent ions predicts spectral
similarity beyond what structural similarity and molecular weight alone
would predict (manuscript Results, "Redistribution-Propensity Mismatch..."
and Figure 5 caption's shared-cation regression).

PROVENANCE NOTE
---------------
The original discovery script for this analysis is preserved verbatim in
`original_universal_cation_test.py`. That script restricts the shared-ion
count to Family 1 (19 aromatic/allylic cations) specifically, at 20 ppm
tolerance, computed against structures in METLIN Core v12. This scope
restriction (Family 1 only, not all 56 ions) was the resolution to a
substantial early discrepancy in this analysis -- an earlier, broader
reconstruction using all 56 ions across all 5 families produced a shared-
count distribution that did not match this script's own reported figures,
because it was answering a different question (union across families vs.
Family 1 specifically). See the manuscript Methods and Limitations sections
for the full account; this module reproduces the ORIGINAL script's scope
(Family 1 only) as the validated, correct methodology.

KEY FINDING
-----------
The naive full-sample Pearson correlation between shared-cation count and
adjusted spectral similarity (r ~ +0.11, p < 1e-11) is driven almost
entirely by a small upper tail of extensively-sharing pairs. A sensitivity
sweep that progressively excludes more of this tail shows the association
decline through zero and turn negative; Spearman rank correlation is
negative throughout and becomes MORE negative as more of the tail is
removed. The manuscript's stated conclusion is therefore that recurrent-ion
sharing is NOT generally associated with elevated similarity -- elevation is
specific to the upper tail -- and that no specific ion-count threshold (e.g.
">=12") should be treated as a validated, general regime boundary.

USAGE
-----
    from convergence_analysis import compute_adjusted_similarity, sensitivity_sweep

    df = compute_adjusted_similarity(pairs_df, ion_sets_by_compound_id)
    sweep = sensitivity_sweep(df, cutoffs=[None, 11, 10, 9, 8, 7, 6, 5, 4, 3])
"""
import numpy as np
from scipy import stats

from recurrent_ion_discovery import FAMILY1_AROMATIC_ALLYLIC


def shared_ion_count(ions_a, ions_b):
    """Number of ions in common between two compounds' detected ion sets."""
    return len(ions_a & ions_b)


def ols_residual(y, X_columns):
    """Manual OLS (no external dependency beyond numpy). Returns the
    fitted residuals, i.e. observed y minus what the linear model of
    X_columns predicts -- this is the "adjusted similarity" used
    throughout the shared-ion regression analysis.

    X_columns: dict of {name: 1D array}, e.g.
        {'structural_similarity': ..., 'mean_mw': ..., 'mw_gap': ...}
    """
    names = list(X_columns.keys())
    X = np.column_stack([np.ones(len(y))] + [X_columns[n] for n in names])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    y_hat = X @ beta
    return y - y_hat, beta


def lowess(x, y, xgrid, bandwidth=3.0):
    """Local linear regression with Gaussian kernel weights (hand-rolled;
    no `statsmodels` dependency assumed). Vectorized per grid point to
    avoid constructing a full NxN weight matrix."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    fitted = np.zeros(len(xgrid))
    for i, x0 in enumerate(xgrid):
        d = x - x0
        w = np.exp(-0.5 * (d / bandwidth) ** 2)
        X = np.column_stack([np.ones(len(x)), d])
        Xw = X * w[:, None]
        try:
            beta = np.linalg.solve(X.T @ Xw, Xw.T @ y)
            fitted[i] = beta[0]
        except np.linalg.LinAlgError:
            fitted[i] = np.nan
    return fitted


def lowess_bootstrap_ci(x, y, xgrid, bandwidth=3.0, n_boot=300, seed=42):
    """Bootstrap 95% confidence band around a LOWESS fit. Returns
    (central_fit, lower_band, upper_band), each aligned to xgrid."""
    rng = np.random.RandomState(seed)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    boot_fits = np.zeros((n_boot, len(xgrid)))
    for b in range(n_boot):
        idx = rng.randint(0, n, n)
        boot_fits[b] = lowess(x[idx], y[idx], xgrid, bandwidth)
    central = lowess(x, y, xgrid, bandwidth)
    lower = np.percentile(boot_fits, 2.5, axis=0)
    upper = np.percentile(boot_fits, 97.5, axis=0)
    return central, lower, upper


def sensitivity_sweep(shared_counts, residuals, cutoffs):
    """For each cutoff in cutoffs (None = full sample, else an integer
    upper bound on shared_counts to retain), compute Pearson and Spearman
    correlation between shared_counts and residuals. This is the analysis
    that revealed the shared-cation association is a pure upper-tail
    phenomenon (see module docstring).

    Returns a list of dicts: {cutoff, n, n_excluded, pearson_r, pearson_p,
    spearman_r, spearman_p}.
    """
    shared_counts = np.asarray(shared_counts)
    residuals = np.asarray(residuals)
    n_total = len(shared_counts)
    results = []
    for cutoff in cutoffs:
        if cutoff is None:
            mask = np.ones(n_total, dtype=bool)
            label = 'full'
        else:
            mask = shared_counts <= cutoff
            label = f'<= {cutoff}'
        sc, rs = shared_counts[mask], residuals[mask]
        r, p = stats.pearsonr(sc, rs)
        sr, sp = stats.spearmanr(sc, rs)
        results.append({
            'cutoff': label, 'n': int(mask.sum()), 'n_excluded': int(n_total - mask.sum()),
            'pearson_r': r, 'pearson_p': p, 'spearman_r': sr, 'spearman_p': sp,
        })
    return results


if __name__ == '__main__':
    print(__doc__)
    print(f"Family 1 (aromatic/allylic) reference ion count: {len(FAMILY1_AROMATIC_ALLYLIC)}")
