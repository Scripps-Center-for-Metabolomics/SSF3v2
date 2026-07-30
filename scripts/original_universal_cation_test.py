import pandas as pd
import numpy as np
from scipy import stats
from mol_core import split_sdf_records, parse_molblock, VALENCE
from rearrangement_classifier import classify_compound, parse_peak_fields

UNIVERSAL_CATIONS = {
    'C3H3+': 39.02293, 'C3H5+': 41.03858, 'C4H5+': 53.03858, 'C4H7+': 55.05423,
    'C5H5+': 65.03858, 'C5H7+': 67.05423, 'C5H9+': 69.06988, 'C6H5+': 77.03858,
    'C6H7+': 79.05423, 'C6H9+': 81.06988, 'C6H11+': 83.08553, 'C7H7+': 91.05423,
    'C7H9+': 93.06988, 'C7H11+': 95.08553, 'C8H9+': 105.06988, 'C8H11+': 107.08553,
    'C8H13+': 109.10118, 'C9H7+': 115.05423, 'C9H11+': 119.08553, 'C10H7+': 127.05423,
}
CATION_PPM_TOL = 20

def match_cation(mz):
    tol = mz * CATION_PPM_TOL * 1e-6
    for name, ref in UNIVERSAL_CATIONS.items():
        if abs(mz - ref) <= tol:
            return name
    return None

with open('/tmp/ssf_pair_ids.txt') as f:
    target_ids = set(line.strip() for line in f if line.strip())

# For each compound, record the SET of universal cation names present anywhere
# in its spectrum (POSITIVE mode only, at any energy -- these are EI/CID low-mass
# aromatic/allylic cations, most relevant in positive mode).
cations_present = {}

for i, (molblock_lines, props) in enumerate(split_sdf_records('METLIN_Core_v12.sdf')):
    metlin_id = props.get('METLIN ID')
    if metlin_id not in target_ids:
        continue
    peak_groups = parse_peak_fields(props)
    found = set()
    for polarity, key, peaks in peak_groups:
        if polarity != 'POSITIVE':
            continue
        for mz, inten in peaks:
            if inten <= 0:
                continue
            name = match_cation(mz)
            if name:
                found.add(name)
    cations_present[metlin_id] = found

print(f"Cation presence computed for {len(cations_present)} compounds")

# Merge into pair data
df = pd.read_csv('/tmp/merged_ssf_rearrangement.csv')
df['id_a'] = df['id_a'].astype(str)
df['id_b'] = df['id_b'].astype(str)

def shared_cation_count(row):
    a = cations_present.get(row['id_a'], set())
    b = cations_present.get(row['id_b'], set())
    return len(a & b)

df['n_shared_cations'] = df.apply(shared_cation_count, axis=1)
df['shares_any_cation'] = (df['n_shared_cations'] > 0).astype(int)

print(f"\nPairs sharing >=1 universal cation: {df['shares_any_cation'].sum()} / {len(df)}")
print(df['n_shared_cations'].value_counts().sort_index())

PRIMARY = 'cosine_20eV_POSITIVE'
sub = df.dropna(subset=[PRIMARY, 'structural_similarity']).copy()

# Model: residual approach. First regress spectral ~ structural_similarity + MW,
# then check whether shares_any_cation predicts the RESIDUAL (i.e., similarity
# above/below what structure+size alone would predict).
def ols(y, X_df):
    X = np.column_stack([np.ones(len(X_df))] + [X_df[c].values for c in X_df.columns])
    yv = y.values
    n, k = X.shape
    beta, _, _, _ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    dof = n - k
    sigma2 = (resid @ resid) / dof
    se = np.sqrt(np.diag(sigma2 * np.linalg.inv(X.T @ X)))
    t = beta / se
    p = 2 * (1 - stats.t.cdf(np.abs(t), dof))
    r2 = 1 - (resid @ resid) / ((yv - yv.mean())**2).sum()
    return dict(beta=beta, se=se, t=t, p=p, r2=r2, n=n, resid=resid, names=['intercept']+list(X_df.columns))

def print_ols(res, label):
    print(f"\n=== {label} (n={res['n']}, R^2={res['r2']:.4f}) ===")
    for name, b, se, t, p in zip(res['names'], res['beta'], res['se'], res['t'], res['p']):
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        print(f"  {name:<28s} beta={b:+.5f}  SE={se:.5f}  t={t:+.2f}  p={p:.4g} {stars}")

# Baseline model
m_base = ols(sub[PRIMARY], sub[['structural_similarity', 'mean_mw', 'mw_gap']])
print_ols(m_base, "Baseline: spectral ~ structural_similarity + MW")

# Full model adding shared-cation indicator
m_cation = ols(sub[PRIMARY], sub[['structural_similarity', 'mean_mw', 'mw_gap', 'shares_any_cation']])
print_ols(m_cation, "+ shares_any_cation (binary: do the two compounds share >=1 universal cation)")

m_ncation = ols(sub[PRIMARY], sub[['structural_similarity', 'mean_mw', 'mw_gap', 'n_shared_cations']])
print_ols(m_ncation, "+ n_shared_cations (count of shared universal cations)")

# Direct comparison: mean residual for cation-sharing pairs vs non-sharing pairs
sub['resid_base'] = m_base['resid']
shares = sub[sub['shares_any_cation'] == 1]['resid_base']
no_share = sub[sub['shares_any_cation'] == 0]['resid_base']
t_res = stats.ttest_ind(shares, no_share, equal_var=False)
print(f"\n=== Direct residual comparison ===")
print(f"Mean residual, pairs SHARING a universal cation (n={len(shares)}): {shares.mean():+.5f}")
print(f"Mean residual, pairs NOT sharing (n={len(no_share)}): {no_share.mean():+.5f}")
print(f"Welch's t-test: t={t_res.statistic:+.3f}, p={t_res.pvalue:.4g}")
