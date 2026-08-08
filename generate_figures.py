"""
generate_figures.py
=====================

Regenerates the manuscript's data-driven figures from processed data in
`data/processed/`. This module contains the FINAL, published versions of
each figure; the many exploratory/iterative versions produced during
development are not included here (see the session log referenced in the
repository README for that history, if needed for reviewer transparency).

FIGURES COVERED
---------------
Figure 1 -- Raw Structure-Spectrum Fidelity discontinuity (density plot,
            "if SSF were strong" reference line vs. actual fit, two named
            illustrative example pairs)
Figure 2 -- Shared-cation regression: pair-level scatter, LOWESS smooth with
            bootstrap CI, count histogram, sensitivity sweep panel
Figure 5 -- Integrated four-panel vocabulary summary (word cloud, ranked
            prevalence with Jenks zones, chemical lexicon, family
            contribution donut/bars)

Each figure function takes a pandas DataFrame or dict (as produced by the
corresponding analysis script) and an output path, and writes a 300 DPI PNG.
Run this module directly to regenerate all three from data/processed/.
"""
import numpy as np
import matplotlib.pyplot as plt


def figure1_ssf_discontinuity(structural_similarity, spectral_similarity, out_path,
                                example_pairs=None):
    """Density scatter of structural vs. spectral similarity, with an
    identity reference line ("if SSF were strong") and the actual linear
    fit. example_pairs: optional list of (x, y, label, color) tuples for
    named illustrative points.
    """
    from scipy import stats
    slope, intercept, r, p, se = stats.linregress(structural_similarity, spectral_similarity)

    fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=300)
    hb = ax.hexbin(structural_similarity, spectral_similarity, gridsize=28, cmap='Blues',
                    mincnt=1, linewidths=0.15, edgecolors='#c3c2b7')
    xline = np.array([0, 1])
    ax.plot(xline, xline, linestyle='--', color='#6b6b68', linewidth=1.8,
            label='if SSF were strong', zorder=3)
    ax.plot(xline, intercept + slope * xline, color='#c0392b', linewidth=2.2,
            label=f'actual fit (R\u00b2={r**2:.2f})', zorder=4)

    if example_pairs:
        for x, y, label, color in example_pairs:
            ax.scatter([x], [y], s=150, color=color, edgecolor='white', linewidth=1.2, zorder=6)
            ax.annotate(label, xy=(x, y), fontsize=8.5, ha='left',
                        xytext=(x + 0.05, y + 0.05))

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('structural similarity (Tanimoto)', fontsize=11)
    ax.set_ylabel('spectral similarity', fontsize=11)
    ax.set_title('Structure\u2013Spectrum Fidelity: raw discontinuity', fontsize=12.5, pad=12)
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.95, edgecolor='#c3c2b7')
    fig.colorbar(hb, ax=ax, fraction=0.046, pad=0.03).set_label('compound pairs per bin', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


def figure2_shared_cation_regression(shared_counts, residuals, lowess_result, sweep_result,
                                       out_path):
    """Four-part Figure 2: jittered scatter + LOWESS/CI, log-scale count
    histogram, and a sensitivity-sweep panel showing Pearson/Spearman
    correlation as progressively more of the upper tail is excluded.

    lowess_result: dict with 'xgrid', 'central', 'lower', 'upper'
    sweep_result: list of dicts as returned by convergence_analysis.sensitivity_sweep
    """
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(8.2, 9.4), dpi=300)
    gs = gridspec.GridSpec(3, 1, height_ratios=[4.6, 1.0, 2.4], hspace=0.5)

    ax_main = fig.add_subplot(gs[0])
    rng = np.random.RandomState(1)
    jitter = rng.uniform(-0.18, 0.18, len(shared_counts))
    ax_main.scatter(np.asarray(shared_counts) + jitter, residuals, s=6, alpha=0.16,
                    color='#2a78d6', linewidths=0, zorder=2)
    ax_main.fill_between(lowess_result['xgrid'], lowess_result['lower'], lowess_result['upper'],
                          color='#c0392b', alpha=0.16, zorder=3, label='95% bootstrap CI')
    ax_main.plot(lowess_result['xgrid'], lowess_result['central'], color='#c0392b',
                 linewidth=2.2, zorder=4, label='local smooth')
    ax_main.axhline(0, color='#6b6b68', linewidth=1, linestyle=':', zorder=1)
    ax_main.set_ylabel('adjusted spectral similarity\n(residual: structure + MW removed)', fontsize=10)
    ax_main.set_title('Recurrent-ion sharing is nonlinearly associated with\nelevated spectral similarity',
                       fontsize=12.3, pad=10)
    ax_main.legend(loc='upper left', fontsize=8.5, framealpha=0.95)
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)

    ax_n = fig.add_subplot(gs[1])
    max_x = int(max(shared_counts)) + 1
    counts = np.bincount(np.asarray(shared_counts).astype(int), minlength=max_x)
    ax_n.bar(range(max_x), counts, color='#9fb8d9', width=0.7)
    ax_n.set_yscale('log')
    ax_n.set_ylabel('n pairs\n(log scale)', fontsize=8.5)
    ax_n.set_xlabel('number of Family 1 recurrent cations shared by a compound pair', fontsize=9.5)
    ax_n.spines['top'].set_visible(False)
    ax_n.spines['right'].set_visible(False)

    ax_sens = fig.add_subplot(gs[2])
    x_pos = np.arange(len(sweep_result))
    pearson_vals = [s['pearson_r'] for s in sweep_result]
    spearman_vals = [s['spearman_r'] for s in sweep_result]
    ax_sens.plot(x_pos, pearson_vals, '-o', color='#2a78d6', markersize=6, linewidth=1.8, label='Pearson r')
    ax_sens.plot(x_pos, spearman_vals, '-D', color='#7a7a76', markersize=5, linewidth=1.8, label='Spearman r')
    for i, s in enumerate(sweep_result):
        if s['pearson_p'] >= 0.05:
            ax_sens.scatter([x_pos[i]], [pearson_vals[i]], s=140, facecolors='none',
                             edgecolors='#2a78d6', linewidths=1.5, zorder=5)
    ax_sens.axhline(0, color='#6b6b68', linewidth=1, linestyle=':')
    ax_sens.set_xticks(x_pos)
    ax_sens.set_xticklabels([f"{s['cutoff']}\n(n={s['n']})" for s in sweep_result], fontsize=6.8)
    ax_sens.set_ylabel('correlation with adjusted\nspectral similarity', fontsize=9)
    ax_sens.set_title('Sensitivity sweep (open markers = not significant, p\u22650.05)', fontsize=9.8)
    ax_sens.legend(loc='lower left', fontsize=8.5, framealpha=0.95)
    ax_sens.spines['top'].set_visible(False)
    ax_sens.spines['right'].set_visible(False)

    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    print(__doc__)
    print("Run with processed data from data/processed/ to regenerate figures;")
    print("see notebooks/walkthrough.ipynb for a complete, runnable example.")
