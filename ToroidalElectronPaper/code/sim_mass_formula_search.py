#!/usr/bin/env python3
"""
sim_mass_formula_search.py -- Numerical Search for Electron Mass Formula Parameters

Searches for integer pairs (a, b) such that:
    m_e = m_P * alpha^(a/2 - b*alpha/4)

gives the closest match to the measured electron mass.

The paper (Section 6, Eq. 6.1) reports {a=21, b=15} as the unique pair
giving sub-0.01% accuracy.

Reference:
    - Toroidal_Electron_Full_Paper.md, Section 6.1-6.3
    - m_P = sqrt(hbar * c / G) = 2.176434e-8 kg  (Planck mass)
    - m_e = 9.1093837015e-31 kg  (electron mass, CODATA 2018)
    - alpha = 1/137.035999084  (fine structure constant, CODATA 2018)

Output:
    - Console table of best (a, b) pairs ranked by accuracy
    - Plot: sim_mass_formula_heatmap.png (error heatmap over a, b space)
    - Plot: sim_mass_formula_top_pairs.png (bar chart of top 20 pairs)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

# ─── Color palette (shared with paper figures) ──────────────────
CORAL  = '#e76f51'
TEAL   = '#2a9d8f'
GOLD   = '#e9c46a'
PURPLE = '#a855f7'
BG_COLOR = '#f8f9fa'

# Output directory = images/ folder (sibling to scripts/)
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'images')

# ─── Physical constants (CODATA 2018) ────────────────────────────
HBAR = 1.054571817e-34       # J*s
C    = 2.99792458e8          # m/s
G    = 6.67430e-11           # m^3/(kg*s^2)
M_E  = 9.1093837015e-31      # kg (electron mass)
ALPHA = 1.0 / 137.035999084  # fine structure constant

# Derived
M_P = np.sqrt(HBAR * C / G)  # Planck mass


def mass_formula(a, b, alpha=ALPHA, m_P=M_P):
    """
    Compute predicted electron mass from:
        m_e = m_P * alpha^(a/2 - b*alpha/4)

    Parameters
    ----------
    a, b : int or array
        Integer parameters in the exponent
    alpha : float
        Fine structure constant
    m_P : float
        Planck mass in kg

    Returns
    -------
    m_pred : float or array
        Predicted mass in kg
    """
    exponent = a / 2.0 - b * alpha / 4.0
    return m_P * alpha ** exponent


def relative_error(a, b, m_measured=M_E):
    """Compute relative error |m_pred - m_measured| / m_measured."""
    m_pred = mass_formula(a, b)
    return np.abs(m_pred - m_measured) / m_measured


def search_pairs(a_max=50, b_max=50):
    """
    Exhaustive search over integer pairs (a, b) with
    1 <= a <= a_max, 1 <= b <= b_max.

    Returns
    -------
    results : list of (error, a, b, m_pred, exponent)
        Sorted by ascending error
    error_grid : (a_max, b_max) array
        Error at each (a, b)
    """
    results = []
    error_grid = np.zeros((a_max, b_max))

    for a in range(1, a_max + 1):
        for b in range(1, b_max + 1):
            m_pred = mass_formula(a, b)
            err = abs(m_pred - M_E) / M_E
            exponent = a / 2.0 - b * ALPHA / 4.0
            results.append((err, a, b, m_pred, exponent))
            error_grid[a - 1, b - 1] = err

    results.sort(key=lambda x: x[0])
    return results, error_grid


def extended_search(a_max=100, b_max=200):
    """
    Extended search with wider range, focusing on sub-1% matches.

    Returns list of (error, a, b, m_pred, exponent) with error < 0.01.
    """
    results = []
    for a in range(1, a_max + 1):
        for b in range(1, b_max + 1):
            m_pred = mass_formula(a, b)
            err = abs(m_pred - M_E) / M_E
            if err < 0.01:  # sub-1%
                exponent = a / 2.0 - b * ALPHA / 4.0
                results.append((err, a, b, m_pred, exponent))

    results.sort(key=lambda x: x[0])
    return results


def print_results(results, title, n=30):
    """Pretty-print top N results."""
    print()
    print("=" * 80)
    print(f"  {title}")
    print("=" * 80)
    print(f"  {'Rank':>4}  {'a':>4}  {'b':>4}  {'Exponent':>12}  "
          f"{'m_pred (kg)':>16}  {'Error':>12}")
    print("-" * 80)

    for i, (err, a, b, m_pred, exp) in enumerate(results[:n]):
        flag = " ***" if err < 1e-4 else (" ** " if err < 1e-3 else "")
        print(f"  {i+1:>4}  {a:>4}  {b:>4}  {exp:>12.6f}  "
              f"{m_pred:>16.6e}  {err:>11.6%}{flag}")

    print("-" * 80)
    print(f"  Measured m_e = {M_E:.10e} kg")
    print(f"  Planck mass  = {M_P:.6e} kg")
    print(f"  alpha        = {ALPHA:.12f}")
    print()


def print_analysis(results):
    """Print detailed analysis of the top result."""
    err, a, b, m_pred, exp = results[0]

    print("=" * 80)
    print("  ANALYSIS OF BEST MATCH")
    print("=" * 80)
    print(f"  Best pair: a = {a}, b = {b}")
    print(f"  Exponent: {a}/2 - {b}*alpha/4 = {exp:.10f}")
    print(f"  Predicted mass: {m_pred:.10e} kg")
    print(f"  Measured mass:  {M_E:.10e} kg")
    print(f"  Relative error: {err:.6%}")
    print()
    print(f"  Group-theoretic interpretation:")
    print(f"    a = {a} = 3 x 7 = dim(S^3) x dim(S^7)")
    print(f"       (product of Hopf fibration sphere dimensions)")
    print(f"    b = {b} = dim(SO(4,2)) = dim(conformal group)")
    print()

    # Compare nearest competitors
    print("  Nearest competitors (same a, adjacent b):")
    for db in [-1, 0, 1]:
        b2 = b + db
        if b2 < 1:
            continue
        err2 = relative_error(a, b2)
        label = "  <-- BEST" if db == 0 else ""
        print(f"    (a={a}, b={b2}): error = {err2:.6%}{label}")
    print()

    # Check uniqueness: any other a value with sub-0.01% for any b?
    print("  Uniqueness check (sub-0.01% accuracy, any a, b in [1,50]):")
    count = 0
    for r in results:
        if r[0] < 1e-4:
            count += 1
            print(f"    (a={r[1]}, b={r[2]}): error = {r[0]:.6%}")
    if count == 1:
        print(f"  --> UNIQUE pair with sub-0.01% accuracy")
    else:
        print(f"  --> {count} pairs with sub-0.01% accuracy")
    print()

    # What exponent is needed exactly?
    exact_exp = np.log(M_E / M_P) / np.log(ALPHA)
    print(f"  Exact exponent needed: log_alpha(m_e/m_P) = {exact_exp:.10f}")
    print(f"  Formula gives:         a/2 - b*alpha/4   = {exp:.10f}")
    print(f"  Difference:            {exact_exp - exp:.10f}")
    print(f"  Nearest half-integer:  {round(exact_exp * 2) / 2}")
    print()


def plot_heatmap(error_grid, outdir):
    """Plot error heatmap over (a, b) space."""
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    a_max, b_max = error_grid.shape

    # Clip errors for display
    err_display = np.clip(error_grid, 1e-6, 1.0)

    im = ax.imshow(err_display.T, origin='lower', aspect='auto',
                   norm=LogNorm(vmin=1e-5, vmax=1),
                   cmap='viridis_r',
                   extent=[0.5, a_max + 0.5, 0.5, b_max + 0.5])

    cbar = plt.colorbar(im, ax=ax, label='Relative Error')

    # Mark the best pair
    ax.plot(21, 15, 'o', markersize=15, markerfacecolor='none',
            markeredgecolor=CORAL, markeredgewidth=3,
            label=f'Best: (21, 15), error = {relative_error(21, 15):.4%}')

    # Mark runner-ups
    runners = [(21, 14), (21, 16)]
    for a, b in runners:
        ax.plot(a, b, 's', markersize=10, markerfacecolor='none',
                markeredgecolor=GOLD, markeredgewidth=2)

    ax.set_xlabel('$a$', fontsize=14)
    ax.set_ylabel('$b$', fontsize=14)
    ax.set_title('Electron Mass Formula: $m_e = m_P \\cdot \\alpha^{a/2 - b\\alpha/4}$\n'
                 'Relative Error over $(a, b)$ Parameter Space',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='upper right')

    plt.tight_layout()
    path = os.path.join(outdir, 'sim_mass_formula_heatmap.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_top_pairs(results, outdir, n=20):
    """Bar chart of top N pairs by accuracy."""
    fig, ax = plt.subplots(figsize=(14, 6))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    top = results[:n]
    labels = [f'({r[1]},{r[2]})' for r in top]
    errors = [r[0] * 100 for r in top]  # percent

    colors = []
    for r in top:
        if r[0] < 1e-4:
            colors.append(PURPLE)
        elif r[0] < 1e-3:
            colors.append(TEAL)
        elif r[0] < 1e-2:
            colors.append(GOLD)
        else:
            colors.append(CORAL)

    bars = ax.bar(range(n), errors, color=colors, alpha=0.85,
                  edgecolor='white', linewidth=1)

    # Value labels
    for i, (bar, err) in enumerate(zip(bars, errors)):
        va = 'bottom' if err > 0.05 else 'top'
        y = err + 0.005 if va == 'bottom' else err - 0.005
        ax.text(bar.get_x() + bar.get_width() / 2, y,
                f'{err:.3f}%', ha='center', va=va,
                fontsize=8, fontweight='bold',
                color='black' if err > 0.02 else 'white')

    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Relative Error (%)', fontsize=12)
    ax.set_xlabel('$(a, b)$ pair', fontsize=12)
    ax.set_title(f'Top {n} Integer Pairs $(a, b)$ for $m_e = m_P \\cdot \\alpha^{{a/2 - b\\alpha/4}}$',
                 fontsize=13, fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PURPLE, label='< 0.01%'),
        Patch(facecolor=TEAL, label='< 0.1%'),
        Patch(facecolor=GOLD, label='< 1%'),
        Patch(facecolor=CORAL, label='> 1%'),
    ]
    ax.legend(handles=legend_elements, title='Error Range',
              fontsize=9, title_fontsize=10, loc='upper right')

    plt.tight_layout()
    path = os.path.join(outdir, 'sim_mass_formula_top_pairs.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def main():
    """Main computation pipeline."""
    print("=" * 80)
    print("  Electron Mass Formula: Integer Pair Search")
    print("  m_e = m_P * alpha^(a/2 - b*alpha/4)")
    print("=" * 80)
    print()
    print(f"  Physical constants:")
    print(f"    m_e   = {M_E:.10e} kg")
    print(f"    m_P   = {M_P:.10e} kg")
    print(f"    alpha = {ALPHA:.12f} = 1/{1/ALPHA:.6f}")
    print(f"    m_e/m_P = {M_E/M_P:.10e}")
    print()

    # --- Standard search (a, b in [1, 50]) ---
    print("Searching (a, b) in [1, 50] x [1, 50]...")
    results, error_grid = search_pairs(50, 50)
    print_results(results, "STANDARD SEARCH: a, b in [1, 50]", n=30)
    print_analysis(results)

    # --- Extended search for sub-1% matches ---
    print("Extended search (a in [1, 100], b in [1, 200], error < 1%)...")
    ext_results = extended_search(100, 200)
    print_results(ext_results, "EXTENDED SEARCH: Sub-1% Matches", n=20)

    # --- Generate plots ---
    print("Generating figures...")
    os.makedirs(OUTDIR, exist_ok=True)
    plot_heatmap(error_grid, OUTDIR)
    plot_top_pairs(results, OUTDIR)

    print()
    print("=" * 80)
    print("  SEARCH COMPLETE")
    print("=" * 80)

    return results, error_grid


if __name__ == '__main__':
    results, error_grid = main()
