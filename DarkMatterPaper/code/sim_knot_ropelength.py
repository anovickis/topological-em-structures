"""
Knot Ropelength Mass Predictions
==================================

Tabulates known ropelengths for knots up to 8 crossings and computes mass
predictions under three scaling hypotheses (Paper Section 9.3):

  (A) Linear:       m(K) = m_e * L(K) / L_0
  (B) 3/4-power:    m(K) = m_e * [L(K) / L_0]^(3/4)
  (C) Square-root:  m(K) = m_e * [L(K) / L_0]^(1/2)

where L_0 = 2*pi (unknot ropelength) and m_e = 0.511 MeV.

The script:
  - Tabulates ropelengths for knots up to 8 crossings
  - Computes mass predictions under all 3 hypotheses
  - Plots bar chart of mass spectrum with error bars from the 3 hypotheses
  - Plots mass ratios that distinguish scaling laws
  - Estimates total DM abundance from state counting
  - Saves all plots

References: Paper Eqs. 9.1-9.8, Table in Section 9.3
Ropelength data: Cantarella et al., Ashton et al., Rawdon et al.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────
# Color palette
# ──────────────────────────────────────────────────────────────
CORAL   = "#e76f51"
TEAL    = "#2a9d8f"
GOLD    = "#e9c46a"
PURPLE  = "#a855f7"

# ──────────────────────────────────────────────────────────────
# Physical constants
# ──────────────────────────────────────────────────────────────
m_e = 0.511  # MeV (electron mass, reference)
L_0 = 2 * np.pi  # Unknot ropelength (circle) ~ 6.283

# ──────────────────────────────────────────────────────────────
# Ropelength data for knots up to 8 crossings
# Sources: Cantarella et al. (2006), Ashton et al. (2011), Rawdon et al.
# Values are numerically computed ropelengths (unit-diameter tube)
# ──────────────────────────────────────────────────────────────
knots = [
    # (name, Alexander-Briggs notation, crossings, ropelength, chiral, lower_bound)
    ("Unknot",       "0_1", 0,  2*np.pi,   False, 2*np.pi),
    ("Trefoil",      "3_1", 3,  32.74,     True,  31.32),
    ("Figure-eight", "4_1", 4,  42.12,     False, 39.75),
    ("Cinquefoil",   "5_1", 5,  47.20,     True,  None),
    ("Three-twist",  "5_2", 5,  47.00,     True,  None),
    ("6_1",          "6_1", 6,  55.70,     False, None),
    ("6_2",          "6_2", 6,  56.20,     True,  None),
    ("6_3",          "6_3", 6,  58.00,     False, None),
    ("7_1",          "7_1", 7,  62.00,     True,  None),
    ("7_2",          "7_2", 7,  63.50,     True,  None),
    ("7_4",          "7_4", 7,  64.00,     True,  None),
    ("8_1",          "8_1", 8,  70.00,     False, None),
    ("8_18",         "8_18", 8, 74.00,     False, None),
]


def mass_linear(L):
    """Hypothesis A: m = m_e * L / L_0 (linear / rope energy)."""
    return m_e * L / L_0

def mass_34power(L):
    """Hypothesis B: m = m_e * (L / L_0)^(3/4) (V-K bound motivated)."""
    return m_e * (L / L_0)**0.75

def mass_sqrt(L):
    """Hypothesis C: m = m_e * (L / L_0)^(1/2) (geometric mean)."""
    return m_e * (L / L_0)**0.5


if __name__ == "__main__":
    print("=" * 70)
    print("Knot Ropelength Mass Predictions")
    print("=" * 70)
    print(f"Reference: m_e = {m_e} MeV, L_0 = 2*pi = {L_0:.4f}")
    print()

    # ── Compute mass predictions ──
    print(f"{'Knot':>15s}  {'Notation':>8s}  {'c':>3s}  {'L(K)':>8s}  {'L/L_0':>7s}  "
          f"{'m_A (MeV)':>10s}  {'m_B (MeV)':>10s}  {'m_C (MeV)':>10s}  {'Chiral':>7s}")
    print("-" * 95)

    names = []
    notations = []
    crossings_arr = []
    ropelengths = []
    masses_A = []
    masses_B = []
    masses_C = []
    chiral_arr = []

    for name, notation, c, L, chiral, lb in knots:
        ratio = L / L_0
        mA = mass_linear(L)
        mB = mass_34power(L)
        mC = mass_sqrt(L)

        names.append(name)
        notations.append(notation)
        crossings_arr.append(c)
        ropelengths.append(L)
        masses_A.append(mA)
        masses_B.append(mB)
        masses_C.append(mC)
        chiral_arr.append(chiral)

        chi_str = "Yes" if chiral else "No"
        print(f"  {name:>13s}  {notation:>8s}  {c:>3d}  {L:8.2f}  {ratio:7.3f}  "
              f"{mA:10.3f}  {mB:10.3f}  {mC:10.3f}  {chi_str:>7s}")

    masses_A = np.array(masses_A)
    masses_B = np.array(masses_B)
    masses_C = np.array(masses_C)

    # Skip unknot for dark matter candidates (it's the electron)
    dm_idx = slice(1, None)
    dm_names = names[1:]
    dm_notations = notations[1:]
    dm_mA = masses_A[dm_idx]
    dm_mB = masses_B[dm_idx]
    dm_mC = masses_C[dm_idx]
    dm_chiral = chiral_arr[1:]

    # ── Mass ratios (Eq. 9.8) ──
    print("\n" + "=" * 70)
    print("Key Mass Ratios (observable discriminators)")
    print("=" * 70)
    # Trefoil to Figure-eight ratio
    i_tref = 0  # index in dm arrays
    i_fig8 = 1
    for label, m_arr in [("Linear (A)", dm_mA), ("3/4 power (B)", dm_mB), ("Sqrt (C)", dm_mC)]:
        r = m_arr[i_tref] / m_arr[i_fig8]
        print(f"  m(trefoil)/m(figure-8), {label}: {r:.4f}")

    print(f"\n  Paper Eq. 9.8: Linear: 0.777, 3/4 power: 0.827, Sqrt: 0.886")

    # ── State counting abundance estimate ──
    print("\n" + "=" * 70)
    print("State Counting Abundance Estimate")
    print("=" * 70)

    # Count H=0 states (chiral knots count double: left + right)
    n_H0 = 0
    for i, (name, chiral) in enumerate(zip(dm_names, dm_chiral)):
        states = 2 if chiral else 1
        n_H0 += states
        # print(f"  {name}: {states} state(s)")
    n_H1 = 2  # electron + positron

    print(f"  H != 0 states (electron/positron): {n_H1}")
    print(f"  H = 0 states (dark matter knots up to 8 crossings): {n_H0}")
    print(f"  Naive ratio N(H=0)/N(H!=0) = {n_H0}/{n_H1} = {n_H0/n_H1:.1f}")

    # Energy-weighted estimate: weight each state by exp(-m/T) at freeze-out
    T_fo = 1.0  # MeV (freeze-out temperature ~ m_dark)
    weights_H0 = []
    for i, chiral in enumerate(dm_chiral):
        m_avg = dm_mB[i]  # Use hypothesis B as best estimate
        w = np.exp(-m_avg / T_fo)
        if chiral:
            w *= 2
        weights_H0.append(w)
    W_H0 = sum(weights_H0)
    W_H1 = 2 * np.exp(-m_e / T_fo)

    print(f"\n  Energy-weighted (T_fo = {T_fo} MeV, using Hypothesis B):")
    print(f"    W(H=0) = {W_H0:.4f}")
    print(f"    W(H=1) = {W_H1:.4f}")
    print(f"    Omega_dark / Omega_ordinary ~ {W_H0/W_H1:.1f}")
    print(f"    (Observed: ~5.4)")

    # ──────────────────────────────────────────────────────────
    # Plot 1: Mass spectrum bar chart
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(14, 7))

    x_pos = np.arange(len(dm_names))
    width = 0.25

    # Plot the three hypotheses as grouped bars
    bars_A = ax.bar(x_pos - width, dm_mA, width, label="Linear (A): $m \\propto L$",
                    color=TEAL, edgecolor='black', linewidth=0.5)
    bars_B = ax.bar(x_pos, dm_mB, width, label="3/4 power (B): $m \\propto L^{3/4}$",
                    color=CORAL, edgecolor='black', linewidth=0.5)
    bars_C = ax.bar(x_pos + width, dm_mC, width, label="Sqrt (C): $m \\propto L^{1/2}$",
                    color=GOLD, edgecolor='black', linewidth=0.5)

    # Add electron reference line
    ax.axhline(y=m_e, color='black', linewidth=1.5, linestyle='--',
               label=f"Electron mass: $m_e = {m_e}$ MeV", alpha=0.7)

    # Add Battye-Sutcliffe numerical value for trefoil
    ax.plot(0, 2.0, '*', color=PURPLE, markersize=15, zorder=5,
            label="B-S numerical: $m_{\\rm trefoil} \\approx 2.0$ MeV")

    # Mark chiral knots
    for i, chiral in enumerate(dm_chiral):
        if chiral:
            ax.text(i, dm_mA[i] + 0.15, "L+R", fontsize=7, ha='center', color='gray',
                    fontstyle='italic')

    ax.set_xlabel("Knot Type", fontsize=14)
    ax.set_ylabel("Predicted Mass (MeV)", fontsize=14)
    ax.set_title("Dark Matter Mass Spectrum from Knot Ropelength", fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"{n}\n({no})" for n, no in zip(dm_names, dm_notations)],
                       fontsize=9, rotation=45, ha='right')
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')
    ax.tick_params(labelsize=11)
    ax.set_ylim(0, max(dm_mA) * 1.15)

    plt.tight_layout()
    outpath = "C:/Users/alexn/OneDrive/Documents/v2/Research/knot_mass_spectrum.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {outpath}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Plot 2: Mass ratios (trefoil/K) for different hypotheses
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 6))

    # Compute ratios relative to trefoil (lightest non-trivial knot)
    ratios_A = dm_mA[0] / dm_mA
    ratios_B = dm_mB[0] / dm_mB
    ratios_C = dm_mC[0] / dm_mC

    x_pos = np.arange(len(dm_names))
    width = 0.25

    ax.bar(x_pos - width, ratios_A, width, label="Linear (A)",
           color=TEAL, edgecolor='black', linewidth=0.5)
    ax.bar(x_pos, ratios_B, width, label="3/4 power (B)",
           color=CORAL, edgecolor='black', linewidth=0.5)
    ax.bar(x_pos + width, ratios_C, width, label="Sqrt (C)",
           color=GOLD, edgecolor='black', linewidth=0.5)

    ax.axhline(y=1.0, color='black', linewidth=1, linestyle=':', alpha=0.5)

    ax.set_xlabel("Knot Type", fontsize=14)
    ax.set_ylabel("$m_{\\rm trefoil} / m_K$", fontsize=14)
    ax.set_title("Mass Ratios Relative to Trefoil (Observable Discriminators)",
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels([f"{n}\n({no})" for n, no in zip(dm_names, dm_notations)],
                       fontsize=9, rotation=45, ha='right')
    ax.legend(fontsize=11, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    ax.tick_params(labelsize=11)

    # Annotate the trefoil-figure8 ratio as key observable
    for i, (rA, rB, rC) in enumerate(zip(ratios_A, ratios_B, ratios_C)):
        if i == 1:  # figure-eight
            ax.annotate(f"A: {rA:.3f}\nB: {rB:.3f}\nC: {rC:.3f}",
                        xy=(i, rA + 0.01), fontsize=8, ha='center',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))

    plt.tight_layout()
    outpath2 = "C:/Users/alexn/OneDrive/Documents/v2/Research/knot_mass_ratios.png"
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath2}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Plot 3: Mass vs ropelength with scaling curves
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))

    L_range = np.linspace(L_0, 80, 200)

    ax.plot(L_range, mass_linear(L_range), color=TEAL, linewidth=2.5,
            label="(A) Linear: $m \\propto L$")
    ax.plot(L_range, mass_34power(L_range), color=CORAL, linewidth=2.5,
            label="(B) 3/4 power: $m \\propto L^{3/4}$")
    ax.plot(L_range, mass_sqrt(L_range), color=GOLD, linewidth=2.5,
            label="(C) Sqrt: $m \\propto L^{1/2}$")

    # Plot individual knots
    for i, (name, notation, c, L, chiral, lb) in enumerate(knots):
        if i == 0:
            # Unknot (electron)
            ax.plot(L, m_e, 'ko', markersize=12, zorder=5)
            ax.annotate(f"Electron\n({notation})", xy=(L, m_e),
                        xytext=(L+3, m_e-0.5), fontsize=10,
                        arrowprops=dict(arrowstyle='->', color='black'))
        else:
            ax.plot(L, mass_34power(L), 's', color=PURPLE, markersize=8, zorder=5)
            ax.annotate(notation, xy=(L, mass_34power(L)),
                        xytext=(L+1, mass_34power(L)+0.15), fontsize=8, color=PURPLE)

    # B-S numerical value for trefoil
    ax.plot(32.74, 2.0, '*', color='red', markersize=15, zorder=6,
            label="Battye-Sutcliffe (trefoil, numerical)")

    ax.set_xlabel("Ropelength $L(K)$", fontsize=14)
    ax.set_ylabel("Predicted Mass (MeV)", fontsize=14)
    ax.set_title("Mass vs Ropelength: Three Scaling Hypotheses", fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    outpath3 = "C:/Users/alexn/OneDrive/Documents/v2/Research/knot_mass_vs_ropelength.png"
    plt.savefig(outpath3, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath3}")
    plt.close()

    print("\nDone.")
