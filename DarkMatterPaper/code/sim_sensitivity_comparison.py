"""
Sensitivity Comparison: MeV Dark Matter Annihilation
=====================================================

Generates a publication-quality comparison plot showing current experimental
upper limits, projected sensitivities, and the framework's predicted dark
matter mass and annihilation cross-section ranges.

Overlaid data:
  - COMPTEL upper limit (~10^{-26} cm^3/s, 1-10 MeV)
  - INTEGRAL/SPI upper limit (~3x10^{-29} cm^3/s, 0.5-8 MeV, NFW GC)
  - Projected COSI sensitivity (~10^{-30} cm^3/s, 0.2-5 MeV)
  - Projected AMEGO-X sensitivity (~10^{-31} cm^3/s, 0.2-10 MeV)
  - Thermal relic target (3x10^{-26} cm^3/s)
  - Framework prediction band (10^{-27} to 10^{-25} cm^3/s, 0.5-4 MeV)
  - Ropelength mass predictions (trefoil ~2.0, figure-8 ~2.3, cinquefoil ~2.5 MeV)

Note: Sensitivity values are approximate order-of-magnitude estimates for the
purposes of this comparison plot.

References: Paper Sections 10-11, Table 2
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patheffects as pe

# ──────────────────────────────────────────────────────────────
# Color palette (consistent with other paper figures)
# ──────────────────────────────────────────────────────────────
CORAL   = "#e76f51"
TEAL    = "#2a9d8f"
GOLD    = "#e9c46a"
PURPLE  = "#a855f7"
GRAY    = "#6b7280"

# Derived colors
TEAL_DARK  = "#1f7a6d"
GOLD_DARK  = "#c9a63a"

# ──────────────────────────────────────────────────────────────
# Experimental / projected limit data (approximate)
# Each entry: (mass_array_MeV, sigmav_array_cm3s, label, style)
# ──────────────────────────────────────────────────────────────

# COMPTEL upper limit: ~10^{-26} cm^3/s for DM -> 2gamma, 1-10 MeV
comptel_mass = np.array([1.0, 10.0])
comptel_sv   = np.array([1e-26, 1e-26])

# INTEGRAL/SPI upper limit: ~3x10^{-29} cm^3/s, narrow lines, 0.5-8 MeV (NFW GC)
integral_mass = np.array([0.5, 8.0])
integral_sv   = np.array([3e-29, 3e-29])

# Projected COSI sensitivity (2027+): ~10^{-30} cm^3/s, 0.2-5 MeV narrow lines
cosi_mass = np.array([0.2, 5.0])
cosi_sv   = np.array([1e-30, 1e-30])

# Projected AMEGO-X sensitivity: ~10^{-31} cm^3/s, 0.2-10 MeV
amego_mass = np.array([0.2, 10.0])
amego_sv   = np.array([1e-31, 1e-31])

# Thermal relic cross-section
thermal_relic_sv = 3e-26  # cm^3/s

# Framework prediction band
pred_mass_lo, pred_mass_hi = 0.5, 4.0   # MeV
pred_sv_lo, pred_sv_hi     = 1e-27, 1e-25  # cm^3/s

# Ropelength mass predictions (3/4-power scaling)
m_trefoil    = 2.0   # MeV
m_figure8    = 2.3   # MeV
m_cinquefoil = 2.5   # MeV


def main():
    print("=" * 70)
    print("Sensitivity Comparison Plot: MeV Dark Matter Annihilation")
    print("=" * 70)
    print()

    # ──────────────────────────────────────────────────────────
    # Create the figure
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 8))

    # --- Framework prediction band (draw first so it's behind everything) ---
    pred_rect = plt.Rectangle(
        (pred_mass_lo, pred_sv_lo),
        pred_mass_hi - pred_mass_lo,
        pred_sv_hi - pred_sv_lo,
        facecolor=CORAL, alpha=0.25, edgecolor=CORAL,
        linewidth=2, linestyle='-', zorder=2,
        label="Framework prediction"
    )
    ax.add_patch(pred_rect)

    # Label the prediction band
    ax.text(
        np.sqrt(pred_mass_lo * pred_mass_hi),
        np.sqrt(pred_sv_lo * pred_sv_hi),
        "Framework\nprediction\nband",
        fontsize=11, fontweight='bold', color=CORAL,
        ha='center', va='center', zorder=6,
        path_effects=[pe.withStroke(linewidth=3, foreground='white')]
    )

    # --- Current experimental upper limits ---
    # COMPTEL
    ax.plot(comptel_mass, comptel_sv, color=TEAL, linewidth=2.5,
            linestyle='-', solid_capstyle='round', zorder=4)
    ax.fill_between(comptel_mass, comptel_sv, 1e-22,
                    color=TEAL, alpha=0.07, zorder=1)
    ax.annotate(
        "COMPTEL\n(upper limit, $2\\gamma$)",
        xy=(4.0, 1e-26), xytext=(5.5, 5e-25),
        fontsize=10, color=TEAL_DARK, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=TEAL_DARK, lw=1.3),
        ha='center', zorder=7
    )

    # INTEGRAL/SPI
    ax.plot(integral_mass, integral_sv, color=TEAL, linewidth=2.5,
            linestyle='--', solid_capstyle='round', zorder=4)
    ax.fill_between(integral_mass, integral_sv, 1e-22,
                    color=TEAL, alpha=0.04, zorder=1)
    ax.annotate(
        "INTEGRAL/SPI\n(NFW, narrow lines)",
        xy=(2.0, 3e-29), xytext=(0.15, 1e-28),
        fontsize=10, color=TEAL_DARK, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=TEAL_DARK, lw=1.3),
        ha='center', zorder=7
    )

    # --- Projected sensitivities ---
    # COSI
    ax.plot(cosi_mass, cosi_sv, color=GOLD_DARK, linewidth=2.5,
            linestyle='-.', solid_capstyle='round', zorder=4)
    ax.annotate(
        "COSI (proj. 2027+)",
        xy=(3.5, 1e-30), xytext=(6.5, 3e-30),
        fontsize=10, color=GOLD_DARK, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=GOLD_DARK, lw=1.3),
        ha='center', zorder=7
    )

    # AMEGO-X
    ax.plot(amego_mass, amego_sv, color=GOLD_DARK, linewidth=2.5,
            linestyle=':', solid_capstyle='round', zorder=4)
    ax.annotate(
        "AMEGO-X (proj.)",
        xy=(5.0, 1e-31), xytext=(7.5, 3e-32),
        fontsize=10, color=GOLD_DARK, fontweight='bold',
        arrowprops=dict(arrowstyle='->', color=GOLD_DARK, lw=1.3),
        ha='center', zorder=7
    )

    # --- Thermal relic target ---
    ax.axhline(y=thermal_relic_sv, color=GRAY, linewidth=1.8,
               linestyle='--', zorder=3, alpha=0.8,
               label=r"Thermal relic: $3\times10^{-26}$ cm$^3$/s")
    ax.text(
        0.11, thermal_relic_sv * 1.5,
        r"Thermal relic $\langle\sigma v\rangle = 3\times10^{-26}$ cm$^3$/s",
        fontsize=9.5, color=GRAY, fontstyle='italic', va='bottom', zorder=7
    )

    # --- Ropelength mass predictions (vertical lines) ---
    knot_masses = [
        (m_trefoil,    r"$3_1$ (trefoil)"),
        (m_figure8,    r"$4_1$ (figure-8)"),
        (m_cinquefoil, r"$5_1$ (cinquefoil)"),
    ]

    # Stagger labels vertically so they don't overlap
    knot_y_positions = [2e-23, 5e-24, 1e-24]

    for i, (m_val, knot_label) in enumerate(knot_masses):
        ax.axvline(
            x=m_val, color=PURPLE, linewidth=1.8,
            linestyle=(0, (5, 3)), alpha=0.8, zorder=3
        )
        y_pos = knot_y_positions[i]
        ax.text(
            m_val * 1.06, y_pos, knot_label,
            fontsize=9.5, color=PURPLE, fontweight='bold',
            rotation=0, va='center', ha='left', zorder=7,
            path_effects=[pe.withStroke(linewidth=2.5, foreground='white')]
        )

    # Single legend entry for the knot lines
    ax.plot([], [], color=PURPLE, linewidth=1.8, linestyle=(0, (5, 3)),
            label="Ropelength mass predictions")

    # ──────────────────────────────────────────────────────────
    # Axes formatting
    # ──────────────────────────────────────────────────────────
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.1, 10.0)
    ax.set_ylim(1e-32, 1e-22)

    ax.set_xlabel(r"Dark matter mass $m_\chi$ (MeV)", fontsize=14)
    ax.set_ylabel(r"$\langle\sigma v\rangle$ (cm$^3$ s$^{-1}$)", fontsize=14)
    ax.set_title(
        "MeV Dark Matter: Sensitivity Comparison and Framework Predictions",
        fontsize=14, fontweight='bold', pad=12
    )

    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=10)
    ax.grid(True, which='major', alpha=0.25, linewidth=0.8)
    ax.grid(True, which='minor', alpha=0.10, linewidth=0.5)

    # Legend -- current limits and projected
    # Build a custom legend with grouped entries
    from matplotlib.lines import Line2D

    legend_elements = [
        Line2D([0], [0], color=TEAL, linewidth=2.5, linestyle='-',
               label="COMPTEL (current)"),
        Line2D([0], [0], color=TEAL, linewidth=2.5, linestyle='--',
               label="INTEGRAL/SPI (current)"),
        Line2D([0], [0], color=GOLD_DARK, linewidth=2.5, linestyle='-.',
               label="COSI (projected)"),
        Line2D([0], [0], color=GOLD_DARK, linewidth=2.5, linestyle=':',
               label="AMEGO-X (projected)"),
        Line2D([0], [0], color=GRAY, linewidth=1.8, linestyle='--',
               label=r"Thermal relic ($3\times10^{-26}$)"),
        plt.Rectangle((0, 0), 1, 1, facecolor=CORAL, alpha=0.25,
                       edgecolor=CORAL, linewidth=1.5,
                       label="Framework prediction band"),
        Line2D([0], [0], color=PURPLE, linewidth=1.8, linestyle=(0, (5, 3)),
               label="Ropelength mass predictions"),
    ]
    ax.legend(
        handles=legend_elements, loc='lower left',
        fontsize=9.5, framealpha=0.92, edgecolor='#cccccc',
        borderpad=0.8, labelspacing=0.5
    )

    # ──────────────────────────────────────────────────────────
    # Annotation: approximate limits disclaimer
    # ──────────────────────────────────────────────────────────
    ax.text(
        0.98, 0.02,
        "Limits are approximate; see text for details.",
        transform=ax.transAxes, fontsize=8.5, color='gray',
        ha='right', va='bottom', fontstyle='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor='#dddddd', alpha=0.9)
    )

    # ──────────────────────────────────────────────────────────
    # Arrow annotations showing exclusion direction
    # ──────────────────────────────────────────────────────────
    # "Excluded" label above the COMPTEL line
    ax.annotate(
        "", xy=(7.5, 3e-26), xytext=(7.5, 3e-24),
        arrowprops=dict(arrowstyle='<-', color=TEAL, lw=1.5, alpha=0.5),
        zorder=3
    )
    ax.text(
        7.8, 3e-25, "excluded", fontsize=9, color=TEAL,
        rotation=90, va='center', alpha=0.6
    )

    # ──────────────────────────────────────────────────────────
    # Save figure
    # ──────────────────────────────────────────────────────────
    plt.tight_layout()
    outpath = "C:/Users/alexn/OneDrive/Documents/v2/Research/sensitivity_comparison.png"
    plt.savefig(outpath, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"Plot saved: {outpath}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Print summary table
    # ──────────────────────────────────────────────────────────
    print()
    print("Summary of plotted limits and predictions:")
    print("-" * 70)
    print(f"  {'Curve':<30s} {'<sigma v> (cm^3/s)':<22s} {'Mass range (MeV)'}")
    print("-" * 70)
    print(f"  {'COMPTEL upper limit':<30s} {'~1e-26':<22s} {'1 - 10'}")
    print(f"  {'INTEGRAL/SPI upper limit':<30s} {'~3e-29':<22s} {'0.5 - 8'}")
    print(f"  {'COSI (projected, 2027+)':<30s} {'~1e-30':<22s} {'0.2 - 5'}")
    print(f"  {'AMEGO-X (projected)':<30s} {'~1e-31':<22s} {'0.2 - 10'}")
    print(f"  {'Thermal relic':<30s} {'3e-26':<22s} {'(all)'}")
    print(f"  {'Framework prediction':<30s} {'1e-27 to 1e-25':<22s} {'0.5 - 4'}")
    print("-" * 70)
    print(f"  Ropelength mass predictions:")
    print(f"    Trefoil (3_1):     {m_trefoil} MeV")
    print(f"    Figure-8 (4_1):    {m_figure8} MeV")
    print(f"    Cinquefoil (5_1):  {m_cinquefoil} MeV")
    print()
    print("Done.")


if __name__ == "__main__":
    main()
