"""
Polarizability Cross-Section and Detection Prospects
=====================================================

Computes the Rayleigh scattering cross-section for a neutral, compact topological
dark matter particle of radius a (Eqs. 10.2-10.4 of the paper):

    sigma_pol = (128 pi^5 / 3) * a^6 / lambda^4

and the geometric cross-section sigma_geom = pi a^2.

The script:
  - Plots sigma_pol vs photon energy from 1 eV to 10 MeV
  - Marks the transition to geometric cross-section
  - Overlays experimental sensitivity limits (XENON, LZ)
  - Shows the Bullet Cluster bound sigma/m < 1 cm^2/g
  - Plots the annihilation cross-section <sigma v> vs f_top for different masses
  - Saves all plots

References: Paper Eqs. 10.1-10.11, 10.30-10.37
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
# Physical constants (CGS)
# ──────────────────────────────────────────────────────────────
hbar    = 1.0546e-27      # erg s
c       = 2.998e10        # cm/s
m_e     = 9.109e-28       # g (electron mass)
m_e_MeV = 0.511           # MeV
eV_to_erg = 1.602e-12     # erg per eV
alpha   = 1.0 / 137.036   # fine structure constant

# Compton wavelength of the electron (natural particle radius)
lambda_C = hbar / (m_e * c)   # ~ 3.86e-11 cm
a = lambda_C                    # particle radius ~ Compton wavelength

print(f"Particle radius a = lambda_C = {a:.3e} cm")
print(f"Geometric cross-section sigma_geom = pi a^2 = {np.pi * a**2:.3e} cm^2")

# ──────────────────────────────────────────────────────────────
# 1. Polarizability cross-section vs photon energy
# ──────────────────────────────────────────────────────────────
def sigma_rayleigh(wavelength_cm, radius_cm):
    """Rayleigh scattering cross-section: (128 pi^5 / 3) a^6 / lambda^4"""
    return (128.0 * np.pi**5 / 3.0) * radius_cm**6 / wavelength_cm**4

def sigma_geometric(radius_cm):
    """Geometric cross-section: pi a^2"""
    return np.pi * radius_cm**2

# Photon energy range: 1 eV to 10 MeV
E_eV = np.logspace(0, 7, 1000)   # eV
E_erg = E_eV * eV_to_erg
wavelength = hbar * c * 2 * np.pi / E_erg   # cm (lambda = 2 pi hbar c / E)

sigma_pol = sigma_rayleigh(wavelength, a)
sigma_geo = sigma_geometric(a)

# The Rayleigh formula is valid only for lambda >> a, i.e., a/lambda << 1
a_over_lambda = a / wavelength

# Effective cross-section: min(sigma_pol, sigma_geo) for the physical cross-section
sigma_eff = np.minimum(sigma_pol, sigma_geo)


if __name__ == "__main__":
    print("=" * 70)
    print("Polarizability Cross-Section Calculation")
    print("=" * 70)
    print()

    # Print key values from the paper
    energies_check = [2.48, 1.24e4, 1.0e6]   # eV: optical (500nm), X-ray (0.1nm), MeV gamma
    labels_check = ["Optical (500 nm)", "X-ray (0.1 nm)", "MeV gamma (1.24 pm)"]
    for Ech, lab in zip(energies_check, labels_check):
        lam = hbar * c * 2 * np.pi / (Ech * eV_to_erg)
        sig = sigma_rayleigh(lam, a)
        ratio = a / lam
        print(f"  {lab:25s}: lambda = {lam:.3e} cm, sigma_pol = {sig:.3e} cm^2, a/lambda = {ratio:.3e}")

    print(f"\n  Geometric limit: sigma_geom = {sigma_geo:.3e} cm^2")
    print()

    # ──────────────────────────────────────────────────────────
    # Plot 1: Cross-section vs photon energy
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 7))

    # Rayleigh cross-section
    ax.loglog(E_eV, sigma_pol, color=TEAL, linewidth=2.5,
              label=r"Rayleigh: $\sigma_{\rm pol} = \frac{128\pi^5}{3}\frac{a^6}{\lambda^4}$")

    # Geometric cross-section
    ax.axhline(y=sigma_geo, color=CORAL, linewidth=2, linestyle='--',
               label=r"Geometric: $\sigma_{\rm geom} = \pi a^2$")

    # Effective (physical) cross-section
    ax.loglog(E_eV, sigma_eff, color='black', linewidth=1.5, linestyle=':',
              label="Effective (physical)", alpha=0.6)

    # Mark the transition energy where sigma_pol = sigma_geo
    # (128 pi^5 / 3) a^6 / lambda^4 = pi a^2
    # lambda_trans^4 = (128 pi^4 / 3) a^4
    # lambda_trans = (128 pi^4 / 3)^(1/4) * a
    lambda_trans = (128.0 * np.pi**4 / 3.0)**0.25 * a
    E_trans_erg = hbar * c * 2 * np.pi / lambda_trans
    E_trans_eV = E_trans_erg / eV_to_erg
    ax.axvline(x=E_trans_eV, color=GOLD, linewidth=1.5, linestyle='-.',
               label=f"Transition: $E \\approx {E_trans_eV/1e6:.2f}$ MeV")

    # Annotate specific photon types
    markers = [
        (2.48, "Optical\n(500 nm)", TEAL),
        (1.24e4, "X-ray\n(0.1 nm)", PURPLE),
        (1e6, "MeV\ngamma", CORAL),
    ]
    for E_mark, label, clr in markers:
        lam_m = hbar * c * 2 * np.pi / (E_mark * eV_to_erg)
        sig_m = sigma_rayleigh(lam_m, a)
        ax.plot(E_mark, min(sig_m, sigma_geo), 'o', color=clr, markersize=10, zorder=5)
        yoff = 3 if sig_m < sigma_geo else 0.1
        ax.annotate(label, xy=(E_mark, min(sig_m, sigma_geo)),
                    xytext=(E_mark*2, min(sig_m, sigma_geo)*yoff),
                    fontsize=10, color=clr, fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color=clr, lw=1.2))

    # Experimental bounds (approximate, for reference)
    # XENON1T: sigma ~ 10^{-46} cm^2 for m ~ 30 GeV (nuclear recoil)
    # LZ: sigma ~ 10^{-48} cm^2 for m ~ 40 GeV
    # These are for nuclear recoil, not polarizability scattering, but shown for context
    ax.axhline(y=1e-46, color='gray', linewidth=1, linestyle='--', alpha=0.5)
    ax.text(2, 3e-46, "XENON1T limit\n($m \\sim 30$ GeV, nuclear recoil)", fontsize=8,
            color='gray', alpha=0.7)
    ax.axhline(y=1e-48, color='gray', linewidth=1, linestyle='--', alpha=0.5)
    ax.text(2, 3e-48, "LZ limit\n($m \\sim 40$ GeV, nuclear recoil)", fontsize=8,
            color='gray', alpha=0.7)

    # Bullet Cluster bound: sigma/m < 1 cm^2/g for self-interaction
    # For m = 1 MeV = 1.78e-27 g: sigma < 1.78e-27 cm^2
    m_dark_g = 1.0e-3 * 1.783e-24  # 1 MeV in grams
    sigma_bullet = 1.0 * m_dark_g   # cm^2
    ax.axhline(y=sigma_bullet, color='red', linewidth=1.5, linestyle='-.',
               alpha=0.6)
    ax.text(2, sigma_bullet * 3,
            f"Bullet Cluster: $\\sigma/m < 1$ cm$^2$/g ($\\sigma < {sigma_bullet:.1e}$ cm$^2$)",
            fontsize=9, color='red', alpha=0.7)

    ax.set_xlabel("Photon Energy (eV)", fontsize=14)
    ax.set_ylabel("Cross-section (cm$^2$)", fontsize=14)
    ax.set_title("Polarizability Cross-Section of Topological Dark Matter ($a = \\lambda_C$)",
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1, 1e7)
    ax.set_ylim(1e-52, 1e-12)
    ax.legend(fontsize=10, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    # Shaded regions
    ax.fill_between([1, E_trans_eV], 1e-52, 1e-12, alpha=0.03, color=TEAL)
    ax.fill_between([E_trans_eV, 1e7], 1e-52, 1e-12, alpha=0.03, color=CORAL)
    ax.text(30, 1e-13, "Rayleigh regime\n($\\lambda \\gg a$)", fontsize=10,
            color=TEAL, ha='center', fontstyle='italic')
    ax.text(5e6, 1e-13, "Geometric\nregime", fontsize=10,
            color=CORAL, ha='center', fontstyle='italic')

    plt.tight_layout()
    outpath = "C:/Users/alexn/OneDrive/Documents/v2/Research/cross_section_polarizability.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Plot 2: Annihilation cross-section <sigma v> vs f_top
    # ──────────────────────────────────────────────────────────
    # From Eq. 10.5: <sigma v> = pi alpha^2 / m^2 * f_top * v
    # With v ~ 10^{-3} c for galactic DM

    fig, ax = plt.subplots(figsize=(10, 7))

    f_top = np.logspace(-10, 0, 200)
    v_gal = 1e-3 * c  # cm/s

    mass_values = [0.5, 1.0, 2.0, 4.0]  # MeV
    colors_mass = [TEAL, CORAL, GOLD, PURPLE]

    for m_val, clr in zip(mass_values, colors_mass):
        m_GeV = m_val * 1e-3
        # sigma = pi alpha^2 / m^2 in natural units -> convert to cm^2
        # m in GeV: 1/m^2 in GeV^{-2}, convert to cm^2: * (hbar c)^2 = 3.894e-28 cm^2 GeV^2
        hbar_c_sq = 3.894e-28  # cm^2 GeV^2
        sigma_base = np.pi * alpha**2 * hbar_c_sq / m_GeV**2  # cm^2
        sigma_v = sigma_base * v_gal * f_top  # cm^3/s

        ax.loglog(f_top, sigma_v, color=clr, linewidth=2,
                  label=f"$m = {m_val}$ MeV")

    # Thermal relic target
    ax.axhline(y=3e-26, color='black', linewidth=2, linestyle='--',
               label="Thermal relic: $3\\times 10^{-26}$ cm$^3$/s")

    # Viable range from the paper (Section 10.4)
    ax.axhspan(1e-27, 1e-25, alpha=0.1, color=GOLD)
    ax.text(1e-8, 5e-26, "Astrophysically\nviable range", fontsize=10,
            color=GOLD, fontweight='bold', fontstyle='italic', ha='center')

    # Instanton suppression
    ax.axhline(y=1e-40, color='gray', linewidth=1, linestyle=':', alpha=0.5)
    ax.text(1e-1, 3e-40, "Effectively stable\n(instanton suppressed)", fontsize=9,
            color='gray', ha='right')

    # Gamma-ray limits
    ax.axhline(y=1e-29, color='red', linewidth=1, linestyle='-.', alpha=0.6)
    ax.text(1e-9, 3e-29, "INTEGRAL/SPI limit (GC, 1 MeV)", fontsize=9,
            color='red', alpha=0.7)

    ax.set_xlabel("Topological suppression factor $f_{\\rm top}$", fontsize=14)
    ax.set_ylabel("$\\langle\\sigma v\\rangle$ (cm$^3$/s)", fontsize=14)
    ax.set_title("Annihilation Cross-Section vs Topological Suppression",
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1e-10, 1)
    ax.set_ylim(1e-42, 1e-14)
    ax.legend(fontsize=11, loc='upper left')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    outpath2 = "C:/Users/alexn/OneDrive/Documents/v2/Research/cross_section_annihilation.png"
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath2}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Print summary of self-interaction vs Bullet Cluster
    # ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("Self-Interaction vs Bullet Cluster Constraint")
    print("=" * 70)
    print(f"Bullet Cluster bound: sigma/m < 1 cm^2/g")
    print(f"For m_dark = 1 MeV = {m_dark_g:.3e} g:")
    print(f"  => sigma < {sigma_bullet:.3e} cm^2")
    print()

    # Contact cross-section (Eq. 10.34)
    r_e = alpha * lambda_C
    sigma_contact = np.pi * r_e**2
    print(f"Contact cross-section (Eq. 10.34): sigma = pi r_e^2 = {sigma_contact:.3e} cm^2")
    print(f"  sigma/m = {sigma_contact / m_dark_g:.3e} cm^2/g")
    print(f"  Ratio to Bullet limit: {sigma_contact / sigma_bullet:.2f}")
    print()

    # Required topological suppression
    f_top_req = sigma_bullet / sigma_contact
    print(f"Required f_top for Bullet Cluster consistency: f_top < {f_top_req:.2f}")
    print(f"(Paper Eq. 10.37: f_top < 0.27 -- moderate suppression)")
    print()

    # Geometric cross-section
    sigma_compton = np.pi * lambda_C**2
    print(f"Compton geometric cross-section: sigma = pi lambda_C^2 = {sigma_compton:.3e} cm^2")
    print(f"  sigma/m = {sigma_compton / m_dark_g:.3e} cm^2/g")
    print(f"  Exceeds Bullet limit by factor: {sigma_compton / sigma_bullet:.0e}")
    print("  (This is why geometric cross-section is the wrong estimate!)")

    print("\nDone.")
