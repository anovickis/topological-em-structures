"""
Gamma-Ray Flux Predictions and Observational Limits
=====================================================

Computes the expected gamma-ray line flux from DM annihilation (Eq. 13.1 / 13.26):

    Phi = <sigma v> / (8 pi m^2) * J

using an NFW dark matter density profile:

    rho(r) = rho_s / [(r/r_s)(1 + r/r_s)^2]

The script:
  - Computes J-factors for different observation cone half-angles (1, 5, 10, 45 deg)
  - Plots Phi vs <sigma v> for different DM masses (0.5, 1, 2, 4 MeV)
  - Overlays COMPTEL, INTEGRAL/SPI sensitivity limits
  - Shows projected COSI and AMEGO-X sensitivity
  - Marks the framework's predicted cross-section range
  - Saves all plots

References: Paper Eqs. 13.1, 13.26-13.34
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ──────────────────────────────────────────────────────────────
# Color palette
# ──────────────────────────────────────────────────────────────
CORAL   = "#e76f51"
TEAL    = "#2a9d8f"
GOLD    = "#e9c46a"
PURPLE  = "#a855f7"

# ──────────────────────────────────────────────────────────────
# Physical constants & NFW parameters
# ──────────────────────────────────────────────────────────────
kpc_to_cm = 3.086e21       # cm per kpc
GeV_to_erg = 1.602e-3      # erg per GeV

# Milky Way NFW profile parameters
rho_s   = 0.184            # GeV/cm^3 (characteristic density for MW)
r_s     = 20.0             # kpc (scale radius)
r_sun   = 8.5              # kpc (solar distance from GC)
rho_local = 0.3            # GeV/cm^3 (local DM density, for normalization check)


def nfw_density(r_kpc):
    """NFW density profile rho(r) in GeV/cm^3."""
    x = r_kpc / r_s
    # Avoid division by zero
    x = np.maximum(x, 1e-6)
    return rho_s / (x * (1.0 + x)**2)


def compute_J_factor(theta_deg, r_sun_kpc=8.5, l_max_kpc=200.0, n_l=2000, n_theta=200):
    """
    Compute the J-factor for a cone of half-angle theta centered on the GC.

    J = integral_{Delta Omega} integral_{l.o.s.} rho^2(r) dl dOmega

    where r = sqrt(l^2 + r_sun^2 - 2 l r_sun cos(psi)) and psi is the angle
    from the GC direction.

    Uses direct numerical integration over the cone.

    Returns J in GeV^2 cm^{-5} (integrated over the solid angle).
    """
    theta_rad = np.deg2rad(theta_deg)
    Delta_Omega = 2.0 * np.pi * (1.0 - np.cos(theta_rad))  # sr

    # For small angles, integrate along the l.o.s. for several psi values
    # and average over the solid angle.
    # We use the simpler approach of averaging over psi.

    psi_arr = np.linspace(0, theta_rad, max(n_theta, 10))
    l_arr = np.linspace(0.01, l_max_kpc, n_l)
    dl = l_arr[1] - l_arr[0]

    J_total = 0.0

    for psi in psi_arr:
        # Distance from GC: r = sqrt(l^2 + r_sun^2 - 2*l*r_sun*cos(psi))
        r = np.sqrt(l_arr**2 + r_sun_kpc**2 - 2.0 * l_arr * r_sun_kpc * np.cos(psi))
        rho_sq = nfw_density(r)**2
        # Line-of-sight integral in GeV^2/cm^6 * kpc
        # Use np.trapezoid (NumPy >=2.0)
        _trapz = np.trapezoid
        los_integral = _trapz(rho_sq, l_arr)  # GeV^2/cm^6 * kpc
        # Convert kpc to cm
        los_integral *= kpc_to_cm  # GeV^2/cm^5

        # Weight by solid angle element: sin(psi) dpsi * 2pi
        # Integral over Omega = integral_0^theta sin(psi) dpsi * 2pi
        if psi > 0:
            J_total += los_integral * np.sin(psi) * 2.0 * np.pi

    # Multiply by dpsi
    dpsi = psi_arr[1] - psi_arr[0] if len(psi_arr) > 1 else theta_rad
    J_total *= dpsi  # GeV^2/cm^5

    return J_total, Delta_Omega


def gamma_ray_flux(sigma_v, m_GeV, J_factor):
    """
    Compute gamma-ray line flux from DM annihilation (Eq. 13.26):
        Phi = <sigma v> / (8 pi m^2) * J
    where sigma_v is in cm^3/s, m is in GeV, J in GeV^2/cm^5.

    Returns Phi in ph/cm^2/s.
    """
    return sigma_v / (8.0 * np.pi * m_GeV**2) * J_factor


# ──────────────────────────────────────────────────────────────
# Main computation
# ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("Gamma-Ray Flux Predictions")
    print("=" * 70)
    print()

    # ── Check NFW normalization ──
    rho_at_sun = nfw_density(r_sun)
    print(f"NFW profile check: rho(r_sun = {r_sun} kpc) = {rho_at_sun:.3f} GeV/cm^3")
    print(f"  (Target local density: {rho_local} GeV/cm^3)")
    # Adjust rho_s to match local density if needed
    rho_s_adjusted = rho_local / (nfw_density(r_sun) / rho_s)
    print(f"  Adjusted rho_s = {rho_s_adjusted:.3f} GeV/cm^3 to match local density")

    # Redefine with adjusted normalization
    rho_s_orig = rho_s
    rho_s = rho_s_adjusted

    # ── Compute J-factors for different cone angles ──
    theta_values = [1.0, 5.0, 10.0, 45.0]  # degrees
    print(f"\nJ-factors (NFW profile, r_s = {r_s} kpc, rho_s = {rho_s:.3f} GeV/cm^3):")
    print(f"{'Theta (deg)':>12s}  {'J (GeV^2/cm^5)':>18s}  {'Delta Omega (sr)':>18s}  {'log10(J)':>10s}")

    J_dict = {}
    for theta in theta_values:
        J, dOmega = compute_J_factor(theta)
        J_dict[theta] = J
        print(f"  {theta:8.1f}       {J:18.4e}       {dOmega:18.4e}       {np.log10(J):8.2f}")

    # Use the 1-degree J-factor as the reference
    J_ref = J_dict[1.0]

    # Cross-check with paper's Eq. 13.29:
    # J(1 deg) ~ 3e20 GeV^2/cm^5 (order of magnitude)
    print(f"\nPaper reference: J(1 deg) ~ 3e20 GeV^2/cm^5")
    print(f"Our calculation: J(1 deg) = {J_ref:.3e} GeV^2/cm^5")
    print(f"  (Factor of {J_ref / 3e20:.1f} from paper value -- depends on profile normalization)")

    # ──────────────────────────────────────────────────────────
    # Plot 1: Flux vs <sigma v> for different DM masses
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 7))

    sigma_v_range = np.logspace(-32, -24, 200)  # cm^3/s
    mass_values = [0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3]   # GeV
    mass_labels = ["0.5 MeV", "1 MeV", "2 MeV", "4 MeV"]
    colors_mass = [TEAL, CORAL, GOLD, PURPLE]

    for m_GeV, m_lab, clr in zip(mass_values, mass_labels, colors_mass):
        phi = gamma_ray_flux(sigma_v_range, m_GeV, J_ref)
        ax.loglog(sigma_v_range, phi, color=clr, linewidth=2.2, label=f"$m = {m_lab}$")

    # Experimental sensitivity limits
    # INTEGRAL/SPI: ~3e-5 ph/cm^2/s (narrow line, GC point source)
    ax.axhline(y=3e-5, color='gray', linewidth=2, linestyle='--', alpha=0.7)
    ax.text(1e-32, 5e-5, "INTEGRAL/SPI ($\\sim 3\\times 10^{-5}$)", fontsize=10,
            color='gray', fontweight='bold')

    # COMPTEL: ~1e-4 ph/cm^2/s/sr (diffuse)
    ax.axhline(y=1e-4, color='gray', linewidth=1.5, linestyle=':', alpha=0.6)
    ax.text(1e-32, 1.5e-4, "COMPTEL (diffuse)", fontsize=9, color='gray')

    # Projected COSI: ~1e-6 ph/cm^2/s
    ax.axhline(y=1e-6, color='blue', linewidth=1.5, linestyle='-.', alpha=0.5)
    ax.text(1e-32, 1.5e-6, "COSI (projected, 2027+)", fontsize=10,
            color='blue', fontstyle='italic')

    # Projected AMEGO-X: ~1e-7 ph/cm^2/s
    ax.axhline(y=1e-7, color='navy', linewidth=1.2, linestyle='-.', alpha=0.4)
    ax.text(1e-32, 1.5e-7, "AMEGO-X (projected)", fontsize=9,
            color='navy', fontstyle='italic')

    # Framework's predicted cross-section range (Section 10.4)
    ax.axvspan(1e-27, 1e-25, alpha=0.1, color=GOLD)
    ax.text(3e-26, 1e-12, "Framework's\npredicted range\n$\\langle\\sigma v\\rangle$",
            fontsize=10, color=GOLD, fontweight='bold', fontstyle='italic',
            ha='center', va='bottom')

    # Thermal relic
    ax.axvline(x=3e-26, color='black', linewidth=1.5, linestyle=':', alpha=0.5)
    ax.text(4e-26, 1e1, "Thermal\nrelic", fontsize=9, color='black', alpha=0.5)

    ax.set_xlabel("$\\langle\\sigma v\\rangle$ (cm$^3$/s)", fontsize=14)
    ax.set_ylabel("Photon flux $\\Phi_\\gamma$ (ph cm$^{-2}$ s$^{-1}$)", fontsize=14)
    ax.set_title("Gamma-Ray Line Flux from DM Annihilation (GC, $\\theta = 1^\\circ$, NFW)",
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1e-32, 1e-24)
    ax.set_ylim(1e-15, 1e2)
    ax.legend(fontsize=12, loc='lower right')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    outpath = "C:/Users/alexn/OneDrive/Documents/v2/Research/gamma_ray_flux_vs_sigmav.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {outpath}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Plot 2: Flux vs DM mass for different observation angles
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(11, 7))

    m_range = np.logspace(-4, -2, 100)  # GeV (0.1 MeV to 10 MeV)
    sigma_v_ref = 1e-27   # cm^3/s (conservative framework value)

    angle_colors = {1.0: CORAL, 5.0: TEAL, 10.0: GOLD, 45.0: PURPLE}

    for theta in theta_values:
        J = J_dict[theta]
        phi = gamma_ray_flux(sigma_v_ref, m_range, J)
        ax.loglog(m_range * 1e3, phi, color=angle_colors[theta], linewidth=2.2,
                  label=f"$\\theta = {theta:.0f}^\\circ$ ($J = {J:.1e}$)")

    # Experimental limits (as horizontal lines)
    ax.axhline(y=3e-5, color='gray', linewidth=2, linestyle='--', alpha=0.7,
               label="INTEGRAL/SPI limit")
    ax.axhline(y=1e-6, color='blue', linewidth=1.5, linestyle='-.', alpha=0.5,
               label="COSI projected")
    ax.axhline(y=1e-7, color='navy', linewidth=1.2, linestyle='-.', alpha=0.4,
               label="AMEGO-X projected")

    # Mark the framework's predicted mass range
    ax.axvspan(0.5, 4.0, alpha=0.08, color=GOLD)
    ax.text(1.5, 1e2, "Predicted\nmass range", fontsize=10, color=GOLD,
            fontweight='bold', fontstyle='italic', ha='center')

    ax.set_xlabel("Dark Matter Mass $m_{\\rm dark}$ (MeV)", fontsize=14)
    ax.set_ylabel("Photon flux $\\Phi_\\gamma$ (ph cm$^{-2}$ s$^{-1}$)", fontsize=14)
    ax.set_title(f"Gamma-Ray Flux vs Mass ($\\langle\\sigma v\\rangle = {sigma_v_ref:.0e}$ cm$^3$/s)",
                 fontsize=14, fontweight='bold')
    ax.set_xlim(0.1, 10)
    ax.set_ylim(1e-10, 1e4)
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    outpath2 = "C:/Users/alexn/OneDrive/Documents/v2/Research/gamma_ray_flux_vs_mass.png"
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath2}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Print predicted fluxes for key parameters
    # ──────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("Predicted Fluxes (theta = 1 deg, NFW)")
    print("=" * 70)
    print(f"{'m (MeV)':>10s}  {'<sv> (cm3/s)':>14s}  {'Phi (ph/cm2/s)':>16s}  {'Detectable?':>15s}")
    print("-" * 60)
    for m_MeV in [0.5, 1.0, 2.0, 4.0]:
        m_GeV = m_MeV * 1e-3
        for sv in [1e-27, 1e-29, 1e-30]:
            phi = gamma_ray_flux(sv, m_GeV, J_ref)
            detectable = "COSI" if phi > 1e-6 else ("AMEGO-X" if phi > 1e-7 else "Below limits")
            print(f"  {m_MeV:8.1f}    {sv:14.1e}    {phi:16.4e}    {detectable}")

    # Paper's estimates (Eqs. 13.32-13.33)
    print(f"\nPaper estimates (m=1 MeV, theta=1 deg):")
    print(f"  <sv> = 1e-27 cm3/s => Phi ~ 1.2e-2 ph/cm2/s (Eq. 13.32)")
    print(f"  <sv> = 1e-30 cm3/s => Phi ~ 1.2e-5 ph/cm2/s (Eq. 13.33)")
    phi_check = gamma_ray_flux(1e-27, 1e-3, J_ref)
    print(f"  Our calculation:     Phi = {phi_check:.3e} ph/cm2/s")

    print("\nDone.")
