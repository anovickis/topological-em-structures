"""
Boltzmann Freeze-out Calculation for Topological Dark Matter
=============================================================

Solves the Boltzmann equation for the comoving number density Y = n/s
(Eq. 11.1 of the paper):

    dY/dx = -(s <sigma v>) / (H x) * (Y^2 - Y_eq^2)

where x = m/T is the dimensionless inverse temperature.

The script:
  - Solves the ODE numerically using scipy.integrate.solve_ivp
  - Uses m_dark = 1 MeV with various <sigma v> from 10^-28 to 10^-24 cm^3/s
  - Computes Omega h^2 for each cross-section
  - Plots Y(x) vs x for different cross-sections
  - Plots Omega h^2 vs <sigma v> with the observed value Omega h^2 = 0.120
  - Marks the thermal relic ("WIMP miracle") cross-section
  - Saves all plots

References: Paper Eqs. 11.1-11.6, Kolb & Turner "The Early Universe"
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.special import kn  # Modified Bessel functions

# ──────────────────────────────────────────────────────────────
# Color palette
# ──────────────────────────────────────────────────────────────
CORAL   = "#e76f51"
TEAL    = "#2a9d8f"
GOLD    = "#e9c46a"
PURPLE  = "#a855f7"

# ──────────────────────────────────────────────────────────────
# Physical constants (CGS / natural units)
# ──────────────────────────────────────────────────────────────
M_Pl       = 1.22e19       # Planck mass in GeV
m_dark_GeV = 1.0e-3        # Dark matter mass in GeV (= 1 MeV)
g_star     = 10.75          # Effective relativistic d.o.f. at ~ 1 MeV (e+-, nu, gamma)
g_star_s   = 10.75          # Entropic d.o.f.
g_dm       = 1              # Internal d.o.f. of the dark matter particle (scalar-like)

# Observed dark matter relic density
Omega_obs_h2 = 0.120

# Entropy density prefactor: s = (2 pi^2 / 45) g_*s T^3
# Hubble: H = sqrt(pi^2 g_* / 90) T^2 / M_Pl  (radiation dominated)


def Y_eq(x, g=1):
    """
    Equilibrium comoving number density Y_eq = n_eq / s.

    For a massive boson with g internal d.o.f.:
        n_eq = g (m T / (2 pi))^(3/2) exp(-m/T)  [non-relativistic]
        s = (2 pi^2 / 45) g_*s T^3

    Y_eq = (45 g) / (4 pi^4 g_*s) * x^(3/2) * e^(-x)    [for x >> 3]

    For x < 3, use the exact Bose-Einstein form:
        Y_eq = (45 / (4 pi^4)) * (g / g_*s) * x^2 * K_2(x)
    where K_2 is the modified Bessel function.
    """
    # Use the Bessel function form which is valid at all x
    return (45.0 / (4.0 * np.pi**4)) * (g / g_star_s) * x**2 * kn(2, x)


def boltzmann_rhs(x, Y_arr, sigma_v):
    """
    Right-hand side of the Boltzmann equation:
        dY/dx = - (lambda / x^2) * (Y^2 - Y_eq^2)

    where lambda = s(m) <sigma v> / H(m) = sqrt(pi/45) g_*s M_Pl m <sigma v> / sqrt(g_*)

    This is derived from:
        s = (2 pi^2 / 45) g_*s m^3 / x^3
        H = sqrt(pi^2 g_* / 90) m^2 / (M_Pl x^2)
        => s / (H x) = (2 pi^2 / 45) g_*s m^3/x^3 / [sqrt(pi^2 g_*/90) m^2/(M_Pl x^2)] / x
                      = sqrt(pi/45) * g_*s * M_Pl * m / (sqrt(g_*) * x^2)
    """
    Y_val = Y_arr[0]

    # The "lambda" parameter (dimensionless measure of interaction strength)
    lam = np.sqrt(np.pi / 45.0) * g_star_s * M_Pl * m_dark_GeV * sigma_v / np.sqrt(g_star)
    # Note: sigma_v is in GeV^-2 for consistency; we convert from cm^3/s below

    Yeq = Y_eq(x, g=g_dm)
    dYdx = -(lam / x**2) * (Y_val**2 - Yeq**2)
    return [dYdx]


def compute_relic_abundance(sigma_v_cm3s, x_start=1.0, x_end=1000.0):
    """
    Solve the Boltzmann equation and return Y_inf, Omega h^2.

    Parameters
    ----------
    sigma_v_cm3s : float
        Thermally averaged annihilation cross-section in cm^3/s.

    Returns
    -------
    x_arr, Y_arr : arrays
        Solution of the Boltzmann equation.
    Omega_h2 : float
        Relic density parameter.
    """
    # Convert sigma_v from cm^3/s to GeV^-2 (natural units)
    # 1 GeV^-2 = 0.3894e-27 cm^2 => 1 cm^2 = 1/(0.3894e-27) GeV^-2
    # sigma*v: [cm^3/s] -> [cm^2 * cm/s] -> [GeV^-2 * (c)]
    # v ~ c, so sigma*v [cm^3/s] ~ sigma [cm^2] * c [cm/s]
    # = sigma [cm^2] * 3e10 [cm/s]
    # sigma [GeV^-2] = sigma [cm^2] / (0.3894e-27 cm^2/GeV^-2)
    # sigma*v [GeV^-2] = sigma_v [cm^3/s] / (0.3894e-27 * 3e10)
    #                   = sigma_v / 1.1682e-17

    hbar_c_cm = 1.9733e-14     # hbar*c in GeV*cm
    conv = hbar_c_cm**2 * 3e10  # = (hbar*c)^2 * c in GeV^2 * cm^3/s ... let me redo properly

    # Correct conversion:
    # In natural units (hbar = c = 1): [sigma v] = GeV^{-2}
    # Physical: sigma [cm^2] = sigma_nat [GeV^-2] * (hbar c)^2 [GeV^2 cm^2]
    # where (hbar c)^2 = (1.9733e-14 GeV cm)^2 = 3.8938e-28 GeV^2 cm^2
    # So sigma_nat = sigma_phys / 3.8938e-28  [GeV^{-2}]
    # sigma*v: (sigma*v)_nat = (sigma*v)_phys / (3.8938e-28 * 3e10)
    #                        = (sigma*v)_phys / 1.1681e-17

    sigma_v_nat = sigma_v_cm3s / 1.1681e-17   # in GeV^{-2}

    # Initial condition: start in equilibrium at x_start
    Y0 = Y_eq(x_start, g=g_dm)

    # Solve ODE
    sol = solve_ivp(
        boltzmann_rhs,
        [x_start, x_end],
        [Y0],
        args=(sigma_v_nat,),
        method='Radau',
        rtol=1e-8,
        atol=1e-12,
        dense_output=True,
        max_step=1.0,
    )

    # Evaluate on a uniform grid for plotting
    x_plot = np.logspace(np.log10(x_start), np.log10(x_end), 500)
    Y_plot = sol.sol(x_plot)[0]

    # Relic abundance
    Y_inf = Y_plot[-1]

    # Omega h^2 = m * s_0 * Y_inf / rho_crit
    # s_0 = 2891.2 cm^{-3}  (present-day entropy density)
    # rho_crit / h^2 = 1.054e-5 GeV/cm^3
    s_0 = 2891.2                # cm^{-3}
    rho_crit_h2 = 1.054e-5     # GeV/cm^3

    Omega_h2 = m_dark_GeV * s_0 * Y_inf / rho_crit_h2

    return x_plot, Y_plot, Omega_h2


# ──────────────────────────────────────────────────────────────
# Main computation
# ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=" * 70)
    print("Boltzmann Freeze-out Calculation")
    print("=" * 70)
    print(f"m_dark = {m_dark_GeV*1e3:.1f} MeV")
    print(f"g_* = {g_star}, g_*s = {g_star_s}, g_DM = {g_dm}")
    print(f"Observed: Omega h^2 = {Omega_obs_h2}")
    print()

    # Cross-sections to scan
    sigma_v_values = [1e-28, 1e-27, 3e-27, 1e-26, 3e-26, 1e-25, 1e-24]
    colors = [PURPLE, CORAL, GOLD, TEAL, 'black', '#264653', '#606060']
    results = []

    print(f"{'<sigma v> (cm3/s)':>22s}  {'Y_inf':>12s}  {'Omega h^2':>12s}  {'Status':>20s}")
    print("-" * 70)

    for sv in sigma_v_values:
        x_arr, Y_arr, Oh2 = compute_relic_abundance(sv)
        results.append((sv, x_arr, Y_arr, Oh2))
        status = ""
        if Oh2 > 0.120 * 1.5:
            status = "OVERPRODUCES"
        elif Oh2 < 0.120 * 0.5:
            status = "UNDERPRODUCES"
        else:
            status = "<-- VIABLE"
        print(f"  {sv:12.2e}         {Y_arr[-1]:12.4e}  {Oh2:12.4f}  {status}")

    # Also compute the exact thermal relic cross-section
    # Omega h^2 ~ 0.120 => <sigma v> ~ 3e-26 cm^3/s
    print()
    print(f"Thermal relic target: <sigma v> = 3e-26 cm^3/s")
    _, _, Oh2_thermal = compute_relic_abundance(3e-26)
    print(f"  => Omega h^2 = {Oh2_thermal:.4f}")

    # ──────────────────────────────────────────────────────────
    # Plot 1: Y(x) vs x for different cross-sections
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 7))

    # Equilibrium curve
    x_eq = np.logspace(0, 3, 500)
    Y_eq_arr = np.array([Y_eq(xi, g=g_dm) for xi in x_eq])
    ax.loglog(x_eq, Y_eq_arr, 'k--', linewidth=2.5, label="$Y_{\\rm eq}$", alpha=0.7)

    for i, (sv, x_arr, Y_arr, Oh2) in enumerate(results):
        label_sv = f"$\\langle\\sigma v\\rangle = {sv:.0e}$"
        label_sv = label_sv.replace("e-0", r"\times 10^{-").replace("e+0", r"\times 10^{")
        label_sv = label_sv.replace("e-", r"\times 10^{-")
        # Simpler labels
        exp = int(np.log10(sv))
        coeff = sv / 10**exp
        if abs(coeff - 1.0) < 0.01:
            label_sv = f"$10^{{{exp}}}$"
        elif abs(coeff - 3.0) < 0.01:
            label_sv = f"$3\\times 10^{{{exp}}}$"
        else:
            label_sv = f"${coeff:.0f}\\times 10^{{{exp}}}$"

        ax.loglog(x_arr, Y_arr, color=colors[i % len(colors)], linewidth=1.8,
                  label=f"<sv> = {label_sv} cm³/s, Ωh² = {Oh2:.3f}")

    ax.set_xlabel("$x = m/T$", fontsize=14)
    ax.set_ylabel("$Y = n/s$ (comoving number density)", fontsize=14)
    ax.set_title("Boltzmann Freeze-out: Topological Dark Matter ($m = 1$ MeV)",
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1, 1000)
    ax.set_ylim(1e-16, 1e-1)
    ax.legend(fontsize=9, loc='lower left')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    # Annotate freeze-out region
    ax.axvspan(15, 30, alpha=0.08, color=GOLD)
    ax.text(22, 5e-2, "Freeze-out\nregion", fontsize=10, ha='center',
            color=GOLD, fontstyle='italic', fontweight='bold')

    plt.tight_layout()
    outpath = "C:/Users/alexn/OneDrive/Documents/v2/Research/relic_abundance_Y_vs_x.png"
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved: {outpath}")
    plt.close()

    # ──────────────────────────────────────────────────────────
    # Plot 2: Omega h^2 vs <sigma v>
    # ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))

    # Fine scan
    sv_scan = np.logspace(-28, -23, 80)
    Oh2_scan = []
    for sv in sv_scan:
        _, _, oh2 = compute_relic_abundance(sv)
        Oh2_scan.append(oh2)
    Oh2_scan = np.array(Oh2_scan)

    ax.loglog(sv_scan, Oh2_scan, color=TEAL, linewidth=2.5, label="Boltzmann freeze-out")

    # Observed value
    ax.axhline(y=Omega_obs_h2, color=CORAL, linewidth=2, linestyle='--',
               label=f"Observed: $\\Omega h^2 = {Omega_obs_h2}$")
    ax.fill_between(sv_scan, Omega_obs_h2 * 0.99, Omega_obs_h2 * 1.01,
                    color=CORAL, alpha=0.2)

    # WIMP miracle cross-section
    ax.axvline(x=3e-26, color=PURPLE, linewidth=1.5, linestyle=':',
               label="Thermal relic: $3\\times 10^{-26}$ cm$^3$/s")

    # Framework's predicted range (from Section 10.4)
    ax.axvspan(1e-27, 1e-25, alpha=0.1, color=GOLD)
    ax.text(3e-26, 5, "Framework's\npredicted range", fontsize=10,
            ha='center', color=GOLD, fontweight='bold', fontstyle='italic')

    # Overproduction / underproduction regions
    ax.fill_between(sv_scan, Omega_obs_h2, 100, alpha=0.05, color='red')
    ax.fill_between(sv_scan, 1e-5, Omega_obs_h2, alpha=0.05, color='blue')
    ax.text(1e-28, 3, "Overproduces DM", fontsize=10, color='red', alpha=0.5)
    ax.text(1e-24, 0.003, "Underproduces DM", fontsize=10, color='blue', alpha=0.5)

    ax.set_xlabel("$\\langle\\sigma v\\rangle$ (cm$^3$/s)", fontsize=14)
    ax.set_ylabel("$\\Omega_{\\rm DM} h^2$", fontsize=14)
    ax.set_title("Relic Abundance vs Annihilation Cross-Section ($m_{\\rm dark} = 1$ MeV)",
                 fontsize=14, fontweight='bold')
    ax.set_xlim(1e-28, 1e-23)
    ax.set_ylim(1e-3, 100)
    ax.legend(fontsize=11, loc='upper right')
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(labelsize=12)

    plt.tight_layout()
    outpath2 = "C:/Users/alexn/OneDrive/Documents/v2/Research/relic_abundance_omega_vs_sigma.png"
    plt.savefig(outpath2, dpi=150, bbox_inches='tight')
    print(f"Plot saved: {outpath2}")
    plt.close()

    print("\nDone.")
