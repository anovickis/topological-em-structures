#!/usr/bin/env python3
"""
sim_coupling_threshold.py -- One-Loop Threshold Matching for FN Couplings

Computes the effective Faddeev-Niemi couplings kappa_2 and kappa_4 arising
from the Cho-Faddeev-Niemi (CFN) decomposition of SU(2)_L x U(1)_Y
electroweak theory, including one-loop corrections from:
  - W-boson loops
  - Higgs boson loops
  - Goldstone boson loops (in R_xi gauge)
  - Top quark loops (dominant fermion contribution)
  - Ghost contributions

The key question: can g_W = 0.65 produce kappa_2 ~ alpha*hbar*c through
radiative corrections, or is the tree-level "factor of 20 discrepancy"
a genuine problem requiring non-perturbative physics?

Physical setup:
  In the CFN decomposition, the SU(2) gauge field decomposes as:
    A_mu^a = C_mu n^a - (1/g) eps^abc n^b d_mu n^c + W_mu^a

  The effective action for the n-field, after integrating out W, H, ghosts:
    S_eff[n] = integral d^4x [kappa_2/2 (d_mu n)^2 + kappa_4/4 (d_mu n x d_nu n)^2]

  Tree level: kappa_2 = 1/g^2 (from Yang-Mills kinetic term)

  One-loop corrections come from:
    (a) W-boson self-energy with n-field legs (gauge sector)
    (b) Higgs/Goldstone loops (scalar sector, modified by EWSB)
    (c) Fermion loops (top quark dominant)
    (d) Ghost loops (Faddeev-Popov)

  Additionally, in the electroweak theory with Higgs VEV v = 246 GeV,
  the Goldstone boson kinetic term |D_mu phi|^2 generates an ADDITIONAL
  tree-level contribution to kappa_2 proportional to v^2.

References:
  [53] Cho, Y.M. (1980). Phys. Rev. D 21, 1080.
  [54] Faddeev & Niemi (1999). Phys. Rev. Lett. 82, 1624.
  Paper: Toroidal_Electron_Full_Paper_R10.md, Section 13.3.1
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

# ─── Color palette ────────────────────────────────────────────────
CORAL  = '#e76f51'
TEAL   = '#2a9d8f'
GOLD   = '#e9c46a'
PURPLE = '#a855f7'
BG_COLOR = '#f8f9fa'

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'images')

# ═══════════════════════════════════════════════════════════════════
#  1. PHYSICAL CONSTANTS (NATURAL UNITS: hbar = c = 1)
# ═══════════════════════════════════════════════════════════════════

# All masses in GeV
M_W   = 80.379       # W boson mass
M_Z   = 91.1876      # Z boson mass
M_H   = 125.25       # Higgs boson mass
M_T   = 172.76       # Top quark mass
M_E   = 0.51100e-3   # Electron mass
V_EW  = 246.22       # Higgs VEV

# Couplings at M_Z scale
ALPHA_EM = 1.0 / 137.036    # Fine structure constant
G_W      = 0.6530           # SU(2)_L coupling at M_Z
G_PRIME  = 0.3492           # U(1)_Y coupling at M_Z
ALPHA_S  = 0.1179           # Strong coupling at M_Z (for comparison)

# Derived
ALPHA_W = G_W**2 / (4 * np.pi)    # SU(2) fine structure constant
THETA_W = np.arctan(G_PRIME / G_W) # Weinberg angle
SIN2_TW = np.sin(THETA_W)**2

# Top Yukawa
Y_T = np.sqrt(2) * M_T / V_EW


# ═══════════════════════════════════════════════════════════════════
#  2. RUNNING COUPLING (ONE-LOOP RG)
# ═══════════════════════════════════════════════════════════════════

def run_g_W(mu, g0=G_W, mu0=M_Z):
    """
    One-loop running of SU(2)_L coupling.

    Beta function: b = -19/6 + n_H/6 + 2*n_f/3
    For SM (1 Higgs doublet, 3 generations): b_2 = -19/6
    (asymptotically free but very slow running)

    g^{-2}(mu) = g^{-2}(mu0) + b/(8 pi^2) * ln(mu/mu0)
    """
    b2 = -19.0 / 6.0  # One-loop beta coefficient for SU(2)_L in SM

    # Avoid log of negative/zero
    if mu <= 0:
        mu = 1e-10

    g2_inv = 1.0 / g0**2 + b2 / (8 * np.pi**2) * np.log(mu / mu0)

    if g2_inv <= 0:
        return np.inf  # Landau pole (never reached for asymptotically free)

    return 1.0 / np.sqrt(g2_inv)


def run_alpha_em(mu, alpha0=ALPHA_EM, mu0=M_Z):
    """
    One-loop running of alpha_EM.

    Below M_W, only charged leptons and light quarks contribute.
    Above M_W, full SM running.
    """
    # Simplified: use the QED beta function with n_f light fermions
    # b_QED = -4/3 * sum_f N_c Q_f^2
    # Below M_W: e, mu, tau (Q=1), u,d,s,c,b (Q=2/3 or 1/3, N_c=3)
    # Sum of N_c Q_f^2 = 3*(1) + 3*(3*(4/9) + 3*(1/9)) = 3 + 5 = 20/3
    # b = -4/3 * 20/3 = -80/9

    b = -80.0 / 9.0
    alpha_inv = 1.0 / alpha0 + b / (2 * np.pi) * np.log(mu / mu0)

    if alpha_inv <= 0:
        return 1.0  # cap
    return 1.0 / alpha_inv


# ═══════════════════════════════════════════════════════════════════
#  3. TREE-LEVEL MATCHING
# ═══════════════════════════════════════════════════════════════════

def kappa2_tree(g):
    """
    Tree-level FN coupling from CFN decomposition.
    kappa_2 = 1/g^2  (dimensionless in natural units)
    """
    return 1.0 / g**2


def kappa4_tree(g):
    """
    Tree-level quartic FN coupling.
    kappa_4 = 1/g^4
    """
    return 1.0 / g**4


# ═══════════════════════════════════════════════════════════════════
#  4. ONE-LOOP CORRECTIONS TO KAPPA_2
# ═══════════════════════════════════════════════════════════════════

def delta_kappa2_gauge(g, mu, M_W_val=M_W):
    """
    W-boson loop contribution to kappa_2.

    The W-boson propagating in a loop with two n-field external legs.
    In the background field method, this is the gauge boson contribution
    to the n-field self-energy.

    For SU(2) pure gauge theory, the one-loop correction to the n-field
    kinetic term is (Faddeev-Niemi, 1999):
        delta_kappa2 = (g^2 / (16 pi^2)) * c_gauge * ln(Lambda/mu)

    With the W-boson mass as the UV cutoff Lambda -> M_W:
        delta_kappa2 = (g^2 / (16 pi^2)) * c_gauge * ln(M_W/mu)

    The coefficient c_gauge for SU(2):
    From the Yang-Mills one-loop effective action in the n-field background,
    the correction to the (dn)^2 operator comes from:
      - W-boson loop with 2 n-field vertices: coefficient = 11/3 (from
        the asymptotic freedom coefficient of SU(2))
    But the sign is POSITIVE (anti-screening for the n-field), meaning
    kappa_2 INCREASES at low energies.
    """
    c_gauge = 22.0 / 3.0  # = 2 * (11/3) for SU(2), from the n-field sector
    # Note: the factor 22/3 comes from the fact that the n-field kinetic
    # term receives contributions from both transverse W polarizations
    # (each contributing 11/3 from the gauge field beta function)

    if mu <= 0 or mu >= M_W_val:
        return 0.0

    return (g**2 / (16 * np.pi**2)) * c_gauge * np.log(M_W_val / mu)


def delta_kappa2_higgs(g, mu, M_H_val=M_H, M_W_val=M_W, v=V_EW):
    """
    Higgs boson loop contribution to kappa_2.

    The Higgs couples to the n-field through the covariant derivative
    |D_mu phi|^2 = ... + (g^2 v^2 / 4) (dn)^2 + (g^2/4)(h^2)(dn)^2 + ...

    The Higgs loop generates a correction:
        delta_kappa2 = (g^2 / (16 pi^2)) * (M_H^2 / M_W^2) * c_H * ln(M_H/mu)

    The coefficient c_H accounts for the quartic Higgs-n coupling.
    In the SM, lambda = M_H^2 / (2 v^2), and the Higgs contribution is:
        c_H = 1/4  (from the single Higgs doublet with 4 real components)

    But more importantly, the Higgs VEV contributes a TREE-LEVEL term:
        kappa_2^(Higgs, tree) = v^2 / (4 * (something with mass scale))
    This is actually absorbed into the definition of the n-field normalization.
    The LOOP correction is what we compute here.
    """
    # Higgs-n coupling from |D phi|^2 vertex
    lambda_H = M_H_val**2 / (2 * v**2)
    c_H = 0.25  # coefficient from the Higgs loop diagram

    if mu <= 0 or mu >= M_H_val:
        return 0.0

    return (g**2 / (16 * np.pi**2)) * c_H * np.log(M_H_val / mu)


def delta_kappa2_goldstone(g, mu, M_W_val=M_W):
    """
    Goldstone boson loop contribution to kappa_2.

    In R_xi gauge, the Goldstone bosons (eaten by W, Z) propagate with
    mass M_G = xi * M_W. In Landau gauge (xi = 0), they are massless;
    in unitary gauge (xi -> inf), they decouple.

    We work in Landau gauge where the Goldstone contribution is maximal.
    The 3 Goldstone bosons (pi^a for a=1,2,3) couple to n through:
        g^2 v^2 / 4 * pi^a pi^b * (n-dependent terms)

    This gives a contribution:
        delta_kappa2 = (3 * g^2) / (16 pi^2) * c_G * ln(M_W/mu)
    where c_G = 1/6 (from the triangle diagram with two Goldstone propagators)
    and the factor of 3 counts the three Goldstone bosons.
    """
    c_G = 1.0 / 6.0
    n_G = 3  # number of Goldstone bosons

    if mu <= 0 or mu >= M_W_val:
        return 0.0

    return n_G * (g**2 / (16 * np.pi**2)) * c_G * np.log(M_W_val / mu)


def delta_kappa2_fermion(g, mu, m_f, N_c=3):
    """
    Fermion loop contribution to kappa_2.

    A fermion of mass m_f and SU(2) coupling g contributes:
        delta_kappa2 = -N_c * (g^2 / (16 pi^2)) * c_f * ln(m_f/mu)

    The coefficient c_f = 2/3 for a Dirac fermion in the fundamental
    of SU(2). The negative sign is from the fermion loop (opposite to bosons).
    N_c = 3 for quarks (color multiplicity), N_c = 1 for leptons.
    """
    c_f = 2.0 / 3.0

    if mu <= 0 or mu >= m_f:
        return 0.0

    # Negative sign: fermions screen
    return -N_c * (g**2 / (16 * np.pi**2)) * c_f * np.log(m_f / mu)


def delta_kappa2_ghost(g, mu, M_W_val=M_W):
    """
    Faddeev-Popov ghost contribution to kappa_2.

    Ghosts contribute with a negative sign and coefficient:
        delta_kappa2 = -(g^2 / (16 pi^2)) * c_ghost * ln(M_W/mu)

    For SU(2): c_ghost = 1/3 (standard ghost contribution to beta function)
    """
    c_ghost = 1.0 / 3.0

    if mu <= 0 or mu >= M_W_val:
        return 0.0

    return -(g**2 / (16 * np.pi**2)) * c_ghost * np.log(M_W_val / mu)


def kappa2_one_loop(mu, g=None, include_higgs_vev=False):
    """
    Total kappa_2 at one-loop order at scale mu.

    Parameters
    ----------
    mu : float
        Renormalization scale (GeV)
    g : float or None
        SU(2) coupling at M_Z (default: G_W)
    include_higgs_vev : bool
        If True, include the tree-level Higgs VEV contribution
        to the n-field kinetic term. This is the key question:
        does the Higgs mechanism provide an additional kappa_2
        beyond 1/g^2?

    Returns
    -------
    kappa2 : float
        Effective kappa_2 at scale mu
    contributions : dict
        Individual contributions broken down
    """
    if g is None:
        g = G_W

    # Run the coupling to the matching scale (take geometric mean of
    # heavy particle masses as matching scale)
    mu_match = np.sqrt(M_W * M_H)  # ~ 100 GeV
    g_at_match = run_g_W(mu_match, g)

    # Tree level
    k2_tree = kappa2_tree(g_at_match)

    # One-loop corrections (evaluated at scale mu)
    dk2_gauge = delta_kappa2_gauge(g_at_match, mu)
    dk2_higgs = delta_kappa2_higgs(g_at_match, mu)
    dk2_gold  = delta_kappa2_goldstone(g_at_match, mu)
    dk2_ghost = delta_kappa2_ghost(g_at_match, mu)

    # Fermion contributions (top quark dominant)
    dk2_top = delta_kappa2_fermion(g_at_match, mu, M_T, N_c=3)
    # Bottom quark (small contribution)
    dk2_bot = delta_kappa2_fermion(g_at_match, mu, 4.18, N_c=3)
    # Tau lepton
    dk2_tau = delta_kappa2_fermion(g_at_match, mu, 1.777, N_c=1)

    dk2_fermions = dk2_top + dk2_bot + dk2_tau

    # Total one-loop
    k2_total = k2_tree + dk2_gauge + dk2_higgs + dk2_gold + dk2_ghost + dk2_fermions

    # Optional: Higgs VEV tree-level contribution
    # In the electroweak theory, |D_mu phi|^2 generates a term
    # (v^2/4) * g^2 * (dn)^2 when the Goldstone modes are expressed
    # through the n-field. This is a tree-level contribution.
    # However, its correct normalization depends on the precise
    # relationship between the Goldstone field and the n-field,
    # which involves the soliton scale.
    k2_higgs_vev = 0.0
    if include_higgs_vev:
        # The Higgs VEV contribution in dimensionless units
        # (normalized to the soliton scale R ~ hbar/(m_e c)):
        # kappa_2^(vev) ~ (v * R)^2 / 4 = (v / m_e)^2 / 4
        # This is HUGE: (246 GeV / 0.511 MeV)^2 / 4 ~ 5.8e10
        # But this overcounts -- the Goldstone modes ARE the W longitudinal
        # components, which are already integrated out in the gauge contribution.
        # The correct procedure is that the Goldstone contribution replaces
        # part of the gauge contribution in unitary gauge.
        # We flag this as an important subtlety.
        k2_higgs_vev = 0.0  # Set to 0: see discussion below

    k2_total += k2_higgs_vev

    contributions = {
        'tree': k2_tree,
        'gauge (W loops)': dk2_gauge,
        'higgs': dk2_higgs,
        'goldstone': dk2_gold,
        'ghost': dk2_ghost,
        'top quark': dk2_top,
        'bottom quark': dk2_bot,
        'tau lepton': dk2_tau,
        'fermions (total)': dk2_fermions,
        'higgs VEV (tree)': k2_higgs_vev,
        'total': k2_total,
    }

    return k2_total, contributions


# ═══════════════════════════════════════════════════════════════════
#  5. COMPARISON WITH REQUIRED VALUE
# ═══════════════════════════════════════════════════════════════════

def required_kappa2():
    """
    What kappa_2 is needed to match the electron?

    From the paper:
        E_min = 192.5 * sqrt(kappa_2 * kappa_4) = m_e c^2
        R = sqrt(kappa_4 / kappa_2) = hbar/(m_e c)

    Solving:
        kappa_4 = kappa_2 * R^2
        kappa_2 * sqrt(kappa_2 * R^2) = (m_e / 192.5)
        kappa_2^(3/2) * R = m_e / 192.5
        kappa_2 = (m_e / (192.5 * R))^(2/3)

    But in natural units (hbar = c = 1):
        R = 1/m_e
        kappa_2 = (m_e * m_e / 192.5)^(2/3) = (m_e^2 / 192.5)^(2/3)

    Wait, let's redo more carefully.
    E_min = 192.5 * sqrt(k2 * k4) = m_e
    R = sqrt(k4/k2) = 1/m_e

    From R: k4 = k2 / m_e^2
    From E: 192.5 * sqrt(k2 * k2/m_e^2) = m_e
             192.5 * k2/m_e = m_e
             k2 = m_e^2 / 192.5

    In natural units (GeV):
        k2 = (0.000511)^2 / 192.5 = 1.36e-9 GeV^2

    But in the paper's convention, kappa_2 is dimensionless (since
    the n-field is dimensionless and dn/dx has dimension 1/length = mass).
    The FN Lagrangian is L = kappa_2/2 (dn)^2 + ... with L having
    dimension mass^4 and (dn)^2 having dimension mass^2, so
    kappa_2 has dimension mass^2 = GeV^2.

    If we want to compare with the tree-level formula kappa_2 = 1/g^2
    (dimensionless), we need to account for the mass scale.

    The resolution: in the CFN decomposition of the Yang-Mills action,
    S = int d^4x (-1/4g^2) F^2, the kinetic term for n is:
        (1/2g^2) (d_mu n)^2
    This is written with the canonical normalization where (d_mu n)^2
    has dimension mass^2 (in natural units), so the coefficient 1/g^2
    is dimensionless. The physical kappa_2 = 1/g^2 IS dimensionless.

    But E = int d^3x [...] gives energy in GeV. And:
        int d^3x (1/2g^2)(dn)^2 has dimension (1/mass^3)(1/mass^2) = 1/mass^5
    That's not energy!

    The issue is that in the FN model as used in the paper, the fields
    are measured in "soliton units" where kappa_2 = kappa_4 = 1, and the
    soliton size is R = sqrt(kappa_4/kappa_2) = 1.

    For the physical matching:
        kappa_2 [in the paper's convention] has NO specific dimension;
        the ratio kappa_2 ~ alpha * hbar * c is a statement about
        the numerical value in appropriate units.

    To make a meaningful comparison, we compute the RATIO:
        kappa_2(computed) / kappa_2(required for electron matching)

    If this ratio is 1, the coupling produces the right electron.
    """
    # The paper states kappa_2 ~ alpha * hbar * c (Eq. 13.16 context)
    # In the tree-level formula: kappa_2 = 1/g^2
    # The "required" value is kappa_2 ~ 1/alpha ~ 137
    # (because if kappa_2 = 1/g^2 and kappa_2 ~ alpha*hbar*c,
    #  then g^2 ~ 1/(alpha*hbar*c) ~ 137 in natural units)
    #
    # But this is precisely the discrepancy! g = sqrt(137) ~ 12, not 0.65.
    #
    # What we want to know: can one-loop corrections bring
    # kappa_2 from 1/g_W^2 ~ 2.35 up to ~ 137?

    k2_tree_gW = 1.0 / G_W**2
    k2_required = 1.0 / ALPHA_EM  # ~ 137

    return k2_required, k2_tree_gW


# ═══════════════════════════════════════════════════════════════════
#  6. ENHANCEMENT FACTOR ANALYSIS
# ═══════════════════════════════════════════════════════════════════

def compute_enhancement_ratio(mu):
    """
    Compute the ratio kappa_2(one-loop) / kappa_2(tree) at scale mu.
    This tells us how much the one-loop corrections enhance kappa_2.
    """
    k2_total, _ = kappa2_one_loop(mu)
    k2_tree = kappa2_tree(G_W)
    return k2_total / k2_tree


def scan_mu_range(mu_min=1e-4, mu_max=100, n_points=500):
    """
    Scan kappa_2 as a function of renormalization scale mu.
    """
    mu_values = np.geomspace(mu_min, mu_max, n_points)
    k2_values = []
    contributions_list = []

    for mu in mu_values:
        k2, contribs = kappa2_one_loop(mu)
        k2_values.append(k2)
        contributions_list.append(contribs)

    return mu_values, np.array(k2_values), contributions_list


# ═══════════════════════════════════════════════════════════════════
#  7. NON-PERTURBATIVE ESTIMATE: DIMENSIONAL TRANSMUTATION
# ═══════════════════════════════════════════════════════════════════

def lambda_su2(g=G_W, mu0=M_Z):
    """
    Estimate the SU(2) confinement scale Lambda_SU2 via dimensional
    transmutation (one-loop formula).

    Lambda = mu0 * exp(-8 pi^2 / (b * g^2(mu0)))

    For SU(2) with SM matter content: b = -19/6
    (Note: b < 0 means asymptotic freedom)

    Lambda is the scale where the coupling becomes strong (g -> inf).
    If Lambda >> m_e, non-perturbative effects dominate the soliton physics.
    """
    b = -19.0 / 6.0
    Lambda = mu0 * np.exp(-8 * np.pi**2 / (abs(b) * g**2))
    return Lambda


def estimate_nonperturbative_kappa2():
    """
    Estimate the non-perturbative contribution to kappa_2.

    By analogy with QCD and the Skyrme model:
      - In QCD, f_pi = 93 MeV is a non-perturbative scale
      - The Skyrme model coupling is e ~ 5.45, f_pi / (2e) ~ 8.5 MeV
      - The pion kinetic term has coefficient f_pi^2 / 4

    For the electroweak theory:
      - The analogue of f_pi is the Higgs VEV: v = 246 GeV
      - The n-field kinetic term from the Higgs sector is ~ v^2
      - But this must be compared in the same units as 1/g^2

    The key insight: in the electroweak theory, the Higgs mechanism
    provides a LARGE contribution to the n-field kinetic energy
    through the Goldstone boson kinetic terms. In unitary gauge,
    this appears as part of the W-boson mass term.

    The total kappa_2 is schematically:
        kappa_2 ~ 1/g^2 + C * (v/Lambda_UV)^2
    where Lambda_UV is the UV cutoff of the FN effective theory.

    If Lambda_UV ~ M_W ~ gv/2, then:
        (v/Lambda_UV)^2 ~ 4/g^2
    and the Higgs contribution is of the same order as the tree-level term.

    But if the relevant UV scale is higher (e.g., GUT scale), the
    contribution could be much larger.
    """
    # Compute Lambda_SU2
    L = lambda_su2()

    results = {
        'Lambda_SU2': L,
        'Lambda_QCD_comparison': 0.200,  # GeV
        'v_EW': V_EW,
        'ratio_v_Lambda': V_EW / L if L > 0 else np.inf,
        'kappa2_tree': 1.0 / G_W**2,
        'kappa2_higgs_estimate': 4.0 / G_W**2,  # from v^2/M_W^2 ~ 4/g^2
    }

    return results


# ═══════════════════════════════════════════════════════════════════
#  8. PLOTTING
# ═══════════════════════════════════════════════════════════════════

def plot_kappa2_vs_mu(mu_values, k2_values, outdir):
    """Plot kappa_2 as function of renormalization scale."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG_COLOR)

    # Left: kappa_2 vs mu
    ax = axes[0]
    ax.set_facecolor('white')
    ax.semilogx(mu_values, k2_values, color=PURPLE, linewidth=2.5,
                label='$\\kappa_2^{\\rm 1-loop}(\\mu)$')
    ax.axhline(1.0 / G_W**2, color=TEAL, linestyle='--', linewidth=1.5,
               label=f'Tree level: $1/g_W^2 = {1/G_W**2:.2f}$')
    ax.axhline(1.0 / ALPHA_EM, color=CORAL, linestyle=':', linewidth=1.5,
               label=f'Required: $1/\\alpha = {1/ALPHA_EM:.1f}$')

    # Mark electron mass scale
    ax.axvline(M_E, color=GOLD, linestyle='-.', linewidth=1, alpha=0.7)
    ax.text(M_E * 1.5, ax.get_ylim()[0] + 0.5, '$m_e$',
            color=GOLD, fontsize=11, fontweight='bold')

    ax.set_xlabel('$\\mu$ (GeV)', fontsize=13)
    ax.set_ylabel('$\\kappa_2(\\mu)$', fontsize=13)
    ax.set_title('$\\kappa_2$ vs. Renormalization Scale', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1e-4, 100)

    # Right: enhancement ratio
    ax = axes[1]
    ax.set_facecolor('white')
    ratio = k2_values / (1.0 / G_W**2)
    ax.semilogx(mu_values, ratio, color=PURPLE, linewidth=2.5)
    ax.axhline(1.0, color=TEAL, linestyle='--', linewidth=1.5,
               label='No correction (ratio = 1)')
    ax.axhline(1.0 / (ALPHA_EM * G_W**2), color=CORAL, linestyle=':', linewidth=1.5,
               label=f'Required ratio: {1/(ALPHA_EM * G_W**2):.1f}')
    ax.axvline(M_E, color=GOLD, linestyle='-.', linewidth=1, alpha=0.7)

    ax.set_xlabel('$\\mu$ (GeV)', fontsize=13)
    ax.set_ylabel('$\\kappa_2^{\\rm 1-loop} / \\kappa_2^{\\rm tree}$', fontsize=13)
    ax.set_title('Enhancement Ratio', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(1e-4, 100)

    plt.tight_layout()
    path = os.path.join(outdir, 'coupling_threshold_kappa2.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_contributions_breakdown(mu_target, outdir):
    """Bar chart of individual contributions at the soliton scale."""
    _, contribs = kappa2_one_loop(mu_target)

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    # Select contributions to plot
    labels = ['Tree\n($1/g_W^2$)', 'W loops', 'Higgs\nloop', 'Goldstone\nloops',
              'Ghost', 'Top\nquark', 'Other\nfermions']
    values = [
        contribs['tree'],
        contribs['gauge (W loops)'],
        contribs['higgs'],
        contribs['goldstone'],
        contribs['ghost'],
        contribs['top quark'],
        contribs['bottom quark'] + contribs['tau lepton'],
    ]

    colors_bar = [PURPLE, TEAL, GOLD, TEAL, CORAL, CORAL, CORAL]

    bars = ax.bar(range(len(labels)), values, color=colors_bar, alpha=0.85,
                  edgecolor='white', linewidth=1.5, width=0.7)

    # Value labels
    for bar, val in zip(bars, values):
        y = bar.get_height()
        va = 'bottom' if y >= 0 else 'top'
        offset = 0.02 if y >= 0 else -0.02
        ax.text(bar.get_x() + bar.get_width() / 2, y + offset,
                f'{val:.4f}', ha='center', va=va,
                fontsize=9, fontweight='bold')

    # Reference lines
    ax.axhline(1.0 / ALPHA_EM, color=CORAL, linestyle=':', linewidth=1.5,
               label=f'Required: $1/\\alpha = {1/ALPHA_EM:.1f}$')
    ax.axhline(0, color='black', linewidth=0.5)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('Contribution to $\\kappa_2$', fontsize=13)
    ax.set_title(f'One-Loop Contributions to $\\kappa_2$ at $\\mu = {mu_target:.1e}$ GeV\n'
                 f'(Electron mass scale)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')

    # Add total annotation
    total = contribs['total']
    ax.text(0.02, 0.95, f'Total $\\kappa_2$ = {total:.4f}\n'
            f'Required = {1/ALPHA_EM:.1f}\n'
            f'Gap factor = {1/(ALPHA_EM * total):.1f}$\\times$',
            transform=ax.transAxes, fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    path = os.path.join(outdir, 'coupling_threshold_breakdown.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_gap_analysis(outdir):
    """
    Plot showing what non-perturbative enhancement factor is needed.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    # Compare different coupling scenarios
    g_values = np.linspace(0.3, 15, 300)
    k2_tree_vals = 1.0 / g_values**2

    ax.semilogy(g_values, k2_tree_vals, color=PURPLE, linewidth=2.5,
                label='Tree level: $\\kappa_2 = 1/g^2$')

    # Required value
    ax.axhline(1.0 / ALPHA_EM, color=CORAL, linestyle=':', linewidth=2,
               label=f'Required: $1/\\alpha \\approx {1/ALPHA_EM:.0f}$')

    # Mark g_W
    ax.axvline(G_W, color=TEAL, linestyle='--', linewidth=1.5, alpha=0.7)
    ax.plot(G_W, 1.0 / G_W**2, 'o', markersize=12, markerfacecolor=TEAL,
            markeredgecolor='white', markeredgewidth=2,
            label=f'$g_W = {G_W}$: $\\kappa_2 = {1/G_W**2:.2f}$')

    # Mark g ~ 12 (the "required" coupling)
    g_req = 1.0 / np.sqrt(ALPHA_EM)
    ax.axvline(g_req, color=CORAL, linestyle='--', linewidth=1.5, alpha=0.7)
    ax.plot(g_req, 1.0 / g_req**2, 's', markersize=12, markerfacecolor=CORAL,
            markeredgecolor='white', markeredgewidth=2,
            label=f'$g \\approx {g_req:.1f}$: matches $\\alpha$')

    # One-loop enhanced kappa_2 at mu = m_e
    k2_1loop, _ = kappa2_one_loop(M_E)
    ax.plot(G_W, k2_1loop, 'D', markersize=12, markerfacecolor=GOLD,
            markeredgecolor='white', markeredgewidth=2,
            label=f'One-loop at $\\mu = m_e$: $\\kappa_2 = {k2_1loop:.2f}$')

    # Shade the "gap"
    ax.fill_between([0.3, 1.0], [k2_1loop, k2_1loop],
                    [1/ALPHA_EM, 1/ALPHA_EM],
                    alpha=0.15, color=CORAL, label='Non-perturbative gap')

    # Annotations
    gap = (1.0 / ALPHA_EM) / k2_1loop
    ax.annotate(f'Gap: {gap:.0f}$\\times$\n(non-perturbative)',
                xy=(0.65, np.sqrt(k2_1loop * 1/ALPHA_EM)),
                fontsize=12, fontweight='bold', color=CORAL,
                ha='center')

    ax.set_xlabel('$g$ (SU(2) coupling)', fontsize=13)
    ax.set_ylabel('$\\kappa_2$', fontsize=13)
    ax.set_title('Coupling Constant Matching: Perturbative vs Required $\\kappa_2$',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.3, 15)
    ax.set_ylim(1e-3, 200)

    plt.tight_layout()
    path = os.path.join(outdir, 'coupling_gap_analysis.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ═══════════════════════════════════════════════════════════════════
#  9. MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    print("=" * 70)
    print("  One-Loop Threshold Matching for FN Couplings")
    print("  CFN Decomposition of Electroweak Theory")
    print("=" * 70)
    print()

    # --- Physical parameters ---
    print("Physical parameters:")
    print(f"  g_W        = {G_W:.4f} (SU(2)_L coupling at M_Z)")
    print(f"  g'         = {G_PRIME:.4f} (U(1)_Y coupling at M_Z)")
    print(f"  alpha_EM   = 1/{1/ALPHA_EM:.3f}")
    print(f"  sin^2(tw)  = {SIN2_TW:.4f}")
    print(f"  M_W        = {M_W:.3f} GeV")
    print(f"  M_Z        = {M_Z:.4f} GeV")
    print(f"  M_H        = {M_H:.2f} GeV")
    print(f"  M_T        = {M_T:.2f} GeV")
    print(f"  v          = {V_EW:.2f} GeV")
    print(f"  m_e        = {M_E:.5e} GeV")
    print(f"  y_t        = {Y_T:.4f}")
    print()

    # --- Running coupling ---
    print("Running coupling g_W(mu):")
    for mu in [M_Z, 10, 1, 0.1, M_E]:
        g = run_g_W(mu)
        print(f"  g_W({mu:.4g} GeV) = {g:.4f}")
    print()

    # --- Tree level ---
    k2_tree = kappa2_tree(G_W)
    k2_required = 1.0 / ALPHA_EM
    print("Tree-level matching:")
    print(f"  kappa_2(tree) = 1/g_W^2 = {k2_tree:.4f}")
    print(f"  kappa_2(required) = 1/alpha = {k2_required:.2f}")
    print(f"  Ratio: {k2_required / k2_tree:.1f}x (this is the 'discrepancy')")
    print()

    # --- One-loop at electron scale ---
    mu_soliton = M_E
    k2_1loop, contribs = kappa2_one_loop(mu_soliton)

    print(f"One-loop matching at mu = {mu_soliton:.4e} GeV (electron scale):")
    print("-" * 55)
    for name, val in contribs.items():
        if name not in ['total', 'fermions (total)', 'higgs VEV (tree)']:
            pct = val / k2_tree * 100
            print(f"  {name:25s}: {val:+.6f}  ({pct:+.2f}% of tree)")
    print("-" * 55)
    print(f"  {'TOTAL':25s}: {k2_1loop:.6f}")
    print(f"  Required:                   {k2_required:.2f}")
    print(f"  Enhancement from 1-loop:    {k2_1loop / k2_tree:.4f}x")
    print(f"  Remaining gap:              {k2_required / k2_1loop:.1f}x")
    print()

    # --- Loop expansion parameter check ---
    g_eff = G_W
    eps_loop = g_eff**2 / (16 * np.pi**2)
    print("Perturbative validity check:")
    print(f"  g_W^2 / (16 pi^2) = {eps_loop:.6f}")
    print(f"  This is SMALL -- perturbation theory is valid for g_W = {G_W}")
    print(f"  One-loop corrections are O({eps_loop:.4f}) relative to tree level")
    print(f"  Total one-loop correction: {(k2_1loop - k2_tree) / k2_tree:.4f} = "
          f"{(k2_1loop - k2_tree) / k2_tree / eps_loop:.2f} * epsilon_loop")
    print()

    # --- The inferred g ~ 12 check ---
    g_inferred = 1.0 / np.sqrt(ALPHA_EM)
    eps_inferred = g_inferred**2 / (16 * np.pi**2)
    print(f"If g = 1/sqrt(alpha) = {g_inferred:.2f}:")
    print(f"  g^2 / (16 pi^2) = {eps_inferred:.4f}")
    print(f"  This is O(1) -- perturbation theory BREAKS DOWN")
    print(f"  The tree-level formula kappa_2 = 1/g^2 is UNRELIABLE at this g")
    print()

    # --- Lambda_SU2 ---
    L = lambda_su2()
    print("Non-perturbative scale (dimensional transmutation):")
    print(f"  Lambda_SU2 = {L:.4e} GeV")
    if L < 1e-20:
        print(f"  This is ASTRONOMICALLY small (< 10^-20 GeV)")
        print(f"  SU(2)_L never becomes strongly coupled at any physical scale")
        print(f"  --> Non-perturbative enhancement within SU(2)_L alone is negligible")
    print()

    # --- Key conclusion ---
    print("=" * 70)
    print("  KEY CONCLUSIONS")
    print("=" * 70)
    print()
    print("  1. One-loop perturbative corrections to kappa_2 are SMALL")
    print(f"     ({(k2_1loop - k2_tree) / k2_tree * 100:.2f}% enhancement)")
    print(f"     because g_W = {G_W} is weakly coupled.")
    print()
    print(f"  2. The required kappa_2 = 1/alpha ~ {k2_required:.0f} is {k2_required/k2_1loop:.0f}x")
    print(f"     larger than the one-loop value {k2_1loop:.4f}.")
    print()
    print(f"  3. Lambda_SU2 = {L:.2e} GeV -- SU(2)_L never confines.")
    print(f"     Non-perturbative enhancement within pure SU(2)_L is negligible.")
    print()
    print("  4. IMPLICATIONS:")
    print("     (a) If the electron is an FN soliton in SU(2)_L, the coupling")
    print("         kappa_2 ~ alpha*hbar*c CANNOT arise from g_W perturbatively")
    print("         or non-perturbatively within SU(2)_L alone.")
    print()
    print("     (b) The mismatch requires EITHER:")
    print("         - A different, strongly coupled SU(2) (BSM physics)")
    print("         - Higgs-sector contributions that modify the matching")
    print("           beyond the simple CFN formula (needs further study)")
    print("         - Treating kappa_2 as a phenomenological parameter")
    print("           (bottom-up approach) rather than deriving it from g_W")
    print()
    print("     (c) The paper's reframing (R10, Section 13.3.1) is correct:")
    print("         the 'discrepancy' is a perturbative formula applied outside")
    print("         its domain -- but the domain is actually PERTURBATIVE for g_W.")
    print("         The issue is that kappa_2 = 1/g^2 is only the tree-level")
    print("         result, and the FULL non-perturbative kappa_2(g) is unknown.")
    print()

    # --- Generate plots ---
    print("Generating figures...")
    os.makedirs(OUTDIR, exist_ok=True)
    mu_values, k2_values, _ = scan_mu_range()
    plot_kappa2_vs_mu(mu_values, k2_values, OUTDIR)
    plot_contributions_breakdown(M_E, OUTDIR)
    plot_gap_analysis(OUTDIR)

    print()
    print("=" * 70)
    print("  COMPUTATION COMPLETE")
    print("=" * 70)

    return k2_1loop, contribs


if __name__ == '__main__':
    k2, contribs = main()
