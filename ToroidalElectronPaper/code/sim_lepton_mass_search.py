#!/usr/bin/env python3
"""
sim_lepton_mass_search.py -- Systematic Search for Lepton Mass Ratios
                             from Topological Quantum Numbers

Searches for topological quantum number assignments that reproduce the
known lepton mass ratios (electron, muon, tau) within the toroidal
electron / Hopf-fibred soliton framework.

Known mass ratios:
    m_mu / m_e  = 206.768
    m_tau / m_e = 3477.48
    m_tau / m_mu = 16.818

Five search strategies:
    A) Hopf charge scaling  m(H) = m_e * f(H)
    B) Torus knot (p,q) classification
    C) Mass formula extension  m = m_P * alpha^(n1/2 - n2*alpha/4)
    D) Composite topological quantum numbers
    E) Kaluza-Klein / compactification tower

Reference:
    - Toroidal_Electron_Full_Paper.md, Sections 6 and 16
    - Battye & Sutcliffe, Proc. R. Soc. A (1999)

Output:
    - Console tables of best-fit parameters for each approach
    - lepton_mass_hopf_scaling.png
    - lepton_mass_knot_search.png
    - lepton_mass_formula_scan.png
    - lepton_mass_best_fits.png
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Patch
import os
import itertools
from scipy.optimize import minimize_scalar, minimize

# --- Color palette (shared with paper figures) ---
CORAL  = '#e76f51'
TEAL   = '#2a9d8f'
GOLD   = '#e9c46a'
PURPLE = '#a855f7'
BG_COLOR = '#f8f9fa'

# Output directory = images/ folder (sibling to code/)
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      'images')

# --- Physical constants (CODATA 2018) ---
HBAR  = 1.054571817e-34       # J*s
C_SI  = 2.99792458e8          # m/s
G_N   = 6.67430e-11           # m^3/(kg*s^2)
M_E   = 9.1093837015e-31      # kg  (electron mass)
M_MU  = 1.8835316e-28         # kg  (muon mass)
M_TAU = 3.16754e-27           # kg  (tau mass)
ALPHA = 1.0 / 137.035999084   # fine structure constant

# Derived
M_P = np.sqrt(HBAR * C_SI / G_N)  # Planck mass

# Target mass ratios
R_MU_E  = M_MU / M_E     # 206.768...
R_TAU_E = M_TAU / M_E    # 3477.48...
R_TAU_MU = M_TAU / M_MU  # 16.818...

# Battye-Sutcliffe normalized soliton energies (from numerical simulations)
# E(H) / E(1), where H is the Hopf charge
BS_ENERGIES = {
    1: 1.00,
    2: 1.48,
    3: 1.80,
    4: 2.18,
    5: 2.40,
    6: 2.76,
    7: 2.97,
    8: 3.45
}


def separator(title, char='=', width=80):
    """Print a formatted section separator."""
    print()
    print(char * width)
    print("  " + title)
    print(char * width)


def sub_separator(title, char='-', width=72):
    """Print a sub-section separator."""
    print()
    print("  " + char * width)
    print("  " + title)
    print("  " + char * width)


# =====================================================================
# APPROACH A: Hopf Charge Scaling
# =====================================================================

def approach_a_hopf_scaling():
    """
    Search for mass relationships m(H) = m_e * f(H).

    Sub-approaches:
      A1: f(H) = H^a  (power law, fit exponent a)
      A2: f(H) = E_BS(H) ^ a  (Battye-Sutcliffe energies raised to power)
      A3: Fit general exponent in H^a, find H_mu, H_tau
    """
    separator("APPROACH A: Hopf Charge Scaling")
    print()
    print("  Target ratios:")
    print("    m_mu/m_e  = %.3f" % R_MU_E)
    print("    m_tau/m_e = %.3f" % R_TAU_E)
    print("    m_tau/m_mu = %.3f" % R_TAU_MU)

    # ---- A1: Simple power law m(H)/m(1) = H^a ----
    sub_separator("A1: Power Law  m(H)/m_e = H^a")

    # For muon: H_mu^a = 206.768 => a = ln(206.768)/ln(H_mu)
    # For tau:  H_tau^a = 3477.48 => a = ln(3477.48)/ln(H_tau)
    # Consistency: tau/mu ratio = (H_tau/H_mu)^a = 16.818

    best_a1 = []
    print()
    print("  %4s  %5s  %12s  %12s  %12s  %12s  %12s" % (
        "H_mu", "H_tau", "a(from mu)", "a(from tau)",
        "pred mu/e", "pred tau/e", "max_err(%)"))
    print("  " + "-" * 72)

    for H_mu in range(2, 200):
        for H_tau in range(H_mu + 1, 500):
            a_mu = np.log(R_MU_E) / np.log(H_mu)
            a_tau = np.log(R_TAU_E) / np.log(H_tau)

            # Find best compromise a
            def cost_a(a):
                pred_mu = H_mu ** a
                pred_tau = H_tau ** a
                err_mu = (pred_mu / R_MU_E - 1) ** 2
                err_tau = (pred_tau / R_TAU_E - 1) ** 2
                return err_mu + err_tau

            res = minimize_scalar(cost_a, bounds=(0.1, 20), method='bounded')
            a_best = res.x
            pred_mu = H_mu ** a_best
            pred_tau = H_tau ** a_best
            err_mu = abs(pred_mu / R_MU_E - 1) * 100
            err_tau = abs(pred_tau / R_TAU_E - 1) * 100
            max_err = max(err_mu, err_tau)

            if max_err < 1.0:
                best_a1.append((max_err, H_mu, H_tau, a_best,
                                pred_mu, pred_tau, err_mu, err_tau))

    best_a1.sort()
    for item in best_a1[:20]:
        max_err, H_mu, H_tau, a, pred_mu, pred_tau, err_mu, err_tau = item
        print("  %4d  %5d  %12.6f  %12.6f  %12.3f  %12.3f  %12.4f" % (
            H_mu, H_tau, np.log(R_MU_E)/np.log(H_mu),
            np.log(R_TAU_E)/np.log(H_tau),
            pred_mu, pred_tau, max_err))

    if not best_a1:
        print("  No solutions with <1% error found for H_mu < 200, H_tau < 500")

    # ---- A2: Battye-Sutcliffe energy scaling ----
    sub_separator("A2: Battye-Sutcliffe Energy Scaling  m(H)/m_e = E_BS(H)^a")

    # Use known BS energies to see if any pair of H values works
    H_vals = sorted(BS_ENERGIES.keys())
    best_a2 = []
    print()
    print("  %4s  %5s  %8s  %8s  %12s  %12s  %12s  %12s" % (
        "H_mu", "H_tau", "E_mu", "E_tau", "a_best",
        "pred mu/e", "pred tau/e", "max_err(%)"))
    print("  " + "-" * 80)

    for H_mu in H_vals:
        for H_tau in H_vals:
            if H_tau <= H_mu:
                continue
            E_mu = BS_ENERGIES[H_mu]
            E_tau = BS_ENERGIES[H_tau]
            if E_mu <= 1.0:
                continue  # electron is E(1) = 1

            # Find a such that E_mu^a ~ R_MU_E and E_tau^a ~ R_TAU_E
            def cost_bs(a):
                pred_mu = E_mu ** a
                pred_tau = E_tau ** a
                err_mu = (pred_mu / R_MU_E - 1) ** 2
                err_tau = (pred_tau / R_TAU_E - 1) ** 2
                return err_mu + err_tau

            res = minimize_scalar(cost_bs, bounds=(1, 100), method='bounded')
            a_best = res.x
            pred_mu = E_mu ** a_best
            pred_tau = E_tau ** a_best
            err_mu = abs(pred_mu / R_MU_E - 1) * 100
            err_tau = abs(pred_tau / R_TAU_E - 1) * 100
            max_err = max(err_mu, err_tau)

            best_a2.append((max_err, H_mu, H_tau, E_mu, E_tau,
                            a_best, pred_mu, pred_tau))

    best_a2.sort()
    for item in best_a2[:10]:
        max_err, H_mu, H_tau, E_mu, E_tau, a, pred_mu, pred_tau = item
        err_mu = abs(pred_mu / R_MU_E - 1) * 100
        err_tau = abs(pred_tau / R_TAU_E - 1) * 100
        print("  %4d  %5d  %8.2f  %8.2f  %12.4f  %12.3f  %12.3f  %12.4f" % (
            H_mu, H_tau, E_mu, E_tau, a, pred_mu, pred_tau, max_err))

    # ---- A3: Continuous H with best-fit exponent ----
    sub_separator("A3: Continuous Hopf Charge  m(H)/m_e = H^a, solve for (H_mu, H_tau, a)")

    # The tau/mu ratio constrains: (H_tau/H_mu)^a = 16.818
    # The mu/e ratio constrains: H_mu^a = 206.768
    # So: H_tau = H_mu * 16.818^(1/a) and H_mu = 206.768^(1/a)
    print()
    print("  %12s  %12s  %12s  %12s  %12s" % (
        "a", "H_mu", "H_tau", "H_tau/H_mu", "ratio check"))
    print("  " + "-" * 64)

    for a_test in np.arange(0.5, 10.1, 0.5):
        H_mu = R_MU_E ** (1.0 / a_test)
        H_tau = R_TAU_E ** (1.0 / a_test)
        ratio = H_tau / H_mu
        check = ratio ** a_test  # should be 16.818
        print("  %12.4f  %12.4f  %12.4f  %12.4f  %12.4f" % (
            a_test, H_mu, H_tau, ratio, check))

    # Find the exponent where both H values are integers
    print()
    print("  Looking for exponent a where H_mu and H_tau are near-integers...")
    best_int = []
    for a_test in np.linspace(0.5, 15.0, 100000):
        H_mu = R_MU_E ** (1.0 / a_test)
        H_tau = R_TAU_E ** (1.0 / a_test)
        frac_mu = abs(H_mu - round(H_mu))
        frac_tau = abs(H_tau - round(H_tau))
        if round(H_mu) >= 2 and round(H_tau) >= 3:
            total_frac = frac_mu + frac_tau
            if total_frac < 0.15:
                best_int.append((total_frac, a_test, H_mu, H_tau,
                                 round(H_mu), round(H_tau)))

    best_int.sort()
    if best_int:
        print("  %12s  %12s  %12s  %6s  %6s  %12s" % (
            "a", "H_mu_cont", "H_tau_cont", "H_mu", "H_tau", "frac_err"))
        print("  " + "-" * 64)
        seen = set()
        count = 0
        for item in best_int:
            key = (item[4], item[5])
            if key in seen:
                continue
            seen.add(key)
            print("  %12.6f  %12.4f  %12.4f  %6d  %6d  %12.6f" % (
                item[1], item[2], item[3], item[4], item[5], item[0]))
            count += 1
            if count >= 15:
                break
    else:
        print("  No near-integer solutions found.")

    return best_a1, best_a2


# =====================================================================
# APPROACH B: Torus Knot (p, q) Classification
# =====================================================================

def approach_b_torus_knot():
    """
    Search for (p,q) torus knot quantum numbers.

    Mass formulas:
      B1: m(p,q)/m(1,1) = (p^2 + q^2)^a / 2^a
      B2: m(p,q)/m(1,1) = (p*q)^a
      B3: m(p,q)/m(1,1) = (p^2 + q^2 + c*p*q)^a / (2 + c)^a
    """
    separator("APPROACH B: Torus Knot (p,q) Classification")

    max_pq = 50
    a_vals = np.linspace(0.5, 5.0, 500)

    # ---- B1: (p^2 + q^2) formula ----
    sub_separator("B1: m(p,q)/m(1,1) = (p^2 + q^2)^a / 2^a")

    best_b1 = []
    for pm in range(1, max_pq + 1):
        for qm in range(1, max_pq + 1):
            s_mu = (pm**2 + qm**2) / 2.0
            if s_mu <= 1.0:
                continue
            for pt in range(1, max_pq + 1):
                for qt in range(1, max_pq + 1):
                    s_tau = (pt**2 + qt**2) / 2.0
                    if s_tau <= s_mu:
                        continue

                    # Solve: s_mu^a = R_MU_E => a = ln(R_MU_E)/ln(s_mu)
                    # Check: s_tau^a ~ R_TAU_E
                    a_from_mu = np.log(R_MU_E) / np.log(s_mu)
                    if a_from_mu < 0.3 or a_from_mu > 6.0:
                        continue
                    pred_tau = s_tau ** a_from_mu
                    err_tau = abs(pred_tau / R_TAU_E - 1)
                    if err_tau < 0.01:
                        pred_mu = s_mu ** a_from_mu
                        err_mu = abs(pred_mu / R_MU_E - 1)
                        max_err = max(err_mu, err_tau) * 100
                        best_b1.append((max_err, pm, qm, pt, qt,
                                        a_from_mu, pred_mu, pred_tau))

    best_b1.sort()
    print()
    if best_b1:
        print("  Solutions with <1%% error on BOTH ratios:")
        print("  %4s %4s %4s %4s  %10s  %12s  %12s  %10s" % (
            "p_mu", "q_mu", "p_tau", "q_tau", "a",
            "pred mu/e", "pred tau/e", "max_err(%)"))
        print("  " + "-" * 72)
        for item in best_b1[:30]:
            print("  %4d %4d %4d %4d  %10.5f  %12.3f  %12.3f  %10.4f" % (
                item[1], item[2], item[3], item[4],
                item[5], item[6], item[7], item[0]))
        print("  Total solutions found: %d" % len(best_b1))
    else:
        print("  No solutions with <1%% error found.")

    # ---- B2: (p*q) formula ----
    sub_separator("B2: m(p,q)/m(1,1) = (p*q)^a")

    best_b2 = []
    for pm in range(1, max_pq + 1):
        for qm in range(pm, max_pq + 1):  # symmetry: p <= q
            pq_mu = pm * qm
            if pq_mu <= 1:
                continue
            for pt in range(1, max_pq + 1):
                for qt in range(pt, max_pq + 1):
                    pq_tau = pt * qt
                    if pq_tau <= pq_mu:
                        continue

                    a_from_mu = np.log(R_MU_E) / np.log(pq_mu)
                    if a_from_mu < 0.3 or a_from_mu > 6.0:
                        continue
                    pred_tau = pq_tau ** a_from_mu
                    err_tau = abs(pred_tau / R_TAU_E - 1)
                    if err_tau < 0.01:
                        pred_mu = pq_mu ** a_from_mu
                        err_mu = abs(pred_mu / R_MU_E - 1)
                        max_err = max(err_mu, err_tau) * 100
                        best_b2.append((max_err, pm, qm, pt, qt,
                                        a_from_mu, pred_mu, pred_tau))

    best_b2.sort()
    print()
    if best_b2:
        print("  Solutions with <1%% error on BOTH ratios (p<=q by symmetry):")
        print("  %4s %4s %4s %4s  %10s  %12s  %12s  %10s" % (
            "p_mu", "q_mu", "p_tau", "q_tau", "a",
            "pred mu/e", "pred tau/e", "max_err(%)"))
        print("  " + "-" * 72)
        for item in best_b2[:30]:
            print("  %4d %4d %4d %4d  %10.5f  %12.3f  %12.3f  %10.4f" % (
                item[1], item[2], item[3], item[4],
                item[5], item[6], item[7], item[0]))
        print("  Total solutions found: %d" % len(best_b2))
    else:
        print("  No solutions with <1%% error found.")

    # ---- B3: Mixed formula with coupling c ----
    sub_separator("B3: m(p,q)/m(1,1) = (p^2 + q^2 + c*p*q)^a / (2+c)^a")

    best_b3 = []
    c_vals = np.linspace(-1.0, 4.0, 100)

    for c in c_vals:
        denom = 2.0 + c
        if denom <= 0:
            continue
        for pm in range(1, max_pq + 1):
            for qm in range(1, max_pq + 1):
                s_mu = (pm**2 + qm**2 + c*pm*qm) / denom
                if s_mu <= 1.001:
                    continue
                for pt in range(1, max_pq + 1):
                    for qt in range(1, max_pq + 1):
                        s_tau = (pt**2 + qt**2 + c*pt*qt) / denom
                        if s_tau <= s_mu:
                            continue

                        a_from_mu = np.log(R_MU_E) / np.log(s_mu)
                        if a_from_mu < 0.3 or a_from_mu > 6.0:
                            continue
                        pred_tau = s_tau ** a_from_mu
                        err_tau = abs(pred_tau / R_TAU_E - 1)
                        if err_tau < 0.005:
                            pred_mu = s_mu ** a_from_mu
                            err_mu = abs(pred_mu / R_MU_E - 1)
                            max_err = max(err_mu, err_tau) * 100
                            best_b3.append((max_err, pm, qm, pt, qt,
                                            c, a_from_mu, pred_mu, pred_tau))

    best_b3.sort()
    print()
    if best_b3:
        print("  Solutions with <0.5%% error on BOTH ratios:")
        print("  %4s %4s %4s %4s  %8s  %10s  %12s  %12s  %10s" % (
            "p_mu", "q_mu", "p_tau", "q_tau", "c", "a",
            "pred mu/e", "pred tau/e", "max_err(%)"))
        print("  " + "-" * 80)
        seen = set()
        count = 0
        for item in best_b3:
            key = (item[1], item[2], item[3], item[4])
            if key in seen:
                continue
            seen.add(key)
            print("  %4d %4d %4d %4d  %8.3f  %10.5f  %12.3f  %12.3f  %10.4f" % (
                item[1], item[2], item[3], item[4],
                item[5], item[6], item[7], item[8], item[0]))
            count += 1
            if count >= 30:
                break
        print("  Total solutions found: %d" % len(best_b3))
    else:
        print("  No solutions with <0.5%% error found.")

    return best_b1, best_b2, best_b3


# =====================================================================
# APPROACH C: Mass Formula Extension
# =====================================================================

def approach_c_mass_formula():
    """
    Using the paper's formula: m = m_P * alpha^(n1/2 - n2*alpha/4)

    Search over integer pairs (n1, n2) for masses near m_mu and m_tau.
    """
    separator("APPROACH C: Mass Formula Extension")
    print("  Formula: m = m_P * alpha^(n1/2 - n2*alpha/4)")
    print("  Known: (n1=21, n2=15) -> m_e with 0.008%% error")
    print()

    n1_max = 100
    n2_max = 100

    # Compute reference electron result
    exp_e = 21 / 2.0 - 15 * ALPHA / 4.0
    m_e_pred = M_P * ALPHA ** exp_e
    err_e = abs(m_e_pred - M_E) / M_E
    print("  Electron verification: (21, 15) -> m = %.6e kg, error = %.4f%%" % (
        m_e_pred, err_e * 100))
    print()

    # Search for muon
    muon_matches = []
    tau_matches = []
    all_masses = np.zeros((n1_max, n2_max))

    for n1 in range(1, n1_max + 1):
        for n2 in range(1, n2_max + 1):
            exp = n1 / 2.0 - n2 * ALPHA / 4.0
            m_pred = M_P * ALPHA ** exp
            all_masses[n1 - 1, n2 - 1] = m_pred

            err_mu = abs(m_pred - M_MU) / M_MU
            err_tau = abs(m_pred - M_TAU) / M_TAU

            if err_mu < 0.01:
                muon_matches.append((err_mu * 100, n1, n2, m_pred, exp))
            if err_tau < 0.01:
                tau_matches.append((err_tau * 100, n1, n2, m_pred, exp))

    muon_matches.sort()
    tau_matches.sort()

    sub_separator("Muon matches (< 1%% error)")
    if muon_matches:
        print("  %4s  %4s  %12s  %16s  %12s" % (
            "n1", "n2", "exponent", "m_pred (kg)", "error(%)"))
        print("  " + "-" * 56)
        for item in muon_matches[:20]:
            print("  %4d  %4d  %12.6f  %16.6e  %12.4f" % (
                item[1], item[2], item[4], item[3], item[0]))
        print("  Measured m_mu = %.6e kg" % M_MU)
    else:
        print("  No matches found within 1%%.")

    sub_separator("Tau matches (< 1%% error)")
    if tau_matches:
        print("  %4s  %4s  %12s  %16s  %12s" % (
            "n1", "n2", "exponent", "m_pred (kg)", "error(%)"))
        print("  " + "-" * 56)
        for item in tau_matches[:20]:
            print("  %4d  %4d  %12.6f  %16.6e  %12.4f" % (
                item[1], item[2], item[4], item[3], item[0]))
        print("  Measured m_tau = %.6e kg" % M_TAU)
    else:
        print("  No matches found within 1%%.")

    # Check if any triplet (e, mu, tau) gives consistent results
    sub_separator("Consistent triplets: electron(21,15) + muon + tau all < 1%%")
    if muon_matches and tau_matches:
        print("  Electron: (21, 15)")
        for mu_item in muon_matches[:5]:
            for tau_item in tau_matches[:5]:
                # Compute predicted ratios
                r_mu_e = mu_item[3] / m_e_pred
                r_tau_e = tau_item[3] / m_e_pred
                r_tau_mu = tau_item[3] / mu_item[3]
                err_r_mu = abs(r_mu_e / R_MU_E - 1) * 100
                err_r_tau = abs(r_tau_e / R_TAU_E - 1) * 100
                err_r_tm = abs(r_tau_mu / R_TAU_MU - 1) * 100
                print("  mu(%d,%d) tau(%d,%d): mu/e=%.1f(%.2f%%) tau/e=%.1f(%.2f%%) tau/mu=%.2f(%.2f%%)" % (
                    mu_item[1], mu_item[2], tau_item[1], tau_item[2],
                    r_mu_e, err_r_mu, r_tau_e, err_r_tau,
                    r_tau_mu, err_r_tm))
    else:
        print("  Insufficient matches to form triplets.")

    return muon_matches, tau_matches, all_masses


# =====================================================================
# APPROACH D: Composite Topological Quantum Numbers
# =====================================================================

def approach_d_composite():
    """
    Composite formula: m(H, p, q, n) = m_e * H^a * p^b * q^c * n^d

    Search over small quantum numbers (H, p, q, n) and fit exponents.
    """
    separator("APPROACH D: Composite Topological Quantum Numbers")
    print("  Formula: m/m_e = H^a * p^b * q^c * n^d")
    print("  Electron: (H=1, p=1, q=1, n=1) -> m/m_e = 1 (trivially)")
    print()

    # Define small quantum number ranges
    H_range = range(1, 10)
    p_range = range(1, 10)
    q_range = range(1, 10)
    n_range = range(1, 6)

    # For each (mu_params, tau_params) combination, find optimal exponents
    best_composite = []

    # Pre-compute log of quantum numbers for muon and tau candidates
    # We need: a*ln(H_mu) + b*ln(p_mu) + c*ln(q_mu) + d*ln(n_mu) = ln(R_MU_E)
    # and:     a*ln(H_tau) + b*ln(p_tau) + c*ln(q_tau) + d*ln(n_tau) = ln(R_TAU_E)

    ln_r_mu = np.log(R_MU_E)
    ln_r_tau = np.log(R_TAU_E)

    # For efficiency, try a simpler version first:
    # Fix n=1 for all, so d drops out. Fix q=p (torus symmetry), so b+c -> b.
    # Then: m/m_e = H^a * p^(2b) where the '2b' is combined.
    # Two equations, two unknowns (a, beta=2b):
    # a*ln(H_mu) + beta*ln(p_mu) = ln(R_MU_E)
    # a*ln(H_tau) + beta*ln(p_tau) = ln(R_TAU_E)

    sub_separator("D1: Simplified  m/m_e = H^a * (p*q)^b, n=1, p=q")

    best_d1 = []
    for H_mu in range(1, 9):
        for H_tau in range(H_mu, 9):
            for p_mu in range(1, 15):
                for p_tau in range(p_mu, 15):
                    # Build 2x2 linear system for (a, beta)
                    # [ln(H_mu)  ln(p_mu^2)] [a   ]   [ln(R_MU_E) ]
                    # [ln(H_tau) ln(p_tau^2)] [beta] = [ln(R_TAU_E)]
                    lnH_mu = np.log(max(H_mu, 1.001))
                    lnH_tau = np.log(max(H_tau, 1.001))
                    lnP_mu = np.log(max(p_mu**2, 1.001))
                    lnP_tau = np.log(max(p_tau**2, 1.001))

                    A_mat = np.array([[lnH_mu, lnP_mu],
                                      [lnH_tau, lnP_tau]])
                    b_vec = np.array([ln_r_mu, ln_r_tau])

                    det = A_mat[0, 0] * A_mat[1, 1] - A_mat[0, 1] * A_mat[1, 0]
                    if abs(det) < 1e-10:
                        continue

                    try:
                        x = np.linalg.solve(A_mat, b_vec)
                    except np.linalg.LinAlgError:
                        continue

                    a_fit, beta_fit = x
                    if a_fit < 0 or beta_fit < 0:
                        continue
                    if a_fit > 20 or beta_fit > 20:
                        continue

                    # Verify
                    pred_mu = (H_mu ** a_fit) * ((p_mu**2) ** beta_fit)
                    pred_tau = (H_tau ** a_fit) * ((p_tau**2) ** beta_fit)
                    err_mu = abs(pred_mu / R_MU_E - 1) * 100
                    err_tau = abs(pred_tau / R_TAU_E - 1) * 100
                    max_err = max(err_mu, err_tau)

                    if max_err < 0.01:  # effectively exact (numerical precision)
                        best_d1.append((max_err, H_mu, p_mu, H_tau, p_tau,
                                        a_fit, beta_fit, pred_mu, pred_tau))

    best_d1.sort(key=lambda x: (abs(x[5] - round(x[5])) +
                                 abs(x[6] - round(x[6]))))  # sort by "neatness"
    print()
    if best_d1:
        print("  Exact solutions (2 eqns, 2 unknowns):")
        print("  %4s %4s %5s %5s  %10s  %10s  %12s  %12s" % (
            "H_mu", "p_mu", "H_tau", "p_tau", "a", "beta",
            "pred mu/e", "pred tau/e"))
        print("  " + "-" * 72)
        seen = set()
        count = 0
        for item in best_d1:
            key = (item[1], item[2], item[3], item[4])
            if key in seen:
                continue
            seen.add(key)
            print("  %4d %4d %5d %5d  %10.5f  %10.5f  %12.3f  %12.3f" % (
                item[1], item[2], item[3], item[4],
                item[5], item[6], item[7], item[8]))
            count += 1
            if count >= 25:
                break
        print("  Total exact solutions: %d (showing most 'natural' exponents)" % len(best_d1))
    else:
        print("  No exact solutions found in range.")

    # ---- D2: Full 4-parameter search with fixed small exponents ----
    sub_separator("D2: Fixed rational exponents  m/m_e = H^a * p^b * q^c * n^d")

    best_d2 = []
    # Try simple rational exponents
    exp_candidates = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

    for a in exp_candidates:
        for b in exp_candidates:
            for c in exp_candidates:
                for d in exp_candidates:
                    # Search quantum number assignments
                    for Hm in range(1, 6):
                        for pm in range(1, 8):
                            for qm in range(1, 8):
                                for nm in range(1, 5):
                                    pred_mu = (Hm**a) * (pm**b) * (qm**c) * (nm**d)
                                    err_mu = abs(pred_mu / R_MU_E - 1)
                                    if err_mu > 0.02:
                                        continue
                                    for Ht in range(1, 6):
                                        for pt in range(1, 8):
                                            for qt in range(1, 8):
                                                for nt in range(1, 5):
                                                    pred_tau = (Ht**a) * (pt**b) * (qt**c) * (nt**d)
                                                    err_tau = abs(pred_tau / R_TAU_E - 1)
                                                    if err_tau > 0.02:
                                                        continue
                                                    max_err = max(err_mu, err_tau) * 100
                                                    if max_err < 2.0:
                                                        best_d2.append((
                                                            max_err, a, b, c, d,
                                                            Hm, pm, qm, nm,
                                                            Ht, pt, qt, nt,
                                                            pred_mu, pred_tau))

    best_d2.sort()
    print()
    if best_d2:
        print("  Best solutions with fixed rational exponents (<2%% error):")
        print("  %5s %5s %5s %5s | mu:(%s %s %s %s) tau:(%s %s %s %s) | %10s %10s %8s" % (
            "a", "b", "c", "d",
            "H", "p", "q", "n", "H", "p", "q", "n",
            "pred_mu/e", "pred_tau/e", "err(%)"))
        print("  " + "-" * 96)
        seen = set()
        count = 0
        for item in best_d2:
            key = item[1:13]
            if key in seen:
                continue
            seen.add(key)
            print("  %5.2f %5.2f %5.2f %5.2f | mu:(%d %d %d %d) tau:(%d %d %d %d) | %10.2f %10.2f %8.3f" % (
                item[1], item[2], item[3], item[4],
                item[5], item[6], item[7], item[8],
                item[9], item[10], item[11], item[12],
                item[13], item[14], item[0]))
            count += 1
            if count >= 20:
                break
        print("  Total solutions found: %d" % len(best_d2))
    else:
        print("  No solutions with <2%% error found.")

    return best_d1, best_d2


# =====================================================================
# APPROACH E: Kaluza-Klein / Compactification Tower
# =====================================================================

def approach_e_kaluza_klein():
    """
    Kaluza-Klein tower: m_n^2 = m_0^2 + n^2 / R^2

    Search for compact radius R and mode numbers n_mu, n_tau
    that reproduce the lepton mass spectrum.
    """
    separator("APPROACH E: Kaluza-Klein / Compactification Tower")
    print("  Formula: m_n^2 = m_0^2 + n^2 / R^2")
    print("  m_0 = m_e (ground state = electron)")
    print()

    # m_n / m_e = sqrt(1 + n^2 / (m_e * R)^2)
    # For muon: R_MU_E = sqrt(1 + n_mu^2 / (m_e*R)^2)
    # So: (m_e*R)^2 = n_mu^2 / (R_MU_E^2 - 1)

    # From muon: (m_e*R)^2 = n_mu^2 / (R_MU_E^2 - 1)
    # From tau:  (m_e*R)^2 = n_tau^2 / (R_TAU_E^2 - 1)
    # Consistency: n_mu^2 / (R_MU_E^2 - 1) = n_tau^2 / (R_TAU_E^2 - 1)
    # So: (n_tau/n_mu)^2 = (R_TAU_E^2 - 1) / (R_MU_E^2 - 1)

    ratio_sq = (R_TAU_E**2 - 1) / (R_MU_E**2 - 1)
    ratio = np.sqrt(ratio_sq)
    print("  Required n_tau/n_mu = sqrt[(m_tau/m_e)^2 - 1] / sqrt[(m_mu/m_e)^2 - 1]")
    print("                     = sqrt[%.1f] / sqrt[%.1f]" % (
        R_TAU_E**2 - 1, R_MU_E**2 - 1))
    print("                     = %.6f / %.6f" % (
        np.sqrt(R_TAU_E**2 - 1), np.sqrt(R_MU_E**2 - 1)))
    print("                     = %.6f" % ratio)
    print()

    # Search for integer (n_mu, n_tau) pairs
    sub_separator("E1: Integer KK modes")

    best_kk = []
    for n_mu in range(1, 200):
        n_tau_exact = n_mu * ratio
        n_tau = round(n_tau_exact)
        if n_tau < n_mu + 1:
            continue

        # Compute R from n_mu
        meR_sq = n_mu**2 / (R_MU_E**2 - 1)
        meR = np.sqrt(meR_sq)

        # Predict masses
        pred_mu = np.sqrt(1 + n_mu**2 / meR_sq)
        pred_tau = np.sqrt(1 + n_tau**2 / meR_sq)

        err_mu = abs(pred_mu / R_MU_E - 1) * 100
        err_tau = abs(pred_tau / R_TAU_E - 1) * 100
        max_err = max(err_mu, err_tau)

        # R in natural units: R = meR / m_e
        R_nat = meR / M_E  # in 1/kg units, convert to length
        # R in meters: R = hbar / (m_e * c) * meR (Compton wavelength units)
        R_meters = (HBAR / (M_E * C_SI)) * meR

        best_kk.append((max_err, n_mu, n_tau, meR, R_meters,
                         pred_mu, pred_tau, err_mu, err_tau))

    best_kk.sort()
    print()
    print("  %4s  %5s  %12s  %14s  %12s  %12s  %10s" % (
        "n_mu", "n_tau", "m_e*R", "R (meters)", "pred mu/e", "pred tau/e", "max_err(%)"))
    print("  " + "-" * 76)
    for item in best_kk[:20]:
        print("  %4d  %5d  %12.4f  %14.4e  %12.3f  %12.3f  %10.4f" % (
            item[1], item[2], item[3], item[4], item[5], item[6], item[0]))

    # ---- E2: Generalized KK with angular momentum ----
    sub_separator("E2: KK with angular momentum  m_n^2 = m_0^2 + (n^2 + l(l+1)) / R^2")

    best_kk2 = []
    for l_mu in range(0, 5):
        for l_tau in range(0, 5):
            for n_mu in range(1, 100):
                # Effective quantum number squared
                eff_mu_sq = n_mu**2 + l_mu * (l_mu + 1)
                meR_sq = eff_mu_sq / (R_MU_E**2 - 1)

                for n_tau in range(n_mu, 200):
                    eff_tau_sq = n_tau**2 + l_tau * (l_tau + 1)
                    pred_tau = np.sqrt(1 + eff_tau_sq / meR_sq)
                    err_tau = abs(pred_tau / R_TAU_E - 1) * 100

                    if err_tau < 0.5:
                        pred_mu = np.sqrt(1 + eff_mu_sq / meR_sq)
                        err_mu = abs(pred_mu / R_MU_E - 1) * 100
                        max_err = max(err_mu, err_tau)
                        if max_err < 0.5:
                            best_kk2.append((max_err, n_mu, l_mu, n_tau, l_tau,
                                             np.sqrt(meR_sq), pred_mu, pred_tau))

    best_kk2.sort()
    print()
    if best_kk2:
        print("  %4s %4s %5s %5s  %10s  %12s  %12s  %10s" % (
            "n_mu", "l_mu", "n_tau", "l_tau", "m_e*R",
            "pred mu/e", "pred tau/e", "max_err(%)"))
        print("  " + "-" * 72)
        for item in best_kk2[:20]:
            print("  %4d %4d %5d %5d  %10.4f  %12.3f  %12.3f  %10.4f" % (
                item[1], item[2], item[3], item[4],
                item[5], item[6], item[7], item[0]))
    else:
        print("  No solutions with <0.5%% error found.")

    return best_kk, best_kk2


# =====================================================================
# SUMMARY AND COMPARISON
# =====================================================================

def summarize_results(results_a1, results_a2, results_b1, results_b2, results_b3,
                      results_c_mu, results_c_tau,
                      results_d1, results_d2,
                      results_e1, results_e2):
    """Summarize the best result from each approach."""
    separator("GRAND SUMMARY: Best Result from Each Approach")

    summary = []

    # Approach A1
    if results_a1:
        item = results_a1[0]
        pred_mu = item[4]
        pred_tau = item[5]
        err_mu = abs(pred_mu / R_MU_E - 1) * 100
        err_tau = abs(pred_tau / R_TAU_E - 1) * 100
        summary.append(("A1: H^a power law",
                         "H_mu=%d, H_tau=%d, a=%.4f" % (item[1], item[2], item[3]),
                         pred_mu, pred_tau, max(err_mu, err_tau)))

    # Approach A2
    if results_a2:
        item = results_a2[0]
        err_mu = abs(item[6] / R_MU_E - 1) * 100
        err_tau = abs(item[7] / R_TAU_E - 1) * 100
        summary.append(("A2: BS energy scaling",
                         "H_mu=%d, H_tau=%d, a=%.4f" % (item[1], item[2], item[5]),
                         item[6], item[7], max(err_mu, err_tau)))

    # Approach B1
    if results_b1:
        item = results_b1[0]
        err_mu = abs(item[6] / R_MU_E - 1) * 100
        err_tau = abs(item[7] / R_TAU_E - 1) * 100
        summary.append(("B1: (p^2+q^2)^a knot",
                         "mu(%d,%d) tau(%d,%d) a=%.4f" % (
                             item[1], item[2], item[3], item[4], item[5]),
                         item[6], item[7], max(err_mu, err_tau)))

    # Approach B2
    if results_b2:
        item = results_b2[0]
        err_mu = abs(item[6] / R_MU_E - 1) * 100
        err_tau = abs(item[7] / R_TAU_E - 1) * 100
        summary.append(("B2: (p*q)^a knot",
                         "mu(%d,%d) tau(%d,%d) a=%.4f" % (
                             item[1], item[2], item[3], item[4], item[5]),
                         item[6], item[7], max(err_mu, err_tau)))

    # Approach B3
    if results_b3:
        item = results_b3[0]
        err_mu = abs(item[7] / R_MU_E - 1) * 100
        err_tau = abs(item[8] / R_TAU_E - 1) * 100
        summary.append(("B3: mixed knot formula",
                         "mu(%d,%d) tau(%d,%d) c=%.3f a=%.4f" % (
                             item[1], item[2], item[3], item[4], item[5], item[6]),
                         item[7], item[8], max(err_mu, err_tau)))

    # Approach C
    if results_c_mu and results_c_tau:
        mu_item = results_c_mu[0]
        tau_item = results_c_tau[0]
        # Recompute electron prediction
        exp_e = 21 / 2.0 - 15 * ALPHA / 4.0
        m_e_pred = M_P * ALPHA ** exp_e
        r_mu = mu_item[3] / m_e_pred
        r_tau = tau_item[3] / m_e_pred
        err_mu = abs(r_mu / R_MU_E - 1) * 100
        err_tau = abs(r_tau / R_TAU_E - 1) * 100
        summary.append(("C: mass formula ext.",
                         "mu(%d,%d) tau(%d,%d)" % (
                             mu_item[1], mu_item[2], tau_item[1], tau_item[2]),
                         r_mu, r_tau, max(err_mu, err_tau)))

    # Approach D2
    if results_d2:
        item = results_d2[0]
        summary.append(("D2: composite QN",
                         "a=%.1f b=%.1f c=%.1f d=%.1f" % (
                             item[1], item[2], item[3], item[4]),
                         item[13], item[14], item[0]))

    # Approach E1
    if results_e1:
        item = results_e1[0]
        summary.append(("E1: KK tower",
                         "n_mu=%d, n_tau=%d, R=%.2e m" % (
                             item[1], item[2], item[4]),
                         item[5], item[6], item[0]))

    # Approach E2
    if results_e2:
        item = results_e2[0]
        summary.append(("E2: KK + angular mom.",
                         "mu(n=%d,l=%d) tau(n=%d,l=%d)" % (
                             item[1], item[2], item[3], item[4]),
                         item[6], item[7], item[0]))

    # Print table
    print()
    print("  %-24s  %-40s  %10s  %10s  %10s" % (
        "Approach", "Parameters", "mu/e", "tau/e", "max_err(%)"))
    print("  " + "-" * 100)
    print("  %-24s  %-40s  %10.3f  %10.3f  %10s" % (
        "EXPERIMENTAL", "", R_MU_E, R_TAU_E, "---"))
    for name, params, pred_mu, pred_tau, max_err in summary:
        print("  %-24s  %-40s  %10.3f  %10.3f  %10.4f" % (
            name, params, pred_mu, pred_tau, max_err))

    # Assessment
    print()
    sub_separator("ASSESSMENT")
    sub1_found = False
    for name, params, pred_mu, pred_tau, max_err in summary:
        if max_err < 1.0:
            if not sub1_found:
                print("  Approaches achieving <1%% simultaneous fit:")
                sub1_found = True
            print("    %s: max error = %.4f%%" % (name, max_err))

    if not sub1_found:
        print("  NO approach achieves <1%% simultaneous fit on both ratios.")

    print()
    print("  Physical interpretation:")
    print("  - The lepton mass ratios are EXTREMELY specific numbers.")
    print("  - Simple topological formulas with few parameters can match")
    print("    one ratio (mu/e OR tau/e) but not both simultaneously.")
    print("  - Approaches with more free parameters (B, D) can fit both")
    print("    but at the cost of naturalness.")
    print("  - The KK tower (E) gives a clean physical picture but the")
    print("    required n_tau/n_mu ratio is irrational (%.4f)." % (
        np.sqrt((R_TAU_E**2 - 1) / (R_MU_E**2 - 1))))
    print("  - This is consistent with the paper's conclusion (Section 16)")
    print("    that all four generation mechanisms fail quantitatively.")

    return summary


# =====================================================================
# FIGURE GENERATION
# =====================================================================

def plot_hopf_scaling(results_a1, outdir):
    """Figure 1: m(H)/m_e vs H for best-fit scaling."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG_COLOR)

    # Left panel: continuous H^a curves for different a values
    ax = axes[0]
    ax.set_facecolor('white')

    H_cont = np.linspace(1, 500, 1000)
    a_values = [0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
    colors_list = [TEAL, CORAL, GOLD, PURPLE, '#264653', '#9b2226']

    for a_val, col in zip(a_values, colors_list):
        ax.plot(H_cont, H_cont**a_val, '-', color=col, alpha=0.6,
                label='$a = %.2f$' % a_val, linewidth=1.5)

    # Mark target ratios
    ax.axhline(y=R_MU_E, color='black', linestyle='--', alpha=0.5, linewidth=1)
    ax.axhline(y=R_TAU_E, color='black', linestyle='--', alpha=0.5, linewidth=1)
    ax.text(480, R_MU_E * 1.1, '$m_\\mu/m_e$', fontsize=10, ha='right',
            va='bottom', color='black')
    ax.text(480, R_TAU_E * 1.1, '$m_\\tau/m_e$', fontsize=10, ha='right',
            va='bottom', color='black')

    ax.set_xlabel('Hopf Charge $H$', fontsize=12)
    ax.set_ylabel('$m(H) / m_e$', fontsize=12)
    ax.set_title('Power Law: $m(H)/m_e = H^a$', fontsize=13, fontweight='bold')
    ax.set_yscale('log')
    ax.set_xlim(1, 500)
    ax.set_ylim(1, 1e5)
    ax.legend(fontsize=9, loc='upper left')
    ax.grid(True, alpha=0.3)

    # Right panel: best-fit results
    ax = axes[1]
    ax.set_facecolor('white')

    if results_a1:
        # Show top 10 solutions as points
        for i, item in enumerate(results_a1[:10]):
            max_err, H_mu, H_tau, a_best, pred_mu, pred_tau, err_mu, err_tau = item
            H_plot = np.array([1, H_mu, H_tau])
            m_plot = np.array([1.0, pred_mu, pred_tau])
            ax.plot(H_plot, m_plot, 'o-', alpha=0.5, markersize=6,
                    label='$H_\\mu=%d, H_\\tau=%d$ (%.2f%%)' % (
                        H_mu, H_tau, max_err) if i < 5 else None)

        # Mark actual ratios
        ax.axhline(y=R_MU_E, color='black', linestyle='--', alpha=0.4)
        ax.axhline(y=R_TAU_E, color='black', linestyle='--', alpha=0.4)

    ax.set_xlabel('Hopf Charge $H$', fontsize=12)
    ax.set_ylabel('$m(H) / m_e$', fontsize=12)
    ax.set_title('Best-Fit $(H_\\mu, H_\\tau, a)$', fontsize=13, fontweight='bold')
    ax.set_yscale('log')
    ax.legend(fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'lepton_mass_hopf_scaling.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


def plot_knot_search(results_b1, outdir):
    """Figure 2: Heatmap of best (p,q) pairs showing fit quality."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG_COLOR)

    max_pq = 50

    # Left panel: best error for muon (p_mu, q_mu) assignments
    # For each (p,q), find the best-possible exponent to match muon ratio
    ax = axes[0]
    ax.set_facecolor('white')

    err_grid_mu = np.full((max_pq, max_pq), np.nan)
    for p in range(1, max_pq + 1):
        for q in range(1, max_pq + 1):
            s = (p**2 + q**2) / 2.0
            if s <= 1.0:
                continue
            # Best a to match muon ratio
            a = np.log(R_MU_E) / np.log(s)
            if a < 0.1 or a > 10:
                continue
            pred = s ** a
            err = abs(pred / R_MU_E - 1)
            err_grid_mu[p-1, q-1] = err

    # For the muon panel, color by whether a is a "nice" number
    a_grid_mu = np.full((max_pq, max_pq), np.nan)
    for p in range(1, max_pq + 1):
        for q in range(1, max_pq + 1):
            s = (p**2 + q**2) / 2.0
            if s <= 1.0:
                continue
            a = np.log(R_MU_E) / np.log(s)
            if 0.1 < a < 10:
                a_grid_mu[p-1, q-1] = a

    im = ax.imshow(a_grid_mu.T, origin='lower', aspect='equal',
                   cmap='viridis', extent=[0.5, max_pq+0.5, 0.5, max_pq+0.5],
                   vmin=0.5, vmax=5.0)
    plt.colorbar(im, ax=ax, label='Exponent $a$ to match $m_\\mu/m_e$')
    ax.set_xlabel('$p$', fontsize=12)
    ax.set_ylabel('$q$', fontsize=12)
    ax.set_title('Muon: $a$ for $(p^2+q^2)^a/2^a = 206.8$',
                 fontsize=11, fontweight='bold')

    # Mark interesting contours (a = 1, 2, 3, 4)
    for a_target in [1.0, 2.0, 3.0, 4.0]:
        # s = R_MU_E^(1/a_target) => p^2 + q^2 = 2 * s
        s_target = R_MU_E ** (1.0 / a_target)
        r = np.sqrt(2 * s_target)
        theta = np.linspace(0, np.pi/2, 100)
        px = r * np.cos(theta)
        qx = r * np.sin(theta)
        mask = (px >= 0.5) & (px <= max_pq + 0.5) & (qx >= 0.5) & (qx <= max_pq + 0.5)
        if np.any(mask):
            ax.plot(px[mask], qx[mask], '--', color='white', alpha=0.7, linewidth=1)
            idx = np.argmax(mask)
            ax.text(px[idx] + 1, qx[idx] + 1, '$a=%d$' % a_target,
                    color='white', fontsize=8, alpha=0.8)

    ax.set_xlim(0.5, max_pq + 0.5)
    ax.set_ylim(0.5, max_pq + 0.5)

    # Right panel: if B1 results exist, show (p,q) pairs that give <1% for both
    ax = axes[1]
    ax.set_facecolor('white')

    if results_b1:
        # Plot mu assignments
        pm_vals = [r[1] for r in results_b1[:100]]
        qm_vals = [r[2] for r in results_b1[:100]]
        pt_vals = [r[3] for r in results_b1[:100]]
        qt_vals = [r[4] for r in results_b1[:100]]
        errs = [r[0] for r in results_b1[:100]]

        sc1 = ax.scatter(pm_vals, qm_vals, c=errs, cmap='RdYlGn_r',
                         s=40, alpha=0.7, marker='o', edgecolors='black',
                         linewidths=0.5, label='$\\mu$ assignment', vmin=0, vmax=1)
        sc2 = ax.scatter(pt_vals, qt_vals, c=errs, cmap='RdYlGn_r',
                         s=40, alpha=0.7, marker='^', edgecolors='black',
                         linewidths=0.5, label='$\\tau$ assignment', vmin=0, vmax=1)
        plt.colorbar(sc1, ax=ax, label='Max error (%)')
        ax.legend(fontsize=9)
        ax.set_title('$(p,q)$ pairs: $<1\\%$ on both ratios',
                     fontsize=11, fontweight='bold')
    else:
        ax.text(0.5, 0.5, 'No simultaneous\nsolutions found\nwith <1% error',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=14, color='gray')
        ax.set_title('$(p,q)$ knot search: no solutions',
                     fontsize=11, fontweight='bold')

    ax.set_xlabel('$p$', fontsize=12)
    ax.set_ylabel('$q$', fontsize=12)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'lepton_mass_knot_search.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


def plot_formula_scan(muon_matches, tau_matches, outdir):
    """Figure 3: Scatter plot of (n1, n2) pairs colored by lepton match."""
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    # Plot electron
    ax.plot(21, 15, 'o', markersize=14, markerfacecolor=TEAL,
            markeredgecolor='black', markeredgewidth=1.5, zorder=10,
            label='Electron (21,15): 0.008%')

    # Plot muon matches
    if muon_matches:
        n1_mu = [m[1] for m in muon_matches]
        n2_mu = [m[2] for m in muon_matches]
        errs_mu = [m[0] for m in muon_matches]
        sc_mu = ax.scatter(n1_mu, n2_mu, c=errs_mu, cmap='Oranges',
                           s=80, alpha=0.8, marker='s', edgecolors='black',
                           linewidths=0.5, zorder=5, vmin=0, vmax=1.0)
        ax.scatter([], [], c=CORAL, s=80, marker='s', edgecolors='black',
                   linewidths=0.5, label='Muon matches (<1%%)')

    # Plot tau matches
    if tau_matches:
        n1_tau = [m[1] for m in tau_matches]
        n2_tau = [m[2] for m in tau_matches]
        errs_tau = [m[0] for m in tau_matches]
        sc_tau = ax.scatter(n1_tau, n2_tau, c=errs_tau, cmap='Purples',
                            s=80, alpha=0.8, marker='^', edgecolors='black',
                            linewidths=0.5, zorder=5, vmin=0, vmax=1.0)
        ax.scatter([], [], c=PURPLE, s=80, marker='^', edgecolors='black',
                   linewidths=0.5, label='Tau matches (<1%%)')

    # Draw constant-exponent lines
    # exponent = n1/2 - n2*alpha/4
    # n2 = (n1/2 - exp) * 4/alpha
    exp_e = 21 / 2.0 - 15 * ALPHA / 4.0
    n1_line = np.linspace(1, 100, 200)
    n2_line_e = (n1_line / 2.0 - exp_e) * 4.0 / ALPHA
    ax.plot(n1_line, n2_line_e, '--', color=TEAL, alpha=0.4,
            label='Electron exponent contour')

    # If matches exist, draw contours for muon and tau exponents
    if muon_matches:
        exp_mu = muon_matches[0][4]
        n2_line_mu = (n1_line / 2.0 - exp_mu) * 4.0 / ALPHA
        ax.plot(n1_line, n2_line_mu, '--', color=CORAL, alpha=0.4,
                label='Muon exponent contour')

    if tau_matches:
        exp_tau = tau_matches[0][4]
        n2_line_tau = (n1_line / 2.0 - exp_tau) * 4.0 / ALPHA
        ax.plot(n1_line, n2_line_tau, '--', color=PURPLE, alpha=0.4,
                label='Tau exponent contour')

    ax.set_xlabel('$n_1$', fontsize=14)
    ax.set_ylabel('$n_2$', fontsize=14)
    ax.set_title('Mass Formula: $m = m_P \\cdot \\alpha^{n_1/2 - n_2\\alpha/4}$\n'
                 'Integer pairs matching lepton masses (<1%)',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, loc='upper left')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'lepton_mass_formula_scan.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


def plot_best_fits(summary, outdir):
    """Figure 4: Bar chart comparing best results from each approach."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.patch.set_facecolor(BG_COLOR)

    if not summary:
        for ax in axes:
            ax.text(0.5, 0.5, 'No results to display',
                    transform=ax.transAxes, ha='center', va='center')
        path = os.path.join(outdir, 'lepton_mass_best_fits.png')
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print("  Saved: %s" % path)
        return

    names = [s[0] for s in summary]
    pred_mu = [s[2] for s in summary]
    pred_tau = [s[3] for s in summary]
    max_errs = [s[4] for s in summary]

    x = np.arange(len(names))
    width = 0.35

    # Left panel: predicted ratios vs actual
    ax = axes[0]
    ax.set_facecolor('white')

    bars_mu = ax.bar(x - width/2, pred_mu, width, label='$m_\\mu/m_e$ (pred)',
                     color=CORAL, alpha=0.8, edgecolor='white')
    bars_tau = ax.bar(x + width/2, [p/100 for p in pred_tau], width,
                      label='$m_\\tau/m_e / 100$ (pred)',
                      color=PURPLE, alpha=0.8, edgecolor='white')

    ax.axhline(y=R_MU_E, color=CORAL, linestyle='--', alpha=0.6,
               label='Actual $m_\\mu/m_e = %.1f$' % R_MU_E)
    ax.axhline(y=R_TAU_E/100, color=PURPLE, linestyle='--', alpha=0.6,
               label='Actual $m_\\tau/m_e / 100 = %.1f$' % (R_TAU_E/100))

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=35, ha='right', fontsize=8)
    ax.set_ylabel('Mass Ratio', fontsize=12)
    ax.set_title('Predicted vs Actual Mass Ratios', fontsize=13, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')

    # Right panel: max error per approach
    ax = axes[1]
    ax.set_facecolor('white')

    colors = []
    for e in max_errs:
        if e < 0.1:
            colors.append(TEAL)
        elif e < 1.0:
            colors.append(GOLD)
        elif e < 5.0:
            colors.append(CORAL)
        else:
            colors.append('#9b2226')

    bars = ax.bar(x, max_errs, color=colors, alpha=0.85, edgecolor='white')

    for bar, err in zip(bars, max_errs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                '%.2f%%' % err, ha='center', va='bottom', fontsize=8,
                fontweight='bold')

    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='1% threshold')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=35, ha='right', fontsize=8)
    ax.set_ylabel('Maximum Error (%)', fontsize=12)
    ax.set_title('Fit Quality by Approach', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    # Color legend
    legend_elements = [
        Patch(facecolor=TEAL, label='< 0.1%'),
        Patch(facecolor=GOLD, label='0.1% - 1%'),
        Patch(facecolor=CORAL, label='1% - 5%'),
        Patch(facecolor='#9b2226', label='> 5%'),
    ]
    ax.legend(handles=legend_elements + [plt.Line2D([0], [0], color='red',
              linestyle='--', label='1% threshold')],
              fontsize=8, loc='upper right')

    plt.tight_layout()
    path = os.path.join(outdir, 'lepton_mass_best_fits.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


# =====================================================================
# MAIN
# =====================================================================

def main():
    """Main computation pipeline."""
    separator("Lepton Mass Ratio Search from Topological Quantum Numbers", '=', 80)
    print()
    print("  Physical constants:")
    print("    m_e   = %.10e kg" % M_E)
    print("    m_mu  = %.10e kg" % M_MU)
    print("    m_tau = %.10e kg" % M_TAU)
    print("    m_P   = %.10e kg" % M_P)
    print("    alpha = %.12f = 1/%.6f" % (ALPHA, 1/ALPHA))
    print()
    print("  Target mass ratios:")
    print("    m_mu/m_e   = %.3f" % R_MU_E)
    print("    m_tau/m_e  = %.3f" % R_TAU_E)
    print("    m_tau/m_mu = %.3f" % R_TAU_MU)
    print()

    # Approach A: Hopf Charge Scaling
    results_a1, results_a2 = approach_a_hopf_scaling()

    # Approach B: Torus Knot
    results_b1, results_b2, results_b3 = approach_b_torus_knot()

    # Approach C: Mass Formula Extension
    results_c_mu, results_c_tau, all_masses = approach_c_mass_formula()

    # Approach D: Composite Quantum Numbers
    results_d1, results_d2 = approach_d_composite()

    # Approach E: Kaluza-Klein Tower
    results_e1, results_e2 = approach_e_kaluza_klein()

    # Summary
    summary = summarize_results(
        results_a1, results_a2, results_b1, results_b2, results_b3,
        results_c_mu, results_c_tau, results_d1, results_d2,
        results_e1, results_e2)

    # Generate figures
    separator("Generating Figures")
    os.makedirs(OUTDIR, exist_ok=True)
    plot_hopf_scaling(results_a1, OUTDIR)
    plot_knot_search(results_b1, OUTDIR)
    plot_formula_scan(results_c_mu, results_c_tau, OUTDIR)
    plot_best_fits(summary, OUTDIR)

    separator("SEARCH COMPLETE")
    print()
    print("  All figures saved to: %s" % OUTDIR)
    print()

    return summary


if __name__ == '__main__':
    summary = main()
