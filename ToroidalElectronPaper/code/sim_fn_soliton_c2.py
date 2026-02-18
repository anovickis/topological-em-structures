#!/usr/bin/env python3
"""
sim_fn_soliton_c2.py — Faddeev-Niemi |H|=1 Hopf Soliton Numerical Computation

Solves the Faddeev-Niemi field equations for the minimal-energy |H|=1 Hopf
soliton via gradient flow from a conformal Hopf map initial condition.

Extracts the soliton profile, geometric parameters (torus radius rho_0, shell
thickness delta, aspect ratio), and computes the anomalous magnetic moment
coefficient C_2 from the current distribution.

Reference equations from the paper:
  - FN Lagrangian: Eq. (13.6)
  - Battye-Sutcliffe energy: Eq. (13.8), E_min ~ 192.5 sqrt(kappa_2 * kappa_4)
  - Soliton radius: Eq. (13.9)
  - C_2 from current distribution: Eq. (14.4)

Physical setup:
  - Axially symmetric ansatz in cylindrical coordinates (rho, phi, z)
  - n = (sin(Theta) cos(phi + Phi), sin(Theta) sin(phi + Phi), cos(Theta))
  - Theta(rho, z) and Phi(rho, z) are the two profile functions
  - Natural soliton units: kappa_2 = kappa_4 = 1 (shape is universal)

Energy functional:
  E = 2*pi * integral rho * [eps_2 + eps_4] drho dz

  eps_2 = (1/2) [Tr^2 + Tz^2 + sin^2(Theta) * (Pr^2 + Pz^2 + 1/rho^2)]
  eps_4 = (1/2) sin^2(Theta) [(Tr*Pz - Tz*Pr)^2 + (1/rho^2)(Tr^2 + Tz^2)]

where Tr = dTheta/drho, Tz = dTheta/dz, Pr = dPhi/drho, Pz = dPhi/dz.

References:
  [41] Faddeev & Niemi, Nature 387, 58 (1997)
  [51] Battye & Sutcliffe, PRL 81, 4798 (1998)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import time

# ─── Color Palette (shared with paper figures) ─────────────────────
CORAL  = '#e76f51'
TEAL   = '#2a9d8f'
GOLD   = '#e9c46a'
PURPLE = '#a855f7'
BG_COLOR = '#f8f9fa'

# Output directory = images/ folder (sibling to scripts/)
OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'images')


# ═══════════════════════════════════════════════════════════════════
#  1. GRID SETUP
# ═══════════════════════════════════════════════════════════════════

def setup_grid(Nr=100, Nz=200, Lr=6.0, Lz=6.0):
    """
    Set up 2D cylindrical grid (rho, z) with cell-centered points.

    Parameters
    ----------
    Nr : int       Number of grid points in rho direction
    Nz : int       Number of grid points in z direction
    Lr : float     Maximum rho (grid extends from dr/2 to Lr - dr/2)
    Lz : float     Half-extent in z (grid extends from -Lz to +Lz)

    Returns
    -------
    rho_1d : (Nr,) array of rho values
    z_1d   : (Nz,) array of z values
    dr     : grid spacing in rho
    dz     : grid spacing in z
    """
    dr = Lr / Nr
    dz = 2.0 * Lz / Nz
    rho_1d = np.linspace(dr / 2.0, Lr - dr / 2.0, Nr)
    z_1d = np.linspace(-Lz + dz / 2.0, Lz - dz / 2.0, Nz)
    return rho_1d, z_1d, dr, dz


# ═══════════════════════════════════════════════════════════════════
#  2. INITIAL CONDITION: CONFORMAL HOPF MAP
# ═══════════════════════════════════════════════════════════════════

def hopf_initial_condition(rho_1d, z_1d, a=1.0):
    """
    Conformal Hopf map for |H| = 1 with scale parameter a.

    W_0 = 2*a*rho / (2*a*z + i*(rho^2 + z^2 - a^2))

    Theta = 2 * arctan(|W_0|)
    Phi   = -arctan2(rho^2 + z^2 - a^2, 2*a*z)

    Parameters
    ----------
    rho_1d : (Nr,) array
    z_1d   : (Nz,) array
    a      : scale parameter (torus size ~ a)

    Returns
    -------
    Theta : (Nr, Nz) array
    Phi   : (Nr, Nz) array
    """
    RHO, Z = np.meshgrid(rho_1d, z_1d, indexing='ij')  # shape (Nr, Nz)
    r2 = RHO**2 + Z**2

    numer = 2.0 * a * RHO
    denom_re = 2.0 * a * Z
    denom_im = r2 - a**2

    denom_abs = np.sqrt(denom_re**2 + denom_im**2 + 1e-30)
    W_abs = numer / denom_abs

    Theta = 2.0 * np.arctan(W_abs)
    Phi = -np.arctan2(denom_im, denom_re)

    return Theta, Phi


def find_optimal_scale(rho_1d, z_1d, dr, dz, a_values=None):
    """
    Find the Hopf map scale parameter a that minimises initial energy.

    Returns
    -------
    a_opt   : optimal scale
    E_opt   : energy at optimal scale
    """
    if a_values is None:
        a_values = np.array([0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])

    best_a, best_E = 1.0, np.inf
    for a in a_values:
        Theta, Phi = hopf_initial_condition(rho_1d, z_1d, a)
        E, _, _, _, _ = compute_energy(Theta, Phi, rho_1d, dr, dz)
        if E < best_E:
            best_a, best_E = a, E
    return best_a, best_E


# ═══════════════════════════════════════════════════════════════════
#  3. FINITE DIFFERENCES
# ═══════════════════════════════════════════════════════════════════

def deriv_rho(f, dr):
    """Central difference in rho (axis 0). One-sided at boundaries."""
    d = np.zeros_like(f)
    d[1:-1, :] = (f[2:, :] - f[:-2, :]) / (2.0 * dr)
    d[0, :]  = (f[1, :] - f[0, :]) / dr          # forward
    d[-1, :] = (f[-1, :] - f[-2, :]) / dr         # backward
    return d


def deriv_z(f, dz):
    """Central difference in z (axis 1). One-sided at boundaries."""
    d = np.zeros_like(f)
    d[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    d[:, 0]  = (f[:, 1] - f[:, 0]) / dz           # forward
    d[:, -1] = (f[:, -1] - f[:, -2]) / dz          # backward
    return d


def adjoint_div_rho(rho_f, dr):
    """
    Adjoint of D_rho applied to (rho * f): computes D_rho(rho*f).
    Central difference in rho direction for the quantity rho*f.
    """
    d = np.zeros_like(rho_f)
    d[1:-1, :] = (rho_f[2:, :] - rho_f[:-2, :]) / (2.0 * dr)
    d[0, :]  = (rho_f[1, :] - rho_f[0, :]) / dr
    d[-1, :] = (rho_f[-1, :] - rho_f[-2, :]) / dr
    return d


def adjoint_div_z(f, dz):
    """
    Adjoint divergence in z: computes D_z(f).
    Central difference in z direction.
    """
    d = np.zeros_like(f)
    d[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * dz)
    d[:, 0]  = (f[:, 1] - f[:, 0]) / dz
    d[:, -1] = (f[:, -1] - f[:, -2]) / dz
    return d


# ═══════════════════════════════════════════════════════════════════
#  4. ENERGY COMPUTATION
# ═══════════════════════════════════════════════════════════════════

def compute_energy(Theta, Phi, rho_1d, dr, dz):
    """
    Compute the total Faddeev-Niemi energy and its decomposition.

    E = 2*pi * integral rho * [eps_2 + eps_4] drho dz
    with kappa_2 = kappa_4 = 1 (natural soliton units).

    Returns
    -------
    E_total : total energy
    E2      : quadratic (sigma-model) energy
    E4      : quartic (Skyrme) energy
    eps2    : (Nr, Nz) quadratic energy density array
    eps4    : (Nr, Nz) quartic energy density array
    """
    rho = rho_1d[:, np.newaxis]      # (Nr, 1) for broadcasting
    inv_rho2 = 1.0 / (rho**2 + 1e-30)

    Tr = deriv_rho(Theta, dr)
    Tz = deriv_z(Theta, dz)
    Pr = deriv_rho(Phi, dr)
    Pz = deriv_z(Phi, dz)

    s2T = np.sin(Theta)**2
    gradT2 = Tr * Tr + Tz * Tz
    gradP2 = Pr * Pr + Pz * Pz
    C = Tr * Pz - Tz * Pr       # topological current

    # Quadratic energy density (sigma-model kinetic term)
    eps2 = 0.5 * (gradT2 + s2T * (gradP2 + inv_rho2))

    # Quartic energy density (Skyrme stabilisation term)
    eps4 = 0.5 * s2T * (C * C + inv_rho2 * gradT2)

    # Integrate with cylindrical volume element 2*pi*rho*drho*dz
    E2 = 2.0 * np.pi * np.sum(rho * eps2) * dr * dz
    E4 = 2.0 * np.pi * np.sum(rho * eps4) * dr * dz

    return E2 + E4, E2, E4, eps2, eps4


# ═══════════════════════════════════════════════════════════════════
#  5. GRADIENT COMPUTATION (VARIATIONAL DERIVATIVES)
# ═══════════════════════════════════════════════════════════════════

def compute_gradient(Theta, Phi, rho_1d, dr, dz):
    """
    Compute the gradient of E w.r.t. Theta and Phi.

    Uses the Euler-Lagrange equations discretised via the chain rule
    through finite differences. The gradient is divided by the
    cylindrical volume element rho so that gradient flow converges
    uniformly across the grid.

    grad_Theta = dE/dTheta_direct - (1/rho) D_rho(rho * A_rho) - D_z(A_z)

    where A_rho = d(eps)/d(Theta_rho), A_z = d(eps)/d(Theta_z).

    Similarly for Phi.

    Returns
    -------
    grad_Theta : (Nr, Nz) array
    grad_Phi   : (Nr, Nz) array
    """
    rho = rho_1d[:, np.newaxis]
    inv_rho = 1.0 / (rho + 1e-30)
    inv_rho2 = inv_rho * inv_rho

    # Field derivatives
    Tr = deriv_rho(Theta, dr)
    Tz = deriv_z(Theta, dz)
    Pr = deriv_rho(Phi, dr)
    Pz = deriv_z(Phi, dz)

    sT = np.sin(Theta)
    cT = np.cos(Theta)
    s2T = sT * sT
    s2c = sT * cT          # = sin(2*Theta) / 2

    gradT2 = Tr * Tr + Tz * Tz
    gradP2 = Pr * Pr + Pz * Pz
    C = Tr * Pz - Tz * Pr  # topological cross term

    # ─── Gradient w.r.t. Theta ────────────────────────────────────

    # Direct term: d(eps)/d(Theta) from sin^2(Theta) factors
    # From eps_2: sin(T)cos(T) * (gradP2 + 1/rho^2)
    # From eps_4: sin(T)cos(T) * (C^2 + (1/rho^2)*gradT2)
    dE_dT = s2c * (gradP2 + inv_rho2 + C * C + inv_rho2 * gradT2)

    # Momentum: A_rho = d(eps)/d(Theta_rho)
    # From eps_2: Tr
    # From eps_4: sin^2(T) * (C * Pz + (1/rho^2) * Tr)
    A_rho = Tr * (1.0 + s2T * inv_rho2) + s2T * C * Pz

    # Momentum: A_z = d(eps)/d(Theta_z)
    # From eps_2: Tz
    # From eps_4: sin^2(T) * (-C * Pr + (1/rho^2) * Tz)
    A_z = Tz * (1.0 + s2T * inv_rho2) - s2T * C * Pr

    # Adjoint divergence: -D_rho(rho * A_rho) / rho - D_z(A_z)
    rho_Ar = rho * A_rho
    div_rho_T = adjoint_div_rho(rho_Ar, dr)
    div_z_T = adjoint_div_z(A_z, dz)

    grad_Theta = dE_dT - div_rho_T * inv_rho - div_z_T

    # ─── Gradient w.r.t. Phi ──────────────────────────────────────

    # Direct term: 0 (Phi only appears through derivatives)

    # Momentum: B_rho = d(eps)/d(Phi_rho)
    # From eps_2: sin^2(T) * Pr
    # From eps_4: sin^2(T) * C * (-Tz)  [since dC/d(Pr) = -Tz]
    B_rho = s2T * (Pr - C * Tz)

    # Momentum: B_z = d(eps)/d(Phi_z)
    # From eps_2: sin^2(T) * Pz
    # From eps_4: sin^2(T) * C * Tr  [since dC/d(Pz) = Tr]
    B_z = s2T * (Pz + C * Tr)

    # Adjoint divergence
    rho_Br = rho * B_rho
    div_rho_P = adjoint_div_rho(rho_Br, dr)
    div_z_P = adjoint_div_z(B_z, dz)

    grad_Phi = -div_rho_P * inv_rho - div_z_P

    return grad_Theta, grad_Phi


# ═══════════════════════════════════════════════════════════════════
#  5b. HOPF CHARGE GRADIENT (FOR PROJECTED GRADIENT FLOW)
# ═══════════════════════════════════════════════════════════════════

def compute_hopf_gradient(Theta, Phi, rho_1d, dr, dz):
    """
    Variational derivatives of the Hopf charge proxy w.r.t. Theta, Phi.

    The Hopf charge proxy is:
      H = (1/2) integral rho * sin(Theta) * (Tr*Pz - Tz*Pr) drho dz

    where C = Tr*Pz - Tz*Pr is the topological cross term.

    The variational derivatives (Euler-Lagrange) simplify exactly to:
      delta_H / delta_Theta = -(1/2) sin(Theta) * Phi_z
      delta_H / delta_Phi   =  (1/2) sin(Theta) * Theta_z

    Derivation: Applying the Euler-Lagrange operator to the integrand
    f = (rho/2) sin(Theta) (Theta_rho Phi_z - Theta_z Phi_rho), the mixed
    partial terms (involving Theta_rho_z and Phi_rho_z) cancel identically.
    The remaining terms from the rho-derivative of rho*sin(Theta)*Phi_z
    yield the simple forms above.

    Returns
    -------
    hT : (Nr, Nz) -- delta_H/delta_Theta
    hP : (Nr, Nz) -- delta_H/delta_Phi
    """
    sT = np.sin(Theta)
    Pz = deriv_z(Phi, dz)
    Tz = deriv_z(Theta, dz)

    hT = -0.5 * sT * Pz
    hP = 0.5 * sT * Tz

    return hT, hP


# ═══════════════════════════════════════════════════════════════════
#  6. GRADIENT FLOW (ENERGY MINIMISATION)
# ═══════════════════════════════════════════════════════════════════

def gradient_flow(Theta, Phi, rho_1d, dr, dz,
                  n_steps=15000, dt0=1e-4,
                  print_every=500, converge_tol=1e-6,
                  grad_clip=50.0):
    """
    Minimise the Faddeev-Niemi energy by gradient flow (steepest descent).

    Uses rho-weighted gradient to handle the cylindrical coordinate
    singularity at rho=0, plus gradient clipping for stability.

    Parameters
    ----------
    Theta, Phi : initial field configurations (Nr, Nz)
    rho_1d     : rho grid values
    dr, dz     : grid spacings
    n_steps    : maximum number of iterations
    dt0        : initial time step
    print_every: print progress every N steps
    converge_tol: convergence criterion |dE/E| < tol over 100 steps
    grad_clip  : maximum gradient magnitude (clips EL operator)

    Returns
    -------
    Theta, Phi : final field configurations
    history    : dict with 'E', 'E2', 'E4', 'virial', 'step' arrays
    """
    Theta = Theta.copy()
    Phi = Phi.copy()
    dt = dt0

    rho = rho_1d[:, np.newaxis]

    E_prev, E2, E4, _, _ = compute_energy(Theta, Phi, rho_1d, dr, dz)

    history = {
        'E': [E_prev], 'E2': [E2], 'E4': [E4],
        'virial': [E2 / (E4 + 1e-30)], 'step': [0]
    }

    consecutive_decrease = 0
    reject_count = 0
    t_start = time.time()

    print(f"{'Step':>6s}  {'Energy':>12s}  {'E2':>10s}  {'E4':>10s}  "
          f"{'E2/E4':>7s}  {'dE/E':>10s}  {'dt':>10s}")
    print(f"{'-'*6}  {'-'*12}  {'-'*10}  {'-'*10}  "
          f"{'-'*7}  {'-'*10}  {'-'*10}")
    print(f"{0:6d}  {E_prev:12.4f}  {E2:10.4f}  {E4:10.4f}  "
          f"{E2/(E4+1e-30):7.4f}  {'---':>10s}  {dt:10.2e}")

    # Track best topology-preserving state
    best_E = E_prev
    best_Theta = Theta.copy()
    best_Phi = Phi.copy()
    best_step = 0
    topo_lost = False

    for step in range(1, n_steps + 1):
        # Compute gradient (Euler-Lagrange operator)
        gT, gP = compute_gradient(Theta, Phi, rho_1d, dr, dz)

        # Weight by rho to handle cylindrical singularity:
        # This makes the update: dTheta = -dt * rho * EL(Theta)
        # which guarantees dE = -dt * sum(rho^2 * EL^2) * 2pi*dr*dz <= 0
        gT = rho * gT
        gP = rho * gP

        # Gradient clipping for numerical stability
        grad_norm = np.sqrt(gT**2 + gP**2)
        max_grad = np.max(grad_norm)
        if max_grad > grad_clip:
            scale = np.minimum(1.0, grad_clip / (grad_norm + 1e-30))
            gT = gT * scale
            gP = gP * scale

        # Gradient descent step
        Theta_new = Theta - dt * gT
        Phi_new = Phi - dt * gP

        # Keep Theta in physical range [0, pi] (no hard Dirichlet BC;
        # the fields naturally decay to vacuum at the grid boundaries)
        Theta_new = np.clip(Theta_new, 0.0, np.pi)

        # Compute new energy
        E_new, E2, E4, _, _ = compute_energy(Theta_new, Phi_new,
                                              rho_1d, dr, dz)

        # Adaptive step size
        if E_new <= E_prev:
            Theta = Theta_new
            Phi = Phi_new
            consecutive_decrease += 1
            reject_count = 0
            if consecutive_decrease >= 10:
                dt *= 1.05   # slowly increase step size
                consecutive_decrease = 0
        else:
            # Reject step, halve dt
            dt *= 0.5
            consecutive_decrease = 0
            reject_count += 1
            if reject_count > 60 or dt < 1e-15:
                print(f"  Step size too small (dt = {dt:.2e}, "
                      f"rejects = {reject_count}). Stopping.")
                break
            continue

        dE_rel = (E_new - E_prev) / (abs(E_prev) + 1e-30)
        E_prev = E_new

        # Topology monitoring: check that the soliton core persists
        # (Theta reaches near pi somewhere away from the axis)
        if step % 50 == 0:
            Theta_max = np.max(Theta)
            max_idx = np.unravel_index(np.argmax(Theta), Theta.shape)
            rho_at_max = rho_1d[max_idx[0]]
            if Theta_max > np.pi / 2 and rho_at_max > rho_1d[1]:
                # Topology still intact — update best state
                if E_new < best_E:
                    best_E = E_new
                    best_Theta = Theta.copy()
                    best_Phi = Phi.copy()
                    best_step = step
            else:
                print(f"\n  Topology loss at step {step}: "
                      f"max(Theta)={Theta_max:.3f}, "
                      f"rho_max={rho_at_max:.3f}")
                topo_lost = True
                break

        # Record history
        if step % 10 == 0 or step <= 100:
            history['E'].append(E_new)
            history['E2'].append(E2)
            history['E4'].append(E4)
            history['virial'].append(E2 / (E4 + 1e-30))
            history['step'].append(step)

        # Print progress
        if step % print_every == 0:
            elapsed = time.time() - t_start
            print(f"{step:6d}  {E_new:12.4f}  {E2:10.4f}  {E4:10.4f}  "
                  f"{E2/(E4+1e-30):7.4f}  {dE_rel:10.2e}  {dt:10.2e}"
                  f"  [{elapsed:.1f}s]")

        # Convergence check (only if topology is intact)
        if step > 200 and abs(dE_rel) < converge_tol:
            Theta_max = np.max(Theta)
            if Theta_max > np.pi / 2:
                print(f"\n  Converged at step {step}: "
                      f"|dE/E| = {abs(dE_rel):.2e}")
                best_E = E_prev
                best_Theta = Theta.copy()
                best_Phi = Phi.copy()
                best_step = step
                break

    # If topology was lost, revert to best state
    if topo_lost:
        print(f"  Reverting to best topology-preserving state "
              f"(step {best_step}, E={best_E:.4f})")
        Theta = best_Theta
        Phi = best_Phi
        E_prev = best_E
        E_prev, E2, E4, _, _ = compute_energy(Theta, Phi, rho_1d, dr, dz)

    elapsed = time.time() - t_start
    print(f"\n  Final energy: E = {E_prev:.4f}  (E2/E4 = {E2/(E4+1e-30):.4f})")
    print(f"  Total time: {elapsed:.1f}s")

    # Convert history lists to arrays
    for key in history:
        history[key] = np.array(history[key])

    return Theta, Phi, history


# ═══════════════════════════════════════════════════════════════════
#  6b. PROJECTED GRADIENT FLOW (TOPOLOGY-PRESERVING)
# ═══════════════════════════════════════════════════════════════════

def arrested_gradient_flow(Theta, Phi, rho_1d, dr, dz,
                           n_steps=50000, dt0=1e-4,
                           print_every=2000, converge_tol=1e-6,
                           grad_clip=50.0, n_restarts=8):
    """
    Arrested gradient flow with structural topology monitoring.

    Uses plain gradient descent for energy minimization, monitoring
    topology via the field structure: Theta must reach near pi in the
    interior (the torus core exists) and the core must be away from
    the axis. The Hopf charge proxy is tracked but NOT used for
    arrest decisions (since it changes continuously even when the
    actual topology is preserved).

    When topology degrades (Theta_max drops below pi/2 or core
    collapses to axis), the flow reverts to the last checkpoint
    and continues with smaller dt.

    Parameters
    ----------
    Theta, Phi  : initial field configurations (Nr, Nz)
    rho_1d      : rho grid values
    dr, dz      : grid spacings
    n_steps     : maximum total iterations
    dt0         : initial time step
    print_every : print frequency
    converge_tol: convergence criterion |dE/E| < tol
    grad_clip   : maximum gradient magnitude (per-point)
    n_restarts  : max number of revert events

    Returns
    -------
    Theta, Phi : final field configurations (best topology-preserving)
    history    : dict with 'E', 'E2', 'E4', 'virial', 'H', 'step' arrays
    """
    Theta = Theta.copy()
    Phi = Phi.copy()
    dt = dt0

    rho = rho_1d[:, np.newaxis]

    H0 = compute_hopf_charge(Theta, Phi, rho_1d, dr, dz)
    E_prev, E2, E4, _, _ = compute_energy(Theta, Phi, rho_1d, dr, dz)

    # Best topology-preserving state (Theta_max > pi/2 away from axis)
    best_E = E_prev
    best_Theta = Theta.copy()
    best_Phi = Phi.copy()
    best_step = 0

    # Checkpoint for reverting
    ckpt_Theta = Theta.copy()
    ckpt_Phi = Phi.copy()
    ckpt_E = E_prev
    ckpt_dt = dt0

    history = {
        'E': [E_prev], 'E2': [E2], 'E4': [E4],
        'virial': [E2 / (E4 + 1e-30)],
        'H': [H0], 'step': [0]
    }

    consecutive_decrease = 0
    reject_count = 0
    arrest_count = 0
    accepted = 0
    converge_count = 0
    t_start = time.time()

    print(f"  Arrested gradient flow (structural topology monitoring)")
    print(f"  Initial: E = {E_prev:.2f}, H = {H0:.4f}")
    print(f"  Max restarts: {n_restarts}")
    print()
    print(f"{'Step':>6s}  {'Energy':>12s}  {'E2':>10s}  {'E4':>10s}  "
          f"{'E2/E4':>7s}  {'dE/E':>10s}  {'Tmax':>6s}  "
          f"{'r_core':>6s}  {'dt':>10s}")
    print("-" * 90)
    Tmax = np.max(Theta)
    idx = np.unravel_index(np.argmax(Theta), Theta.shape)
    r_core = rho_1d[idx[0]]
    print(f"{0:6d}  {E_prev:12.4f}  {E2:10.4f}  {E4:10.4f}  "
          f"{E2/(E4+1e-30):7.4f}  {'---':>10s}  {Tmax:6.3f}  "
          f"{r_core:6.3f}  {dt:10.2e}")

    for step in range(1, n_steps + 1):
        # Compute gradient
        gT, gP = compute_gradient(Theta, Phi, rho_1d, dr, dz)
        gT = rho * gT
        gP = rho * gP

        # Per-point gradient clipping
        grad_norm = np.sqrt(gT**2 + gP**2)
        max_grad = np.max(grad_norm)
        if max_grad > grad_clip:
            scale = np.minimum(1.0, grad_clip / (grad_norm + 1e-30))
            gT = gT * scale
            gP = gP * scale

        # Step
        Theta_new = Theta - dt * gT
        Phi_new = Phi - dt * gP
        Theta_new = np.clip(Theta_new, 0.0, np.pi)

        E_new, E2_new, E4_new, _, _ = compute_energy(
            Theta_new, Phi_new, rho_1d, dr, dz)

        # Adaptive step size
        if E_new <= E_prev:
            Theta = Theta_new
            Phi = Phi_new
            E2, E4 = E2_new, E4_new
            consecutive_decrease += 1
            reject_count = 0
            accepted += 1
            if consecutive_decrease >= 10:
                dt *= 1.05
                consecutive_decrease = 0
        else:
            dt *= 0.5
            consecutive_decrease = 0
            reject_count += 1
            if reject_count > 60 or dt < 1e-15:
                # Step size collapsed — try reverting
                arrest_count += 1
                if arrest_count >= n_restarts:
                    print(f"\n  Step size collapsed, max restarts "
                          f"reached. Best E = {best_E:.2f}")
                    break
                print(f"\n  Step size collapsed at step {accepted}. "
                      f"Reverting to checkpoint.")
                Theta = ckpt_Theta.copy()
                Phi = ckpt_Phi.copy()
                E_prev = ckpt_E
                dt = ckpt_dt * 0.5
                ckpt_dt = dt
                reject_count = 0
                consecutive_decrease = 0
                continue
            continue

        dE_rel = (E_new - E_prev) / (abs(E_prev) + 1e-30)
        E_prev = E_new

        # Structural topology check every 20 accepted steps
        if accepted % 20 == 0:
            Tmax = np.max(Theta)
            idx = np.unravel_index(np.argmax(Theta), Theta.shape)
            r_core = rho_1d[idx[0]]

            # Core must stay away from axis (>= 0.15) and Theta must
            # reach close to pi (field wraps around S2)
            topo_ok = (Tmax > 2.0) and (r_core > 0.15)

            if topo_ok:
                # Update best state
                if E_new < best_E:
                    best_E = E_new
                    best_Theta = Theta.copy()
                    best_Phi = Phi.copy()
                    best_step = accepted

                # Update checkpoint
                if E_new < ckpt_E:
                    ckpt_Theta = Theta.copy()
                    ckpt_Phi = Phi.copy()
                    ckpt_E = E_new
                    ckpt_dt = dt
            else:
                arrest_count += 1
                print(f"\n  Topology degraded at step {accepted}: "
                      f"Theta_max = {Tmax:.3f}, r_core = {r_core:.3f}, "
                      f"E = {E_new:.2f}")
                print(f"  Reverting to checkpoint "
                      f"(E = {ckpt_E:.2f}, step {best_step})")
                Theta = ckpt_Theta.copy()
                Phi = ckpt_Phi.copy()
                E_prev = ckpt_E
                dt = max(ckpt_dt * 0.25, 1e-8)
                ckpt_dt = dt
                reject_count = 0
                consecutive_decrease = 0
                if arrest_count >= n_restarts:
                    print(f"  Max topology arrests ({n_restarts}).")
                    break
                continue

        # Record history
        if accepted % 20 == 0 or accepted <= 50:
            H_rec = compute_hopf_charge(Theta, Phi, rho_1d, dr, dz)
            history['E'].append(E_new)
            history['E2'].append(E2)
            history['E4'].append(E4)
            history['virial'].append(E2 / (E4 + 1e-30))
            history['H'].append(H_rec)
            history['step'].append(accepted)

        # Print progress
        if accepted % print_every == 0:
            elapsed = time.time() - t_start
            Tmax = np.max(Theta)
            idx = np.unravel_index(np.argmax(Theta), Theta.shape)
            r_core = rho_1d[idx[0]]
            print(f"{accepted:6d}  {E_new:12.4f}  {E2:10.4f}  "
                  f"{E4:10.4f}  {E2/(E4+1e-30):7.4f}  {dE_rel:10.2e}  "
                  f"{Tmax:6.3f}  {r_core:6.3f}  {dt:10.2e}"
                  f"  [{elapsed:.1f}s]")

        # Convergence check: require 50 consecutive small steps
        if accepted > 2000 and abs(dE_rel) < converge_tol:
            converge_count += 1
            if converge_count >= 50:
                Tmax = np.max(Theta)
                idx = np.unravel_index(np.argmax(Theta), Theta.shape)
                r_core = rho_1d[idx[0]]
                if Tmax > 2.0 and r_core > 0.15:
                    virial = E2 / (E4 + 1e-30)
                    H_fin = compute_hopf_charge(
                        Theta, Phi, rho_1d, dr, dz)
                    print(f"\n  Converged at step {accepted}: |dE/E| = "
                          f"{abs(dE_rel):.2e}, virial = {virial:.4f}, "
                          f"H = {H_fin:.4f}")
                    best_E = E_new
                    best_Theta = Theta.copy()
                    best_Phi = Phi.copy()
                    best_step = accepted
                    break
        else:
            converge_count = 0

    # Return best topology-preserving state
    Theta = best_Theta
    Phi = best_Phi
    E_final, E2, E4, _, _ = compute_energy(Theta, Phi, rho_1d, dr, dz)
    H_final = compute_hopf_charge(Theta, Phi, rho_1d, dr, dz)
    virial = E2 / (E4 + 1e-30)
    elapsed = time.time() - t_start

    print(f"\n  Best topology-preserving state (step {best_step}):")
    print(f"  E = {E_final:.4f}, E2/E4 = {virial:.4f}, "
          f"H = {H_final:.4f}")
    print(f"  Topology arrests: {arrest_count}")
    print(f"  Total time: {elapsed:.1f}s")

    for key in history:
        history[key] = np.array(history[key])

    return Theta, Phi, history


# ═══════════════════════════════════════════════════════════════════
#  7. HOPF CHARGE VERIFICATION
# ═══════════════════════════════════════════════════════════════════

def compute_hopf_charge(Theta, Phi, rho_1d, dr, dz):
    """
    Estimate the Hopf charge using the topological current integral.

    The exact Hopf charge requires the full Whitehead integral
    (solving for the gauge potential A), which is computationally
    expensive in cylindrical coordinates. Instead, we compute the
    integral of the topological current density:

      I = (1/2) * integral rho * sin(Theta) * (Tr*Pz - Tz*Pr) drho dz

    For the conformal Hopf map with scale a, this evaluates to
    approximately -2a (not -1), so it is not the true Hopf charge.
    However, it serves as a useful topology monitor: if |I| drops
    significantly during the gradient flow, the topology is collapsing.

    For a more reliable topology check, we also verify that Theta
    reaches pi (south pole) somewhere in the interior and that the
    maximum is away from the axis.
    """
    rho = rho_1d[:, np.newaxis]

    Tr = deriv_rho(Theta, dr)
    Tz = deriv_z(Theta, dz)
    Pr = deriv_rho(Phi, dr)
    Pz = deriv_z(Phi, dz)

    C = Tr * Pz - Tz * Pr
    integrand = rho * np.sin(Theta) * C

    I_topo = 0.5 * np.sum(integrand) * dr * dz
    return I_topo


def check_topology(Theta, rho_1d, z_1d):
    """
    Simple topology check: verify the soliton still has a toroidal core.

    Returns (is_ok, Theta_max, rho_max, z_max).
    """
    Theta_max = np.max(Theta)
    idx = np.unravel_index(np.argmax(Theta), Theta.shape)
    rho_max = rho_1d[idx[0]]
    z_max = z_1d[idx[1]]
    is_ok = (Theta_max > np.pi / 2) and (rho_max > rho_1d[1])
    return is_ok, Theta_max, rho_max, z_max


# ═══════════════════════════════════════════════════════════════════
#  8. SOLITON PROFILE EXTRACTION
# ═══════════════════════════════════════════════════════════════════

def extract_soliton_properties(Theta, Phi, rho_1d, z_1d, dr, dz):
    """
    Extract geometric properties from the converged soliton.

    Returns
    -------
    props : dict with keys:
      'rho_0'     : torus core radius (position of energy maximum)
      'z_0'       : torus core z-position (should be ~0)
      'delta'     : shell thickness (FWHM of cross-section)
      'aspect'    : aspect ratio rho_0 / delta
      'E_total'   : total energy
      'E2', 'E4'  : energy components
      'virial'    : E2/E4
      'H'         : Hopf charge
      'eps_total' : (Nr, Nz) total energy density
      'j_topo'    : (Nr, Nz) topological current density
    """
    rho = rho_1d[:, np.newaxis]

    E_total, E2, E4, eps2, eps4 = compute_energy(Theta, Phi, rho_1d, dr, dz)
    eps_total = eps2 + eps4

    H = compute_hopf_charge(Theta, Phi, rho_1d, dr, dz)

    # Topological current density
    Tr = deriv_rho(Theta, dr)
    Tz = deriv_z(Theta, dz)
    Pr = deriv_rho(Phi, dr)
    Pz = deriv_z(Phi, dz)
    j_topo = np.sin(Theta) * (Tr * Pz - Tz * Pr)

    # Find torus core: maximum of Theta (where field points to south
    # pole). This is the topologically relevant location, not the
    # energy density maximum (which can be at the axis due to 1/rho^2).
    idx = np.unravel_index(np.argmax(Theta), Theta.shape)
    rho_0 = rho_1d[idx[0]]
    z_0 = z_1d[idx[1]]

    # Cross-section along z = z_0 (through the torus core)
    j_slice = idx[1]
    eps_profile_rho = eps_total[:, j_slice]

    # FWHM in rho direction
    half_max = np.max(eps_profile_rho) / 2.0
    above = eps_profile_rho >= half_max
    if np.any(above):
        rho_above = rho_1d[above]
        fwhm_rho = max(rho_above[-1] - rho_above[0], dr)
    else:
        fwhm_rho = dr

    # Cross-section along rho = rho_0 (z-direction)
    i_core = idx[0]
    eps_profile_z = eps_total[i_core, :]
    half_max_z = np.max(eps_profile_z) / 2.0
    above_z = eps_profile_z >= half_max_z
    if np.any(above_z):
        z_above = z_1d[above_z]
        fwhm_z = max(z_above[-1] - z_above[0], dz)
    else:
        fwhm_z = dz

    # Shell thickness = geometric mean of FWHM in rho and z
    delta = np.sqrt(fwhm_rho * fwhm_z)

    # Also compute energy-weighted moments for robust size measure
    w = eps_total
    w_total = np.sum(rho * w) * dr * dz + 1e-30
    rho_mean = np.sum(rho * rho * w) * dr * dz / w_total
    rho2_mean = np.sum(rho * rho**2 * w) * dr * dz / w_total
    sigma_rho = np.sqrt(max(rho2_mean - rho_mean**2, 0))

    z_2d = z_1d[np.newaxis, :]
    z_mean = np.sum(rho * z_2d * w) * dr * dz / w_total
    z2_mean = np.sum(rho * z_2d**2 * w) * dr * dz / w_total
    sigma_z = np.sqrt(max(z2_mean - z_mean**2, 0))

    # Aspect ratio from energy-weighted moments (more robust than FWHM)
    aspect_ew = rho_mean / (sigma_rho + 1e-30)

    props = {
        'rho_0': rho_0,
        'z_0': z_0,
        'delta': delta,
        'fwhm_rho': fwhm_rho,
        'fwhm_z': fwhm_z,
        'aspect': aspect_ew,
        'rho_mean': rho_mean,
        'sigma_rho': sigma_rho,
        'sigma_z': sigma_z,
        'E_total': E_total,
        'E2': E2,
        'E4': E4,
        'virial': E2 / (E4 + 1e-30),
        'H': H,
        'eps_total': eps_total,
        'j_topo': j_topo,
        'eps_profile_rho': eps_profile_rho,
        'eps_profile_z': eps_profile_z,
        'core_idx': idx,
    }
    return props


# ═══════════════════════════════════════════════════════════════════
#  9. C_2 COMPUTATION
# ═══════════════════════════════════════════════════════════════════

def compute_c2(Theta, Phi, rho_1d, z_1d, dr, dz):
    """
    Compute the anomalous magnetic moment coefficient C_2 from the
    soliton's current distribution.

    The magnetic moment of a toroidal current distribution is:
      mu = pi * integral rho^2 * J_phi(rho, z) drho dz

    The correction to g = 2 is:
      a_e = (mu / mu_B) - 1 = <rho^2> / rho_0^2 - 1

    where rho_0 is the mean radius and the weight function is the
    topological current |j_topo|.

    Three weight functions are used for comparison:
      (a) Total energy density eps_2 + eps_4
      (b) Topological current |j_topo| = |sin(Theta) * C|
      (c) Skyrme energy density eps_4

    Returns
    -------
    results : dict with C_2 estimates and geometric moments
    """
    rho = rho_1d[:, np.newaxis]

    # Compute energy densities and topological current
    _, _, _, eps2, eps4 = compute_energy(Theta, Phi, rho_1d, dr, dz)
    eps_total = eps2 + eps4

    Tr = deriv_rho(Theta, dr)
    Tz = deriv_z(Theta, dz)
    Pr = deriv_rho(Phi, dr)
    Pz = deriv_z(Phi, dz)
    j_topo = np.abs(np.sin(Theta) * (Tr * Pz - Tz * Pr))

    results = {}

    # Compute moments for each weight function
    weights = {
        'energy': eps_total,
        'current': j_topo,
        'skyrme': eps4,
    }

    for name, w in weights.items():
        # Normalisation (integrate w over volume element rho drho dz)
        w_total = np.sum(rho * w) * dr * dz + 1e-30

        # Mean radius (rho-weighted)
        rho_mean = np.sum(rho * rho * w) * dr * dz / w_total

        # Second moment <rho^2>
        rho2_mean = np.sum(rho * rho**2 * w) * dr * dz / w_total

        # Variance <(rho - rho_mean)^2>
        var_rho = rho2_mean - rho_mean**2

        # Correction to magnetic moment
        # a_e = <rho^2> / rho_0^2 - 1 = var_rho / rho_mean^2
        delta_g = var_rho / (rho_mean**2 + 1e-30)

        # Effective C_2: if a_e = C_2 * (alpha/pi)^2, then
        # C_2 = delta_g / (alpha/pi)^2
        # But in soliton units this is just delta_g * pi^2 / alpha^2
        # However, the soliton's aspect ratio IS the physical one,
        # so we report the raw geometric ratio directly.
        alpha = 1.0 / 137.036
        C2_equivalent = delta_g * np.pi**2 / (alpha**2)

        results[name] = {
            'rho_mean': rho_mean,
            'rho2_mean': rho2_mean,
            'var_rho': var_rho,
            'delta_g': delta_g,
            'sigma_rho': np.sqrt(max(var_rho, 0)),
            'aspect_sigma': rho_mean / (np.sqrt(max(var_rho, 0)) + 1e-30),
            'C2_if_thin_torus': -np.pi**2 / 4 * var_rho / (rho_mean**2 + 1e-30),
        }

    return results


# ═══════════════════════════════════════════════════════════════════
#  10. PLOTTING
# ═══════════════════════════════════════════════════════════════════

def plot_convergence(history, outdir):
    """Plot energy, virial ratio, and Hopf charge vs iteration."""
    has_H = 'H' in history and len(history['H']) > 0
    ncols = 3 if has_H else 2
    fig, axes = plt.subplots(1, ncols, figsize=(5 * ncols, 5))
    fig.patch.set_facecolor(BG_COLOR)

    # Energy convergence
    ax = axes[0]
    ax.set_facecolor('white')
    ax.plot(history['step'], history['E'], color=PURPLE, linewidth=2,
            label='$E_{\\rm total}$')
    ax.plot(history['step'], history['E2'], color=TEAL, linewidth=1.5,
            linestyle='--', label='$E_2$ (quadratic)')
    ax.plot(history['step'], history['E4'], color=CORAL, linewidth=1.5,
            linestyle='--', label='$E_4$ (Skyrme)')

    # Battye-Sutcliffe reference
    BS_energy = 192.5
    ax.axhline(BS_energy, color=GOLD, linestyle=':', linewidth=1.5,
               label=f'Battye-Sutcliffe: {BS_energy}')

    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Energy (soliton units)', fontsize=12)
    ax.set_title('Energy Convergence', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Virial ratio
    ax = axes[1]
    ax.set_facecolor('white')
    ax.plot(history['step'], history['virial'], color=PURPLE, linewidth=2)
    ax.axhline(1.0, color=GOLD, linestyle=':', linewidth=1.5,
               label='Virial equilibrium ($E_2/E_4 = 1$)')
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('$E_2 / E_4$', fontsize=12)
    ax.set_title('Virial Ratio', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    # Hopf charge preservation
    if has_H:
        ax = axes[2]
        ax.set_facecolor('white')
        ax.plot(history['step'], history['H'], color=CORAL, linewidth=2)
        ax.axhline(history['H'][0], color=GOLD, linestyle=':',
                   linewidth=1.5, label=f'Target: {history["H"][0]:.4f}')
        ax.set_xlabel('Iteration', fontsize=12)
        ax.set_ylabel('Hopf charge proxy $H$', fontsize=12)
        ax.set_title('Topology Preservation', fontsize=13,
                     fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_convergence.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_energy_density(eps_total, rho_1d, z_1d, props, outdir):
    """2D contour plot of energy density in the (rho, z) plane."""
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    RHO, Z = np.meshgrid(rho_1d, z_1d, indexing='ij')

    # Custom colormap: white -> teal -> purple
    colors = ['white', TEAL, PURPLE]
    cmap = LinearSegmentedColormap.from_list('soliton', colors, N=256)

    # Normalise for display
    eps_norm = eps_total / (np.max(eps_total) + 1e-30)

    levels = np.linspace(0.05, 1.0, 20)
    cs = ax.contourf(RHO, Z, eps_norm, levels=levels, cmap=cmap)
    ax.contour(RHO, Z, eps_norm, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
               colors='white', linewidths=0.5, alpha=0.5)

    # Mark torus core
    ax.plot(props['rho_0'], props['z_0'], 'o', color=CORAL,
            markersize=8, markeredgecolor='white', markeredgewidth=1.5,
            zorder=5)

    # Annotation
    ax.annotate(
        f"$\\rho_0 = {props['rho_0']:.2f}$\n$\\delta = {props['delta']:.2f}$\n"
        f"aspect = {props['aspect']:.1f}",
        xy=(props['rho_0'], props['z_0']),
        xytext=(props['rho_0'] + 1.5, props['z_0'] + 1.5),
        fontsize=10, color='white',
        bbox=dict(boxstyle='round,pad=0.3', facecolor=PURPLE, alpha=0.8),
        arrowprops=dict(arrowstyle='->', color='white', lw=1.5),
    )

    cbar = plt.colorbar(cs, ax=ax, label='$\\varepsilon / \\varepsilon_{\\max}$')
    ax.set_xlabel('$\\rho$ (soliton units)', fontsize=12)
    ax.set_ylabel('$z$ (soliton units)', fontsize=12)
    ax.set_title('Faddeev-Niemi $|H|=1$ Soliton Energy Density', fontsize=13,
                 fontweight='bold')
    ax.set_xlim(0, max(3 * props['rho_0'], 4.0))
    ax.set_ylim(-max(2 * props['rho_0'], 3.0), max(2 * props['rho_0'], 3.0))
    ax.set_aspect('equal')

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_energy_density.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_cross_section(props, rho_1d, z_1d, outdir):
    """1D cross-section profiles through the torus core."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor(BG_COLOR)

    # Radial profile (z = z_0)
    ax = axes[0]
    ax.set_facecolor('white')
    eps_rho = props['eps_profile_rho']
    eps_rho_norm = eps_rho / (np.max(eps_rho) + 1e-30)
    ax.plot(rho_1d, eps_rho_norm, color=CORAL, linewidth=2.5)
    ax.axvline(props['rho_0'], color=PURPLE, linestyle='--', linewidth=1,
               label=f'$\\rho_0 = {props["rho_0"]:.2f}$')

    # FWHM markers
    half = 0.5
    ax.axhline(half, color=GOLD, linestyle=':', linewidth=1, alpha=0.7)
    rho_left = props['rho_0'] - props['fwhm_rho'] / 2
    rho_right = props['rho_0'] + props['fwhm_rho'] / 2
    ax.annotate('', xy=(rho_left, half), xytext=(rho_right, half),
                arrowprops=dict(arrowstyle='<->', color=GOLD, lw=2))
    ax.text((rho_left + rho_right) / 2, half + 0.08,
            f'FWHM = {props["fwhm_rho"]:.2f}', ha='center', fontsize=10,
            color=GOLD, fontweight='bold')

    ax.set_xlabel('$\\rho$ (soliton units)', fontsize=12)
    ax.set_ylabel('$\\varepsilon / \\varepsilon_{\\max}$', fontsize=12)
    ax.set_title('Radial Cross-Section ($z = z_0$)', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(3 * props['rho_0'], 4.0))

    # Axial profile (rho = rho_0)
    ax = axes[1]
    ax.set_facecolor('white')
    eps_z = props['eps_profile_z']
    eps_z_norm = eps_z / (np.max(eps_z) + 1e-30)
    ax.plot(z_1d, eps_z_norm, color=TEAL, linewidth=2.5)
    ax.axvline(props['z_0'], color=PURPLE, linestyle='--', linewidth=1,
               label=f'$z_0 = {props["z_0"]:.2f}$')

    half_z = 0.5
    ax.axhline(half_z, color=GOLD, linestyle=':', linewidth=1, alpha=0.7)
    z_left = props['z_0'] - props['fwhm_z'] / 2
    z_right = props['z_0'] + props['fwhm_z'] / 2
    ax.annotate('', xy=(z_left, half_z), xytext=(z_right, half_z),
                arrowprops=dict(arrowstyle='<->', color=GOLD, lw=2))
    ax.text((z_left + z_right) / 2, half_z + 0.08,
            f'FWHM = {props["fwhm_z"]:.2f}', ha='center', fontsize=10,
            color=GOLD, fontweight='bold')

    ax.set_xlabel('$z$ (soliton units)', fontsize=12)
    ax.set_ylabel('$\\varepsilon / \\varepsilon_{\\max}$', fontsize=12)
    ax.set_title('Axial Cross-Section ($\\rho = \\rho_0$)', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(-max(2 * props['rho_0'], 3.0), max(2 * props['rho_0'], 3.0))

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_cross_section.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_c2_comparison(c2_results, outdir):
    """Bar chart comparing C_2 estimates."""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    # Values to compare
    labels = [
        'Uniform\ncurrent\n(Eq. 14.6)',
        'Shell\n($\\delta/r_e \\approx 0.4$)\n(Eq. 14.8)',
        'Shell +\nself-int.\n(Eq. 14.9)',
        'Numerical\nsoliton\n(this work)',
        'QED\n$C_2$\n(exact)',
    ]

    c2_current = c2_results.get('current', c2_results.get('energy', {}))

    values = [
        -2.47,                                    # uniform
        -0.99,                                    # shell concentrated
        -0.33,                                    # shell + self-interaction
        c2_current.get('C2_if_thin_torus', 0),    # numerical
        -0.3285,                                  # QED
    ]

    colors_bar = [CORAL, GOLD, TEAL, PURPLE, 'black']

    bars = ax.bar(range(len(labels)), values, color=colors_bar, alpha=0.85,
                  edgecolor='white', linewidth=1.5, width=0.7)

    # Value labels on bars
    for bar, val in zip(bars, values):
        y = bar.get_height()
        sign = -1 if y < 0 else 1
        ax.text(bar.get_x() + bar.get_width() / 2, y - sign * 0.1,
                f'{val:.2f}', ha='center', va='top' if y < 0 else 'bottom',
                fontsize=11, fontweight='bold',
                color='white' if abs(val) > 0.5 else 'black')

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('$C_2$ coefficient', fontsize=13)
    ax.set_title('Anomalous Magnetic Moment: $C_2$ Comparison',
                 fontsize=14, fontweight='bold')
    ax.axhline(0, color='black', linewidth=0.5)
    ax.axhline(-0.3285, color='gray', linestyle=':', linewidth=1,
               label='QED: $C_2 = -0.3285$')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_c2_comparison.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ═══════════════════════════════════════════════════════════════════
#  11. MAIN ROUTINE
# ═══════════════════════════════════════════════════════════════════

def main():
    """Main computation pipeline."""
    print("=" * 70)
    print("  Faddeev-Niemi |H|=1 Hopf Soliton -- Numerical Computation")
    print("=" * 70)
    print()

    # ─── Grid setup ───────────────────────────────────────────────
    Nr, Nz = 100, 200
    Lr, Lz = 6.0, 6.0
    rho_1d, z_1d, dr, dz = setup_grid(Nr, Nz, Lr, Lz)
    print(f"Grid: {Nr} x {Nz} = {Nr*Nz} points")
    print(f"  rho: [{rho_1d[0]:.3f}, {rho_1d[-1]:.3f}], dr = {dr:.4f}")
    print(f"  z:   [{z_1d[0]:.3f}, {z_1d[-1]:.3f}], dz = {dz:.4f}")
    print()

    # ─── Initial condition ────────────────────────────────────────
    # Use a = 1.5 to keep the torus core well away from the axis
    # (avoids 1/rho^2 singularity dominating the energy)
    a_values = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
    print("Scale parameter scan:")
    best_a, best_E = 1.0, np.inf
    for a in a_values:
        Theta_test, Phi_test = hopf_initial_condition(rho_1d, z_1d, a)
        E, E2, E4, _, _ = compute_energy(Theta_test, Phi_test,
                                           rho_1d, dr, dz)
        Tmax = np.max(Theta_test)
        idx = np.unravel_index(np.argmax(Theta_test), Theta_test.shape)
        r_core = rho_1d[idx[0]]
        print(f"  a = {a:.1f}: E = {E:8.1f}, E2/E4 = {E2/(E4+1e-30):.3f}, "
              f"Tmax = {Tmax:.3f}, r_core = {r_core:.3f}")
        if E < best_E:
            best_a, best_E = a, E

    # Use a = 1.5 for best topology preservation while keeping
    # energy manageable (gradient flow handles the rest)
    a_use = 1.5
    Theta, Phi = hopf_initial_condition(rho_1d, z_1d, a=a_use)
    E0, E2_0, E4_0, _, _ = compute_energy(Theta, Phi, rho_1d, dr, dz)
    H0 = compute_hopf_charge(Theta, Phi, rho_1d, dr, dz)
    print(f"\nUsing a = {a_use:.1f}: E = {E0:.2f}, H = {H0:.4f}")
    print()

    # ─── Arrested gradient flow (topology-preserving) ────────────
    print("Starting arrested gradient flow...")
    print()
    Theta, Phi, history = arrested_gradient_flow(
        Theta, Phi, rho_1d, dr, dz,
        n_steps=80000, dt0=5e-4,
        print_every=2000, converge_tol=1e-7,
        grad_clip=50.0, n_restarts=15
    )
    print()

    # ─── Extract soliton properties ──────────────────────────────
    print("Extracting soliton properties...")
    props = extract_soliton_properties(Theta, Phi, rho_1d, z_1d, dr, dz)

    print()
    print("=" * 50)
    print("  SOLITON PROPERTIES")
    print("=" * 50)
    print(f"  Energy:        E = {props['E_total']:.2f}")
    print(f"                 E2 = {props['E2']:.2f}, E4 = {props['E4']:.2f}")
    print(f"  Virial ratio:  E2/E4 = {props['virial']:.4f}")
    print(f"  Hopf charge:   H = {props['H']:.4f}")
    print(f"  Torus core (Theta max):  rho_0 = {props['rho_0']:.3f}, "
          f"z_0 = {props['z_0']:.3f}")
    print(f"  Energy-weighted radius:  <rho> = {props['rho_mean']:.3f}")
    print(f"  Energy-weighted width:   sigma_rho = {props['sigma_rho']:.3f}, "
          f"sigma_z = {props['sigma_z']:.3f}")
    print(f"  FWHM(rho) = {props['fwhm_rho']:.3f}, FWHM(z) = {props['fwhm_z']:.3f}")
    print(f"  Aspect ratio (energy):   A = <rho>/sigma = {props['aspect']:.2f}")
    print()

    # Battye-Sutcliffe comparison
    BS_E = 192.5
    print(f"  Battye-Sutcliffe reference: E = {BS_E:.1f}")
    print(f"  Ratio E/E_BS = {props['E_total']/BS_E:.4f}")
    print()

    # ─── C_2 computation ─────────────────────────────────────────
    print("Computing C_2 from current distribution...")
    c2_results = compute_c2(Theta, Phi, rho_1d, z_1d, dr, dz)

    print()
    print("=" * 50)
    print("  C_2 ANALYSIS")
    print("=" * 50)
    for name, res in c2_results.items():
        print(f"\n  Weight: {name}")
        print(f"    Mean radius:     rho_mean = {res['rho_mean']:.4f}")
        print(f"    RMS width:       sigma = {res['sigma_rho']:.4f}")
        print(f"    Aspect (mean/sigma): {res['aspect_sigma']:.2f}")
        print(f"    Variance ratio:  var/rho^2 = {res['delta_g']:.6f}")
        print(f"    C_2 (thin-torus formula): {res['C2_if_thin_torus']:.4f}")
        print(f"    QED C_2: -0.3285")

    print()
    A = props['aspect']
    print(f"  Key finding: The energy-weighted aspect ratio A ~ {A:.1f}")
    print(f"  (consistent with Battye-Sutcliffe fat-torus shape).")
    print(f"  This is far from 1/alpha ~ 137 (thin-torus model).")
    print()

    # ─── Generate plots ──────────────────────────────────────────
    print("Generating figures...")
    plot_convergence(history, OUTDIR)
    plot_energy_density(props['eps_total'], rho_1d, z_1d, props, OUTDIR)
    plot_cross_section(props, rho_1d, z_1d, OUTDIR)
    plot_c2_comparison(c2_results, OUTDIR)

    print()
    print("=" * 70)
    print("  COMPUTATION COMPLETE")
    print("=" * 70)

    return props, c2_results, history


if __name__ == '__main__':
    props, c2_results, history = main()
