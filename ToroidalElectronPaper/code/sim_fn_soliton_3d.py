#!/usr/bin/env python3
"""
sim_fn_soliton_3d.py -- 3D Faddeev-Niemi |H|=1 Hopf Soliton Solver
(Arrested Damped Newton Flow on Full 3D Cartesian Grid)

Solves for the minimum-energy configuration of the Faddeev-Niemi nonlinear
sigma model with Hopf charge H=1, targeting the Battye-Sutcliffe result
E ~ 192.5 in soliton units (kappa_2 = kappa_4 = 1).

Unlike sim_fn_soliton_c2.py (which uses a 2D cylindrical ansatz), this script
works on a full 3D Cartesian grid, eliminating the 1/rho^2 axis singularity
and allowing the soliton to relax without axial-symmetry constraints.

Physics:
  The FN model has a unit-vector field n(x,y,z) on S^2, with energy
    E = integral [eps_2 + eps_4] d^3x
  where
    eps_2 = (1/2) sum_{i,a} (d_i n_a)^2         (sigma-model / Dirichlet)
    eps_4 = (1/2) sum_{i<j} [n . (d_i n x d_j n)]^2  (Skyrme stabilisation)

  The Hopf charge H = (1/4pi^2) integral A . B d^3x where
    F_ij = n . (d_i n x d_j n),  B_i = (1/2) eps_ijk F_jk,  curl A = B.

Algorithm:
  Arrested damped Newton flow (momentum-based gradient descent with energy
  monitoring -- reject steps that increase energy, reset velocity, halve dt).

References:
  [41] Faddeev & Niemi, Nature 387, 58 (1997)
  [51] Battye & Sutcliffe, PRL 81, 4798 (1998)
  [51] Battye & Sutcliffe, Proc. R. Soc. Lond. A 455, 4305 (1999)
"""

import numpy as np
try:
    from scipy.ndimage import map_coordinates as _map_coordinates
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import sys
import time
import argparse

# ---- Color palette (shared with paper figures) --------------------------
CORAL  = '#e76f51'
TEAL   = '#2a9d8f'
GOLD   = '#e9c46a'
PURPLE = '#a855f7'
BG_COLOR = '#f8f9fa'

# Output directories
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(os.path.dirname(SCRIPT_DIR), 'images')

# =========================================================================
#  1. GRID SETUP
# =========================================================================

def setup_grid(N=101, L=6.0):
    """
    Create a 3D cell-centred Cartesian grid on [-L, L]^3.

    Parameters
    ----------
    N : int     Number of grid points along each axis.
    L : float   Half-extent of the domain.

    Returns
    -------
    x1d : (N,) array of cell-centre coordinates along each axis.
    h   : float  Grid spacing.
    """
    h = 2.0 * L / N
    x1d = np.linspace(-L + h / 2.0, L - h / 2.0, N)
    return x1d, h


# =========================================================================
#  2. INITIAL CONDITION -- HOPF MAP
# =========================================================================

def hopf_map_initial(x1d, a=1.5):
    """
    Standard Hopf map R^3 -> S^2 giving Hopf charge H = 1.

    For point (x, y, z) with r^2 = x^2 + y^2 + z^2, define
      n_1 = 2(x*z + y*a) / (r^2 + a^2)
      n_2 = 2(y*z - x*a) / (r^2 + a^2)
      n_3 = (x^2 + y^2 - z^2 - a^2) / (r^2 + a^2)

    Boundary: n -> (0, 0, -1) as r -> infinity.

    Parameters
    ----------
    x1d : (N,) array of grid coordinates.
    a   : float  Scale parameter (controls torus size).

    Returns
    -------
    n : (3, N, N, N) array with |n| = 1 everywhere.
    """
    N = len(x1d)
    # Build coordinate grids -- indexing='ij' gives shape (Nx, Ny, Nz)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing='ij')
    r2 = X**2 + Y**2 + Z**2
    denom = r2 + a**2  # > 0 everywhere

    n = np.empty((3, N, N, N), dtype=np.float64)
    n[0] = 2.0 * (X * Z + Y * a) / denom
    n[1] = 2.0 * (Y * Z - X * a) / denom
    n[2] = (X**2 + Y**2 - Z**2 - a**2) / denom

    # Renormalise (should already be unit, but enforce numerically)
    norm = np.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
    norm = np.maximum(norm, 1e-30)
    n[0] /= norm
    n[1] /= norm
    n[2] /= norm

    return n


# =========================================================================
#  3. FINITE DIFFERENCE DERIVATIVES (VECTORISED)
# =========================================================================

def partial_x(f, h):
    """Central difference along axis 0, one-sided at boundaries."""
    d = np.empty_like(f)
    d[1:-1, :, :] = (f[2:, :, :] - f[:-2, :, :]) / (2.0 * h)
    d[0, :, :]    = (f[1, :, :] - f[0, :, :]) / h       # forward
    d[-1, :, :]   = (f[-1, :, :] - f[-2, :, :]) / h      # backward
    return d


def partial_y(f, h):
    """Central difference along axis 1, one-sided at boundaries."""
    d = np.empty_like(f)
    d[:, 1:-1, :] = (f[:, 2:, :] - f[:, :-2, :]) / (2.0 * h)
    d[:, 0, :]    = (f[:, 1, :] - f[:, 0, :]) / h
    d[:, -1, :]   = (f[:, -1, :] - f[:, -2, :]) / h
    return d


def partial_z(f, h):
    """Central difference along axis 2, one-sided at boundaries."""
    d = np.empty_like(f)
    d[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * h)
    d[:, :, 0]    = (f[:, :, 1] - f[:, :, 0]) / h
    d[:, :, -1]   = (f[:, :, -1] - f[:, :, -2]) / h
    return d


def laplacian(f, h):
    """7-point Laplacian (second-order central differences, 3D)."""
    lap = np.zeros_like(f)
    ih2 = 1.0 / (h * h)
    # Interior points
    lap[1:-1, :, :] += (f[2:, :, :] - 2.0 * f[1:-1, :, :] + f[:-2, :, :]) * ih2
    lap[:, 1:-1, :] += (f[:, 2:, :] - 2.0 * f[:, 1:-1, :] + f[:, :-2, :]) * ih2
    lap[:, :, 1:-1] += (f[:, :, 2:] - 2.0 * f[:, :, 1:-1] + f[:, :, :-2]) * ih2
    # Boundary: use one-sided second derivative
    # x-boundaries
    lap[0, :, :]  += (f[2, :, :] - 2.0 * f[1, :, :] + f[0, :, :]) * ih2
    lap[-1, :, :] += (f[-3, :, :] - 2.0 * f[-2, :, :] + f[-1, :, :]) * ih2
    # y-boundaries
    lap[:, 0, :]  += (f[:, 2, :] - 2.0 * f[:, 1, :] + f[:, 0, :]) * ih2
    lap[:, -1, :] += (f[:, -3, :] - 2.0 * f[:, -2, :] + f[:, -1, :]) * ih2
    # z-boundaries
    lap[:, :, 0]  += (f[:, :, 2] - 2.0 * f[:, :, 1] + f[:, :, 0]) * ih2
    lap[:, :, -1] += (f[:, :, -3] - 2.0 * f[:, :, -2] + f[:, :, -1]) * ih2
    return lap


def compute_all_derivs(n, h):
    """
    Compute all 9 first-derivative arrays d_i n_a.

    Parameters
    ----------
    n : (3, N, N, N) field array.
    h : float grid spacing.

    Returns
    -------
    dn : (3, 3, N, N, N) array where dn[a, i] = d_i n_a.
         a = field component (0,1,2), i = spatial direction (0=x,1=y,2=z).
    """
    N = n.shape[1]
    dn = np.empty((3, 3, N, N, N), dtype=np.float64)
    partials = [partial_x, partial_y, partial_z]
    for a in range(3):
        for i in range(3):
            dn[a, i] = partials[i](n[a], h)
    return dn


# =========================================================================
#  4. ENERGY COMPUTATION
# =========================================================================

def compute_Cij(n, dn):
    """
    Compute the three independent cross-product terms C_ij for i<j.

    C_ij = n . (d_i n x d_j n)
         = sum_abc eps_abc n_a (d_i n_b)(d_j n_c)

    Expanding:
      C_ij = n_0 (d_i n_1 d_j n_2 - d_i n_2 d_j n_1)
           + n_1 (d_i n_2 d_j n_0 - d_i n_0 d_j n_2)
           + n_2 (d_i n_0 d_j n_1 - d_i n_1 d_j n_0)

    Returns
    -------
    C_xy, C_xz, C_yz : each (N, N, N) arrays.
    """
    # (i,j) = (0,1) = (x,y)
    C_xy = (n[0] * (dn[1, 0] * dn[2, 1] - dn[2, 0] * dn[1, 1])
          + n[1] * (dn[2, 0] * dn[0, 1] - dn[0, 0] * dn[2, 1])
          + n[2] * (dn[0, 0] * dn[1, 1] - dn[1, 0] * dn[0, 1]))

    # (i,j) = (0,2) = (x,z)
    C_xz = (n[0] * (dn[1, 0] * dn[2, 2] - dn[2, 0] * dn[1, 2])
          + n[1] * (dn[2, 0] * dn[0, 2] - dn[0, 0] * dn[2, 2])
          + n[2] * (dn[0, 0] * dn[1, 2] - dn[1, 0] * dn[0, 2]))

    # (i,j) = (1,2) = (y,z)
    C_yz = (n[0] * (dn[1, 1] * dn[2, 2] - dn[2, 1] * dn[1, 2])
          + n[1] * (dn[2, 1] * dn[0, 2] - dn[0, 1] * dn[2, 2])
          + n[2] * (dn[0, 1] * dn[1, 2] - dn[1, 1] * dn[0, 2]))

    return C_xy, C_xz, C_yz


def compute_energy(n, h):
    """
    Compute the total Faddeev-Niemi energy and decomposition.

    E = integral [eps_2 + eps_4] d^3x,  kappa_2 = kappa_4 = 1.

    eps_2 = (1/2) sum_{i,a} (d_i n_a)^2
    eps_4 = (1/2) (C_xy^2 + C_xz^2 + C_yz^2)

    Returns
    -------
    E_total, E2, E4 : float  Total, sigma-model, and Skyrme energies.
    eps2, eps4       : (N,N,N) energy density arrays.
    dn               : (3,3,N,N,N) derivative array (reused by caller).
    """
    dn = compute_all_derivs(n, h)

    # eps_2 = (1/2) sum_{a,i} (d_i n_a)^2
    eps2 = np.zeros_like(n[0])
    for a in range(3):
        for i in range(3):
            eps2 += dn[a, i]**2
    eps2 *= 0.5

    # eps_4 = (1/2)(C_xy^2 + C_xz^2 + C_yz^2)
    C_xy, C_xz, C_yz = compute_Cij(n, dn)
    eps4 = 0.5 * (C_xy**2 + C_xz**2 + C_yz**2)

    dV = h**3
    E2 = np.sum(eps2) * dV
    E4 = np.sum(eps4) * dV
    return E2 + E4, E2, E4, eps2, eps4, dn


# =========================================================================
#  5. GRADIENT COMPUTATION (PROJECTED ONTO S^2 TANGENT PLANE)
# =========================================================================

def compute_gradient(n, h):
    """
    Compute the projected gradient of E with respect to n.

    The functional derivative has two parts:

    (a) Sigma-model term:  g2_a = -Laplacian(n_a)

    (b) Skyrme term:
        g4_a = sum_{i<j} C_ij (d_i n x d_j n)_a
             - sum_i d_i { sum_{j != i} sgn(i,j) C_ij [n x d_j n]_a }

        where sgn(i,j) = +1 if i<j, -1 if i>j.

    The total unprojected gradient is g_a = g2_a + g4_a.
    We then project onto the S^2 tangent plane:
        g_a  ->  g_a - (g . n) n_a

    Returns
    -------
    grad : (3, N, N, N) projected gradient array.
    """
    dn = compute_all_derivs(n, h)

    # ---- Sigma-model gradient: g2_a = -Laplacian(n_a) -------------------
    grad = np.empty_like(n)
    for a in range(3):
        grad[a] = -laplacian(n[a], h)

    # ---- Skyrme gradient -------------------------------------------------
    C_xy, C_xz, C_yz = compute_Cij(n, dn)
    # Pack C into a structure: C[(i,j)] with i<j
    # (0,1)=xy, (0,2)=xz, (1,2)=yz
    C_arr = {(0, 1): C_xy, (0, 2): C_xz, (1, 2): C_yz}

    partials = [partial_x, partial_y, partial_z]

    # Part 1: Direct term from n_a appearing in C_ij
    # For each (i,j) pair, (d_i n x d_j n)_a is the a-th component of the
    # cross product of d_i n and d_j n.
    # (d_i n x d_j n)_a = eps_{abc} (d_i n_b)(d_j n_c)
    for (i, j), Cij in C_arr.items():
        # Cross product components: eps_{abc} d_i n_b d_j n_c
        cross_0 = dn[1, i] * dn[2, j] - dn[2, i] * dn[1, j]
        cross_1 = dn[2, i] * dn[0, j] - dn[0, i] * dn[2, j]
        cross_2 = dn[0, i] * dn[1, j] - dn[1, i] * dn[0, j]
        grad[0] += Cij * cross_0
        grad[1] += Cij * cross_1
        grad[2] += Cij * cross_2

    # Part 2: Divergence term from d_i n_a appearing in C_ij
    # We need: -sum_i d_i { sum_{j != i} sgn(i,j) C_ij [n x d_j n]_a }
    #
    # For each spatial direction i, collect all j != i:
    #   flux_i_a = sum_{j != i} sgn(i,j) C_ij [n x d_j n]_a
    #
    # Then subtract d_i(flux_i_a) from grad[a].
    #
    # [n x d_j n]_a = eps_{abc} n_b (d_j n_c)

    for i in range(3):
        # Accumulate flux for direction i
        flux = np.zeros((3, n.shape[1], n.shape[2], n.shape[3]),
                        dtype=np.float64)
        for j in range(3):
            if j == i:
                continue
            # Determine sign and get C_ij
            if i < j:
                sgn = 1.0
                Cij = C_arr[(i, j)]
            else:
                sgn = -1.0
                Cij = C_arr[(j, i)]

            # [n x d_j n]_a = eps_{abc} n_b d_j n_c
            cross_0 = n[1] * dn[2, j] - n[2] * dn[1, j]
            cross_1 = n[2] * dn[0, j] - n[0] * dn[2, j]
            cross_2 = n[0] * dn[1, j] - n[1] * dn[0, j]

            flux[0] += sgn * Cij * cross_0
            flux[1] += sgn * Cij * cross_1
            flux[2] += sgn * Cij * cross_2

        # Add divergence d_i(flux_i_a) -- the EL equation gives:
        #   delta E4 / delta n_a = Part1 + sum_i d_i(flux_i_a)
        # where flux_i_a = sum_{j!=i} sgn(i,j) C_ij [n x d_j n]_a.
        # Derivation: d(eps4)/d(d_i n_beta) = -flux_i_beta, so the EL
        # divergence term -d_i(d(eps4)/d(d_i n_beta)) = +d_i(flux_i_beta).
        for a in range(3):
            grad[a] += partials[i](flux[a], h)

    # ---- Project onto tangent plane of S^2 --------------------------------
    # g_a -> g_a - (g . n) n_a
    gdotn = grad[0] * n[0] + grad[1] * n[1] + grad[2] * n[2]
    for a in range(3):
        grad[a] -= gdotn * n[a]

    return grad


# =========================================================================
#  6. HOPF CHARGE COMPUTATION (FFT-BASED WHITEHEAD INTEGRAL)
# =========================================================================

def compute_hopf_charge(n, h):
    """
    Compute the Hopf invariant via the Whitehead formula:

      H = (1 / 4 pi^2) integral A . B d^3x

    where B_i = (1/2) eps_{ijk} F_{jk},  F_{ij} = n . (d_i n x d_j n),
    and A is found by solving curl A = B in Fourier space (Coulomb gauge).

    In Fourier space:  A_hat(k) = i k x B_hat(k) / |k|^2.

    Returns
    -------
    H : float   The Hopf charge (should be close to 1 for the Hopf map).
    """
    dn = compute_all_derivs(n, h)
    C_xy, C_xz, C_yz = compute_Cij(n, dn)

    # B_i = (1/2) eps_{ijk} F_{jk}
    # F_{01}=C_xy, F_{02}=C_xz, F_{12}=C_yz (and antisymmetric)
    # B_0 = (1/2)(eps_{012} F_{12} + eps_{021} F_{21}) = (1/2)(F_{12} - F_{21}) = F_{12} = C_yz
    # B_1 = (1/2)(eps_{102} F_{02} + eps_{120} F_{20}) = (1/2)(-F_{02} + F_{20}) = -C_xz
    # B_2 = (1/2)(eps_{201} F_{01} + eps_{210} F_{10}) = (1/2)(F_{01} - F_{10}) = C_xy
    # Actually more carefully with the standard convention:
    # B_x = (1/2)(eps_{xyz} F_{yz} + eps_{xzy} F_{zy}) = F_{yz} = C_yz
    # B_y = (1/2)(eps_{yxz} F_{xz} + eps_{yzx} F_{zx}) = -F_{xz} = -C_xz
    # Wait -- let me be precise.
    # eps_{ijk}: B_i = (1/2) sum_{j<k or j>k} eps_{ijk} F_{jk}
    # B_0 = eps_{012} F_{12} + eps_{021} F_{21} = 1*C_yz + (-1)*(-C_yz) = 2 C_yz -> /2 = C_yz
    # B_1 = eps_{102} F_{02} + eps_{120} F_{20} = (-1)*C_xz + (1)*(-C_xz) = -2 C_xz -> /2 = -C_xz
    # B_2 = eps_{201} F_{01} + eps_{210} F_{10} = (1)*C_xy + (-1)*(-C_xy) = 2 C_xy -> /2 = C_xy
    Bx = C_yz
    By = -C_xz
    Bz = C_xy

    # Fourier transforms
    Bx_hat = np.fft.fftn(Bx)
    By_hat = np.fft.fftn(By)
    Bz_hat = np.fft.fftn(Bz)

    # Wave vectors
    N = n.shape[1]
    kx = np.fft.fftfreq(N, d=h) * 2.0 * np.pi
    ky = np.fft.fftfreq(N, d=h) * 2.0 * np.pi
    kz = np.fft.fftfreq(N, d=h) * 2.0 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0  # avoid division by zero (zero mode)

    # A_hat = i k x B_hat / |k|^2  (Coulomb gauge: k . A = 0)
    # (k x B)_x = ky Bz - kz By
    # (k x B)_y = kz Bx - kx Bz
    # (k x B)_z = kx By - ky Bx
    Ax_hat = 1j * (KY * Bz_hat - KZ * By_hat) / k2
    Ay_hat = 1j * (KZ * Bx_hat - KX * Bz_hat) / k2
    Az_hat = 1j * (KX * By_hat - KY * Bx_hat) / k2

    # Zero mode: set to zero (gauge choice)
    Ax_hat[0, 0, 0] = 0.0
    Ay_hat[0, 0, 0] = 0.0
    Az_hat[0, 0, 0] = 0.0

    # Inverse transform
    Ax = np.real(np.fft.ifftn(Ax_hat))
    Ay = np.real(np.fft.ifftn(Ay_hat))
    Az = np.real(np.fft.ifftn(Az_hat))

    # Hopf charge: H = (1/4pi^2) integral A . B d^3x
    dV = h**3
    integrand = Ax * Bx + Ay * By + Az * Bz
    H = np.sum(integrand) * dV / (4.0 * np.pi**2)

    return H


def compute_hopf_charge_proxy(n, h):
    """
    Quick proxy for monitoring topology changes during flow.

    Q = (1/4pi) integral n . (d_x n x d_y n) d^3x

    This integrates the z-plane mapping degree and is cheaper than the
    full FFT-based Whitehead integral.

    Returns
    -------
    Q : float   Topology proxy (not the exact Hopf charge).
    """
    dx_n0 = partial_x(n[0], h)
    dx_n1 = partial_x(n[1], h)
    dx_n2 = partial_x(n[2], h)
    dy_n0 = partial_y(n[0], h)
    dy_n1 = partial_y(n[1], h)
    dy_n2 = partial_y(n[2], h)

    # n . (d_x n x d_y n)
    cross_0 = dx_n1 * dy_n2 - dx_n2 * dy_n1
    cross_1 = dx_n2 * dy_n0 - dx_n0 * dy_n2
    cross_2 = dx_n0 * dy_n1 - dx_n1 * dy_n0
    integrand = n[0] * cross_0 + n[1] * cross_1 + n[2] * cross_2

    dV = h**3
    Q = np.sum(integrand) * dV / (4.0 * np.pi)
    return Q


# =========================================================================
#  7. ARRESTED DAMPED NEWTON FLOW
# =========================================================================

def arrested_flow(n, h, n_steps=6000, dt0=5e-4, gamma=0.5,
                  print_every=100, converge_tol=1e-7, converge_window=50):
    """
    Minimise the FN energy using normalised gradient descent with
    topology-preserving backtracking line search.

    Algorithm:
      For each step:
        1. Compute projected gradient g (tangent to S^2)
        2. Normalise: direction d = g / max(|g|) so max change = dt (radians)
        3. Backtracking: try n - dt*d, halving dt up to 8x until E decreases
           AND virial ratio doesn't blow up AND Q proxy stays stable
        4. Accept/reject with adaptive dt
        5. Convergence: |dE/E| < tol for converge_window steps

    Parameters
    ----------
    n       : (3, N, N, N) initial field configuration.
    h       : float grid spacing.
    n_steps : int maximum iterations.
    dt0     : float initial max angular step (radians).
    gamma   : float (unused, kept for interface compatibility).
    print_every : int progress print frequency.
    converge_tol : float convergence threshold for |dE/E|.
    converge_window : int number of consecutive converged steps required.

    Returns
    -------
    n       : (3, N, N, N) final field (best topology-preserved energy).
    history : dict with 'E', 'E2', 'E4', 'virial', 'step', 'grad_norm'.
    """
    n = n.copy()
    N = n.shape[1]
    dt = dt0
    dt_max = dt0 * 3.0
    dt_min = 1e-7

    # Compute initial energy
    E_prev, E2, E4, _, _, _ = compute_energy(n, h)
    virial = E2 / (E4 + 1e-30)

    # Initial topology proxy
    Q_init = compute_hopf_charge_proxy(n, h)
    Q_current = Q_init
    virial_init = virial
    virial_max = virial * 1.10  # allow 10% virial increase from initial

    # Track best configuration
    n_best = n.copy()
    E_best = E_prev

    history = {
        'E': [E_prev], 'E2': [E2], 'E4': [E4],
        'virial': [virial], 'step': [0], 'grad_norm': [0.0]
    }

    consec_converge = 0
    consec_accept = 0
    t_start = time.time()

    print("  Normalised gradient descent (topology-preserving)")
    print("  N = %d, h = %.4f, dt0 = %.2e, virial_cap = %.2f (strict)"
          % (N, h, dt0, virial_init))
    print("  Initial energy: E = %.4f (E2 = %.2f, E4 = %.2f)" % (E_prev, E2, E4))
    print("  Initial Q_proxy = %.4f" % Q_init)
    print()
    header = ("  %6s  %12s  %10s  %10s  %7s  %10s  %10s"
              % ("Step", "Energy", "E2", "E4", "E2/E4", "|grad|", "dt"))
    print(header)
    print("  " + "-" * (len(header) - 2))

    rejects = 0
    Q_check_interval = 50  # check full Q proxy every N steps

    for step in range(1, n_steps + 1):
        # 1. Compute projected gradient
        grad = compute_gradient(n, h)
        gnorm = np.sqrt(np.sum(grad**2) * h**3)

        # 2. Normalise gradient by max pointwise magnitude
        grad_mag = np.sqrt(grad[0]**2 + grad[1]**2 + grad[2]**2)
        gmax = grad_mag.max()
        if gmax < 1e-15:
            elapsed = time.time() - t_start
            print()
            print("  Gradient vanished at step %d. |grad| = %.2e" % (step, gnorm))
            print("  Final E = %.4f, Time: %.1fs" % (E_prev, elapsed))
            break

        d = grad / gmax  # max |d| = 1 at any point

        # 3. Backtracking line search with topology guard
        accepted = False
        dt_try = dt
        for ls in range(10):
            n_new = n - dt_try * d
            # Renormalise to S^2
            norm = np.sqrt(n_new[0]**2 + n_new[1]**2 + n_new[2]**2)
            norm = np.maximum(norm, 1e-30)
            for a in range(3):
                n_new[a] /= norm

            E_new, E2_new, E4_new, _, _, _ = compute_energy(n_new, h)

            if E_new < E_prev:
                # Cheap topology guard: virial ratio cap
                v_new = E2_new / (E4_new + 1e-30)
                if v_new > virial_max:
                    dt_try *= 0.5
                    continue

                # Periodic expensive topology check
                if step % Q_check_interval == 0:
                    Q_new = compute_hopf_charge_proxy(n_new, h)
                    if abs(Q_new - Q_init) > 0.08 * abs(Q_init):
                        dt_try *= 0.5
                        continue
                    Q_current = Q_new

                accepted = True
                break
            dt_try *= 0.5

        if accepted:
            n = n_new
            dE_rel = (E_new - E_prev) / (abs(E_prev) + 1e-30)
            E_prev = E_new
            E2, E4 = E2_new, E4_new
            virial = E2 / (E4 + 1e-30)
            dt = dt_try

            # Ratchet virial cap downward only -- once virial decreases,
            # don't let it go back up
            if virial < virial_max:
                virial_max = virial

            # Track best
            if E_new < E_best:
                n_best = n.copy()
                E_best = E_new

            # Grow dt after consecutive accepts
            consec_accept += 1
            if consec_accept >= 10:
                dt = min(dt * 1.03, dt_max)
                consec_accept = 0

            # Record history
            if step % 10 == 0 or step <= 20:
                history['E'].append(E_prev)
                history['E2'].append(E2)
                history['E4'].append(E4)
                history['virial'].append(virial)
                history['step'].append(step)
                history['grad_norm'].append(gnorm)

            # Progress output
            if step % print_every == 0:
                elapsed = time.time() - t_start
                print("  %6d  %12.4f  %10.4f  %10.4f  %7.4f  %10.2e  %10.2e  [%.1fs]"
                      % (step, E_prev, E2, E4, virial, gnorm, dt, elapsed))

            # Convergence check
            if abs(dE_rel) < converge_tol:
                consec_converge += 1
                if consec_converge >= converge_window:
                    elapsed = time.time() - t_start
                    print()
                    print("  Converged at step %d: |dE/E| = %.2e for %d steps"
                          % (step, abs(dE_rel), converge_window))
                    print("  Final E = %.4f, E2/E4 = %.4f, |grad| = %.2e"
                          % (E_prev, virial, gnorm))
                    print("  Time: %.1fs" % elapsed)
                    break
            else:
                consec_converge = 0

        else:
            rejects += 1
            consec_accept = 0
            dt = max(dt_try, dt_min)

            if dt <= dt_min:
                dt = dt0 * 0.1  # reset

            if rejects > n_steps // 3:
                elapsed = time.time() - t_start
                print()
                print("  Too many rejects (%d). E = %.4f, Time: %.1fs"
                      % (rejects, E_prev, elapsed))
                break

    else:
        elapsed = time.time() - t_start
        print()
        print("  Reached max steps (%d). E = %.4f, E2/E4 = %.4f"
              % (n_steps, E_prev, virial))
        print("  Time: %.1fs" % elapsed)

    # Final topology check
    H_final = compute_hopf_charge(n_best, h)
    Q_final = compute_hopf_charge_proxy(n_best, h)
    print("  Final topology: H = %.4f, Q = %.4f (init: %.4f)" %
          (H_final, Q_final, Q_init))

    # Return best configuration
    if E_best < E_prev:
        n = n_best

    # Convert to arrays
    for key in history:
        history[key] = np.array(history[key])

    return n, history


# =========================================================================
#  7a. DERRICK RESCALING + TOPOLOGY-SAFE GRADIENT FLOW
# =========================================================================

def derrick_rescale(n, h, alpha=0.1):
    """
    Partial Derrick rescaling:  n(x) -> n(lam*x) with
      lam_opt = sqrt(E2/E4)
      lam_applied = 1 + alpha * (lam_opt - 1)

    Uses cubic interpolation (scipy) or nearest-neighbor fallback.
    Returns (n_new, lam_applied, E_new, E2_new, E4_new).
    """
    E, E2, E4, _, _, _ = compute_energy(n, h)
    lam_opt = np.sqrt(E2 / (E4 + 1e-30))
    lam = 1.0 + alpha * (lam_opt - 1.0)

    if abs(lam - 1.0) < 0.005:
        return n, 1.0, E, E2, E4

    N = n.shape[1]
    c = (N - 1) / 2.0
    ix = np.arange(N, dtype=np.float64)
    IX, IY, IZ = np.meshgrid(ix, ix, ix, indexing='ij')

    IX_s = c + (IX - c) * lam
    IY_s = c + (IY - c) * lam
    IZ_s = c + (IZ - c) * lam

    coords = np.array([IX_s.ravel(), IY_s.ravel(), IZ_s.ravel()])

    n_new = np.empty_like(n)
    if HAS_SCIPY:
        for a in range(3):
            n_new[a] = _map_coordinates(
                n[a], coords, order=3, mode='nearest'
            ).reshape(N, N, N)
    else:
        # Nearest-neighbor fallback
        ci = np.clip(np.round(coords).astype(int), 0, N - 1)
        for a in range(3):
            n_new[a] = n[a][ci[0], ci[1], ci[2]].reshape(N, N, N)

    # Renormalise to S^2
    nrm = np.sqrt(np.sum(n_new**2, axis=0, keepdims=True))
    n_new /= np.maximum(nrm, 1e-30)

    E_n, E2_n, E4_n, _, _, _ = compute_energy(n_new, h)
    return n_new, lam, E_n, E2_n, E4_n


def compute_derrick_direction(n, h, x1d):
    """
    Infinitesimal Derrick flow direction: d_D = (x . nabla) n.

    Under n -> n + eps * d_D:
      E2 -> E2(1 - eps) + O(eps^2)
      E4 -> E4(1 + eps) + O(eps^2)
      dE = eps * (E4 - E2)  (decreases E when E2 > E4)

    This avoids interpolation entirely -- uses only existing derivatives.
    Returns d_D projected onto T(S^2) and normalised per-point.
    """
    N = n.shape[1]
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing='ij')
    dn = compute_all_derivs(n, h)

    d = np.zeros_like(n)
    for a in range(3):
        d[a] = X * dn[a, 0] + Y * dn[a, 1] + Z * dn[a, 2]

    # Project onto tangent plane of S^2
    gdotn = d[0] * n[0] + d[1] * n[1] + d[2] * n[2]
    for a in range(3):
        d[a] -= gdotn * n[a]

    return d


def topology_safe_flow(n, h, x1d, n_steps=5000, dt0=1e-4,
                       derrick_weight=0.0,
                       topo_check_every=5, topo_tol=0.03,
                       print_every=100):
    """
    Gradient flow with:
      - Infinitesimal Derrick flow mixed into the gradient for virial control
      - Frequent topology monitoring via Q proxy and periodic full H check
      - Adaptive dt based on energy decrease + topology preservation

    The combined direction is:
      d = grad_EL - mu * d_Derrick
    where mu is proportional to (virial - 1) * derrick_weight.
    When virial > 1 (E2 > E4), the Derrick component pushes toward
    compression (increasing E4, decreasing E2).

    Parameters
    ----------
    n              : (3, N, N, N) initial field.
    h              : float grid spacing.
    x1d            : (N,) coordinate array.
    n_steps        : int maximum iterations.
    dt0            : float initial angular step size.
    derrick_weight : float strength of Derrick correction (0 = off).
    topo_check_every : int steps between Q proxy checks.
    topo_tol       : float max allowed |Q-Q0|/|Q0| per topo check.
    print_every    : int progress output frequency.

    Returns
    -------
    n       : (3, N, N, N) best field (topology preserved).
    history : dict.
    """
    n = n.copy()
    N = n.shape[1]
    dt = dt0
    dt_max = dt0 * 5.0
    dt_min = 1e-8

    E, E2, E4, _, _, _ = compute_energy(n, h)
    Q0 = compute_hopf_charge_proxy(n, h)
    H0 = compute_hopf_charge(n, h)
    virial = E2 / (E4 + 1e-30)

    n_best = n.copy()
    E_best = E

    history = {
        'E': [E], 'E2': [E2], 'E4': [E4],
        'virial': [virial], 'step': [0], 'grad_norm': [0.0],
        'H': [H0], 'dt': [dt]
    }

    t_start = time.time()
    consec_accept = 0

    print("  Topology-safe gradient flow + Derrick direction")
    print("  N=%d, h=%.4f, dt0=%.2e, derrick_w=%.2f" %
          (N, h, dt0, derrick_weight))
    print("  Topo check: every %d steps, tol=%.3f" % (topo_check_every, topo_tol))
    print("  Initial: E=%.2f (E2=%.2f E4=%.2f E2/E4=%.2f)" % (E, E2, E4, virial))
    print("  Initial: H=%.4f Q=%.4f" % (H0, Q0))
    print()
    print("  %6s  %10s  %8s  %8s  %7s  %10s  %10s  %6s" %
          ("Step", "E", "E2", "E4", "E2/E4", "|grad|", "dt", "Q"))
    print("  " + "-" * 75)

    for step in range(1, n_steps + 1):
        # 1. Compute EL gradient
        grad = compute_gradient(n, h)
        gnorm = np.sqrt(np.sum(grad**2) * h**3)

        # 2. Optionally mix in Derrick direction
        if derrick_weight > 0 and virial > 1.2:
            d_D = compute_derrick_direction(n, h, x1d)
            mu = derrick_weight * (virial - 1.0)
            # Subtract Derrick from gradient: we step n -= dt*d,
            # and Derrick direction should be ADDED to n (n += eps*d_D),
            # so we subtract it from d.
            d_combined = grad - mu * d_D
        else:
            d_combined = grad

        # 3. Normalise by max pointwise magnitude
        d_mag = np.sqrt(d_combined[0]**2 + d_combined[1]**2 + d_combined[2]**2)
        dmax = d_mag.max()
        if dmax < 1e-15:
            print("  Direction vanished at step %d" % step)
            break
        d = d_combined / dmax

        # 4. Line search with energy + topology checks
        accepted = False
        dt_try = dt
        for ls in range(12):
            n_try = n - dt_try * d
            nrm = np.sqrt(n_try[0]**2 + n_try[1]**2 + n_try[2]**2)
            nrm = np.maximum(nrm, 1e-30)
            for a in range(3):
                n_try[a] /= nrm

            E_try, E2_try, E4_try, _, _, _ = compute_energy(n_try, h)
            if E_try >= E:
                dt_try *= 0.5
                continue

            # Topology check
            if step % topo_check_every == 0:
                Q_try = compute_hopf_charge_proxy(n_try, h)
                if abs(Q_try - Q0) > topo_tol * abs(Q0):
                    dt_try *= 0.5
                    continue

            accepted = True
            break

        if accepted:
            n = n_try
            E, E2, E4 = E_try, E2_try, E4_try
            virial = E2 / (E4 + 1e-30)
            dt = dt_try
            consec_accept += 1
            if consec_accept >= 8:
                dt = min(dt * 1.05, dt_max)
                consec_accept = 0
            if E < E_best:
                E_best = E
                n_best = n.copy()
        else:
            consec_accept = 0
            dt = max(dt_try, dt_min)

        # Full Hopf check periodically
        if step % (topo_check_every * 5) == 0:
            H_now = compute_hopf_charge(n, h)
            history['H'].append(H_now)
            if abs(H_now - H0) > 0.15 * abs(H0):
                print("  ** Topology lost at step %d: H=%.4f (was %.4f), stopping" %
                      (step, H_now, H0))
                n = n_best
                break

        # Logging
        if step % 10 == 0:
            history['E'].append(E)
            history['E2'].append(E2)
            history['E4'].append(E4)
            history['virial'].append(virial)
            history['step'].append(step)
            history['grad_norm'].append(gnorm)
            history['dt'].append(dt)

        if step % print_every == 0:
            elapsed = time.time() - t_start
            Q_now = compute_hopf_charge_proxy(n, h)
            print("  %6d  %10.4f  %8.2f  %8.2f  %7.3f  %10.2e  %10.2e  %6.3f  [%.0fs]"
                  % (step, E, E2, E4, virial, gnorm, dt, Q_now, elapsed))

    elapsed = time.time() - t_start
    H_final = compute_hopf_charge(n_best, h)
    Q_final = compute_hopf_charge_proxy(n_best, h)
    print()
    print("  Done in %.1fs (%.1f min)" % (elapsed, elapsed / 60))
    print("  Best E=%.4f, H=%.4f, Q=%.4f (init H=%.4f Q=%.4f)" %
          (E_best, H_final, Q_final, H0, Q0))

    for key in history:
        history[key] = np.array(history[key])
    return n_best, history


# =========================================================================
#  7b. PERTURBATION AND CHECKPOINT UTILITIES
# =========================================================================

def perturb_field(n, amplitude=0.15, fraction=0.05, rng=None):
    """
    Apply random S^2-preserving perturbations to a subset of grid points.

    At each selected point, apply a small SO(3) rotation to the n-vector.
    This pushes the configuration into a different energy basin while
    preserving the unit-vector constraint |n| = 1 everywhere.

    Parameters
    ----------
    n         : (3, N, N, N) field array.
    amplitude : float   Max rotation angle in radians.
    fraction  : float   Fraction of grid points to perturb (0 to 1).
    rng       : numpy Generator or None.

    Returns
    -------
    n_new : (3, N, N, N) perturbed field with |n| = 1 everywhere.
    """
    if rng is None:
        rng = np.random.default_rng()

    n_new = n.copy()
    N = n.shape[1]
    total_pts = N ** 3

    # Select random subset of points to perturb
    n_perturb = max(1, int(fraction * total_pts))
    flat_idx = rng.choice(total_pts, size=n_perturb, replace=False)

    # Convert flat indices to 3D indices
    ix = flat_idx // (N * N)
    iy = (flat_idx % (N * N)) // N
    iz = flat_idx % N

    # Generate random rotation axes (unit vectors on S^2)
    phi = rng.uniform(0, 2 * np.pi, n_perturb)
    cos_theta = rng.uniform(-1, 1, n_perturb)
    sin_theta = np.sqrt(1.0 - cos_theta**2)
    ax = sin_theta * np.cos(phi)
    ay = sin_theta * np.sin(phi)
    az = cos_theta

    # Generate random rotation angles
    angles = rng.uniform(-amplitude, amplitude, n_perturb)

    # Apply Rodrigues' rotation formula at each selected point
    # v_rot = v cos(theta) + (k x v) sin(theta) + k (k.v)(1 - cos(theta))
    for idx in range(n_perturb):
        i, j, k = ix[idx], iy[idx], iz[idx]
        v = np.array([n_new[0, i, j, k], n_new[1, i, j, k], n_new[2, i, j, k]])
        kk = np.array([ax[idx], ay[idx], az[idx]])
        ct = np.cos(angles[idx])
        st = np.sin(angles[idx])
        k_cross_v = np.cross(kk, v)
        k_dot_v = np.dot(kk, v)
        v_rot = v * ct + k_cross_v * st + kk * k_dot_v * (1.0 - ct)
        # Renormalise
        v_rot /= np.linalg.norm(v_rot) + 1e-30
        n_new[0, i, j, k] = v_rot[0]
        n_new[1, i, j, k] = v_rot[1]
        n_new[2, i, j, k] = v_rot[2]

    return n_new


def save_checkpoint(n, h, x1d, energy, hopf_charge, filepath):
    """Save field configuration and metadata to .npz checkpoint."""
    np.savez_compressed(filepath,
                        n=n, h=np.array([h]), x1d=x1d,
                        energy=np.array([energy]),
                        hopf_charge=np.array([hopf_charge]))
    print("  Checkpoint saved: %s (E=%.4f, H=%.4f)" %
          (filepath, energy, hopf_charge))


def load_checkpoint(filepath):
    """
    Load field configuration from .npz checkpoint.

    Returns
    -------
    n : (3, N, N, N) field array.
    h : float grid spacing.
    x1d : (N,) coordinate array.
    energy : float last known energy.
    hopf_charge : float last known Hopf charge.
    """
    data = np.load(filepath)
    n = data['n']
    h = float(data['h'][0])
    x1d = data['x1d']
    energy = float(data['energy'][0])
    hopf_charge = float(data['hopf_charge'][0])
    print("  Checkpoint loaded: %s (E=%.4f, H=%.4f, N=%d)" %
          (filepath, energy, hopf_charge, n.shape[1]))
    return n, h, x1d, energy, hopf_charge


# =========================================================================
#  7c. BASIN-HOPPING MONTE CARLO
# =========================================================================

def basin_hopping(x1d, h, n_basins=15, a_range=(0.3, 3.0),
                  short_steps=400, perturb_rounds=5,
                  perturb_amplitude=0.15, perturb_fraction=0.05,
                  dt0=5e-4, checkpoint_dir=None, seed=42):
    """
    Monte Carlo basin-hopping for the FN soliton energy minimisation.

    Strategy:
      1. Scan Hopf map scale parameter 'a' over a_range — different
         values of 'a' start in different energy basins.
      2. For each 'a', run short arrested_flow, then apply random
         S^2-preserving perturbations and repeat.
      3. Track the best (lowest energy, topology-preserved) configuration
         across all basins.
      4. Final long relaxation from the best starting point found.

    Parameters
    ----------
    x1d         : (N,) grid coordinates.
    h           : float grid spacing.
    n_basins    : int number of Hopf map scale values to try.
    a_range     : (float, float) range of scale parameter a.
    short_steps : int gradient steps per short relaxation.
    perturb_rounds : int number of perturb-relax cycles per basin.
    perturb_amplitude : float max rotation angle for perturbations.
    perturb_fraction  : float fraction of points to perturb.
    dt0         : float initial time step for gradient descent.
    checkpoint_dir : str or None directory for checkpoint files.
    seed        : int random seed for reproducibility.

    Returns
    -------
    n_best    : (3, N, N, N) best field configuration found.
    E_best    : float best energy.
    history   : dict with basin-hopping log.
    """
    rng = np.random.default_rng(seed)
    N = len(x1d)

    if checkpoint_dir is None:
        checkpoint_dir = OUTDIR
    os.makedirs(checkpoint_dir, exist_ok=True)

    a_values = np.linspace(a_range[0], a_range[1], n_basins)
    best_n = None
    best_E = np.inf
    best_H = 0.0
    best_a = 0.0

    basin_log = {
        'a': [], 'E_initial': [], 'E_final': [],
        'H_final': [], 'accepted': []
    }

    print("=" * 70)
    print("  BASIN-HOPPING MONTE CARLO")
    print("  N = %d, basins = %d, a in [%.2f, %.2f]"
          % (N, n_basins, a_range[0], a_range[1]))
    print("  Short steps = %d, perturb rounds = %d"
          % (short_steps, perturb_rounds))
    print("  Perturb: amplitude = %.3f rad, fraction = %.3f"
          % (perturb_amplitude, perturb_fraction))
    print("=" * 70)
    print()

    t_total_start = time.time()

    for ib, a in enumerate(a_values):
        print("-" * 60)
        print("  Basin %d/%d: a = %.3f" % (ib + 1, n_basins, a))
        print("-" * 60)

        # Initial condition from Hopf map with this scale parameter
        n = hopf_map_initial(x1d, a=a)
        E0, E2, E4, _, _, _ = compute_energy(n, h)
        H0 = compute_hopf_charge_proxy(n, h)
        print("  Initial: E = %.2f (E2=%.1f, E4=%.1f), Q = %.4f"
              % (E0, E2, E4, H0))

        basin_log['a'].append(a)
        basin_log['E_initial'].append(E0)

        # Short gradient descent
        n, hist = arrested_flow(
            n, h, n_steps=short_steps, dt0=dt0,
            print_every=max(1, short_steps // 5),
            converge_tol=1e-7, converge_window=30
        )
        E_after = hist['E'][-1]
        print("  After short flow: E = %.4f" % E_after)

        # Perturb-relax cycles
        for pr in range(perturb_rounds):
            n_pert = perturb_field(n, amplitude=perturb_amplitude,
                                   fraction=perturb_fraction, rng=rng)

            # Check topology wasn't destroyed by perturbation
            Q_pert = compute_hopf_charge_proxy(n_pert, h)
            if abs(Q_pert - H0) > 0.15 * abs(H0):
                print("    Perturb %d: Q drifted (%.4f -> %.4f), skip"
                      % (pr + 1, H0, Q_pert))
                continue

            # Short relaxation after perturbation
            n_pert, hist_pert = arrested_flow(
                n_pert, h, n_steps=short_steps // 2, dt0=dt0,
                print_every=max(1, short_steps // 4),
                converge_tol=1e-7, converge_window=20
            )
            E_pert = hist_pert['E'][-1]
            print("    Perturb %d: E = %.4f (delta = %+.4f)"
                  % (pr + 1, E_pert, E_pert - E_after))

            # Accept if energy decreased
            if E_pert < E_after:
                n = n_pert
                E_after = E_pert
                print("    -> Accepted (lower energy)")
            else:
                # Metropolis-like: accept with small probability
                # for basin escape (kT ~ 5% of current energy)
                kT = 0.05 * abs(E_after)
                dE = E_pert - E_after
                prob = np.exp(-dE / (kT + 1e-30))
                if rng.random() < prob:
                    n = n_pert
                    E_after = E_pert
                    print("    -> Accepted (Metropolis, p=%.4f)" % prob)
                else:
                    print("    -> Rejected")

        # Final topology check for this basin
        H_basin = compute_hopf_charge(n, h)
        basin_log['E_final'].append(E_after)
        basin_log['H_final'].append(H_basin)

        topology_ok = abs(H_basin - 1.0) < 0.2
        basin_log['accepted'].append(topology_ok and E_after < best_E)

        print("  Basin %d result: E = %.4f, H = %.4f, topo_ok = %s"
              % (ib + 1, E_after, H_basin, topology_ok))

        if topology_ok and E_after < best_E:
            best_n = n.copy()
            best_E = E_after
            best_H = H_basin
            best_a = a
            print("  *** NEW BEST: E = %.4f (a = %.3f) ***" % (best_E, a))

            # Save checkpoint
            ckpt_path = os.path.join(checkpoint_dir,
                                      'basin_best_N%d.npz' % N)
            save_checkpoint(best_n, h, x1d, best_E, best_H, ckpt_path)

        print()

    # ---- Phase 2: Long relaxation from best basin ----
    if best_n is None:
        print("  WARNING: No basin preserved topology! Using last.")
        best_n = n
        best_E = E_after

    print("=" * 70)
    print("  PHASE 2: Long relaxation from best basin")
    print("  Starting E = %.4f (from a = %.3f, H = %.4f)"
          % (best_E, best_a, best_H))
    print("=" * 70)
    print()

    # Long relaxation: 4x the short steps
    long_steps = short_steps * 4
    best_n, hist_final = arrested_flow(
        best_n, h, n_steps=long_steps, dt0=dt0,
        print_every=max(1, long_steps // 30),
        converge_tol=1e-8, converge_window=50
    )
    best_E = hist_final['E'][-1]
    best_H = compute_hopf_charge(best_n, h)

    # Save final checkpoint
    ckpt_path = os.path.join(checkpoint_dir, 'basin_final_N%d.npz' % N)
    save_checkpoint(best_n, h, x1d, best_E, best_H, ckpt_path)

    elapsed_total = time.time() - t_total_start
    print()
    print("=" * 70)
    print("  BASIN-HOPPING COMPLETE")
    print("  Best energy: E = %.4f (target: 192.5)" % best_E)
    print("  Best Hopf charge: H = %.4f" % best_H)
    print("  Best scale param: a = %.3f" % best_a)
    print("  Total time: %.1f s (%.1f min)" %
          (elapsed_total, elapsed_total / 60.0))
    print("=" * 70)

    # Convert log to arrays
    for key in basin_log:
        basin_log[key] = np.array(basin_log[key])
    basin_log['history_final'] = hist_final

    return best_n, best_E, basin_log


def plot_basin_hopping(basin_log, outdir):
    """Plot basin-hopping results: energy vs scale parameter."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor(BG_COLOR)

    # Panel 1: Energy vs scale parameter a
    ax = axes[0]
    ax.set_facecolor('white')
    a_arr = basin_log['a']
    E_init = basin_log['E_initial']
    E_final = basin_log['E_final']
    accepted = basin_log['accepted']

    ax.scatter(a_arr, E_init, color=GOLD, s=40, alpha=0.5,
               label='Initial (Hopf map)', zorder=3)
    ax.scatter(a_arr[accepted], E_final[accepted], color=TEAL, s=60,
               marker='o', label='Final (topo OK)', zorder=4)
    ax.scatter(a_arr[~accepted], E_final[~accepted], color=CORAL, s=40,
               marker='x', label='Final (topo lost)', zorder=4)

    BS_E = 192.5
    ax.axhline(BS_E, color=PURPLE, linestyle=':', linewidth=2,
               label='BS minimum: %.1f' % BS_E)
    ax.set_xlabel('Hopf map scale $a$', fontsize=12)
    ax.set_ylabel('Energy (soliton units)', fontsize=12)
    ax.set_title('Basin-Hopping: Energy vs Scale', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel 2: Final convergence history
    ax = axes[1]
    ax.set_facecolor('white')
    hist = basin_log['history_final']
    steps = hist['step']
    mask = steps > 0
    if np.sum(mask) > 1:
        ax.plot(steps[mask], hist['E'][mask], color=PURPLE, linewidth=2,
                label='$E_{\\rm total}$')
        ax.plot(steps[mask], hist['E2'][mask], color=TEAL, linewidth=1.5,
                linestyle='--', label='$E_2$')
        ax.plot(steps[mask], hist['E4'][mask], color=CORAL, linewidth=1.5,
                linestyle='--', label='$E_4$')
    ax.axhline(BS_E, color=GOLD, linestyle=':', linewidth=2,
               label='BS: %.1f' % BS_E)
    ax.set_xscale('log')
    ax.set_xlabel('Iteration (final relaxation)', fontsize=12)
    ax.set_ylabel('Energy', fontsize=12)
    ax.set_title('Final Relaxation from Best Basin', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_3d_basin_hopping.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


# =========================================================================
#  7d. NONLINEAR CONJUGATE GRADIENT (POLAK-RIBIERE ON S^2)
# =========================================================================

def cg_flow(n, h, x1d, n_steps=5000, dt0=5e-4,
            restart_every=100, print_every=100,
            topo_check_every=50, topo_tol=0.03,
            global_H_tol=0.20,
            checkpoint_every=500, checkpoint_dir=None):
    """
    Minimise FN energy using nonlinear CG (Polak-Ribiere) on S^2.

    Uses Riemannian CG with vector transport via projection to T_n(S^2).
    Converges much faster than steepest descent for smooth energy landscapes.

    The CG update at step k is:
      1. Compute projected gradient g_k in T_{n_k}(S^2)
      2. Transport g_{k-1} and d_{k-1} to T_{n_k}(S^2) via projection
      3. beta_PR = max(0, <g_k, g_k - g_{k-1}^T> / <g_{k-1}, g_{k-1}>)
      4. d_k = -g_k + beta_PR * d_{k-1}^T
      5. Line search: n_{k+1} = n_k + alpha * d_k / |d_k|_max
      6. Renormalise to S^2

    Topology monitoring (dual-tolerance trust region):
      - Per-chunk: H compared against CHUNK START (topo_tol, default 3%)
        If exceeded, revert to chunk start and halve dt.
      - Global: H compared against INITIAL H0 (global_H_tol, default 20%)
        If exceeded, revert to best-ever and stop.
      - This uses H (reliable) instead of Q proxy (unreliable for FN solitons).

    Parameters
    ----------
    n              : (3, N, N, N) initial field on S^2.
    h              : float grid spacing.
    x1d            : (N,) coordinate array (for checkpoints).
    n_steps        : int maximum iterations.
    dt0            : float initial step size (max angular change per step).
    restart_every  : int CG restart interval (reset to steepest descent).
    print_every    : int progress output frequency.
    topo_check_every : int steps between topology (H) checks.
    topo_tol       : float max per-chunk |dH/H_chunk| drift. If exceeded, revert.
    global_H_tol   : float max global |dH/H0| drift. If exceeded, stop.
    checkpoint_every : int steps between checkpoint saves.
    checkpoint_dir : str or None directory for checkpoints.

    Returns
    -------
    n_best  : (3, N, N, N) best field (lowest E with topology preserved).
    history : dict with convergence data.
    """
    n = n.copy()
    Ng = n.shape[1]
    dV = h**3

    # Initial state
    E, E2, E4, _, _, _ = compute_energy(n, h)
    Q0 = compute_hopf_charge_proxy(n, h)
    H0 = compute_hopf_charge(n, h)
    virial = E2 / (E4 + 1e-30)

    # Initial gradient (projected to S^2 tangent plane)
    grad = compute_gradient(n, h)
    g_norm_sq = np.sum(grad**2) * dV

    # Initial CG direction = steepest descent (negative gradient)
    d = -grad.copy()

    n_best = n.copy()
    E_best = E

    history = {
        'E': [E], 'E2': [E2], 'E4': [E4],
        'virial': [virial], 'step': [0],
        'grad_norm': [np.sqrt(g_norm_sq)],
        'H': [H0], 'dt': [dt0], 'beta': [0.0]
    }

    t_start = time.time()
    consec_accept = 0
    dt = dt0
    dt_max = dt0 * 10.0
    dt_min = 1e-8
    rejects_total = 0
    cg_resets = 0

    print("  Nonlinear Conjugate Gradient (Polak-Ribiere) on S^2")
    print("  N=%d, h=%.4f, dt0=%.2e, CG restart every %d"
          % (Ng, h, dt0, restart_every))
    print("  Topology: per-chunk tol=%.0f%% every %d steps, "
          "global tol=%.0f%%"
          % (topo_tol * 100, topo_check_every, global_H_tol * 100))
    print("  Initial: E=%.2f (E2=%.2f E4=%.2f E2/E4=%.2f)"
          % (E, E2, E4, virial))
    print("  Initial: H=%.4f Q=%.4f" % (H0, Q0))
    print()
    hdr = ("  %6s  %10s  %8s  %8s  %7s  %10s  %10s  %6s"
           % ("Step", "E", "E2", "E4", "E2/E4", "|grad|", "dt", "beta"))
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    # Chunk-based topology trust region state
    n_chunk = n.copy()
    E_chunk = E
    grad_chunk = grad.copy()
    H_chunk = H0  # Hopf charge at chunk start (for per-chunk comparison)
    topo_reverts = 0
    max_topo_reverts = 30  # stop after this many topology reversions

    for step in range(1, n_steps + 1):
        # Normalise direction for step-size control
        d_mag = np.sqrt(d[0]**2 + d[1]**2 + d[2]**2)
        d_max = d_mag.max()
        if d_max < 1e-15:
            print("  Direction vanished at step %d" % step)
            break
        d_hat = d / d_max  # max |d_hat| = 1 at any point

        # Backtracking line search (pure energy descent, no topology check)
        accepted = False
        dt_try = dt
        for ls in range(12):
            n_try = n + dt_try * d_hat
            # Renormalise to S^2
            nrm = np.sqrt(n_try[0]**2 + n_try[1]**2 + n_try[2]**2)
            nrm = np.maximum(nrm, 1e-30)
            for a in range(3):
                n_try[a] /= nrm

            E_try, E2_try, E4_try, _, _, _ = compute_energy(n_try, h)
            if E_try < E:
                accepted = True
                break
            dt_try *= 0.5

        if not accepted:
            rejects_total += 1
            consec_accept = 0
            dt = max(dt_try, dt_min)
            # Reset CG to steepest descent
            d = -grad.copy()
            cg_resets += 1
            if dt <= dt_min:
                dt = dt0 * 0.1  # unstick
            if rejects_total > n_steps // 3:
                print("\n  Too many rejects (%d), stopping" % rejects_total)
                break
            continue

        # Accept step
        n = n_try
        E, E2, E4 = E_try, E2_try, E4_try
        virial = E2 / (E4 + 1e-30)
        dt = dt_try

        # Grow dt after consecutive accepts (conservative to preserve topology)
        consec_accept += 1
        if consec_accept >= 10:
            dt = min(dt * 1.03, dt_max)
            consec_accept = 0

        # New gradient (projected to S^2 tangent plane)
        grad_new = compute_gradient(n, h)
        g_new_norm_sq = np.sum(grad_new**2) * dV

        # Polak-Ribiere beta with vector transport
        # Transport old gradient to new tangent plane of S^2
        gdotn = (grad[0] * n[0] + grad[1] * n[1] + grad[2] * n[2])
        grad_T = grad.copy()
        for a in range(3):
            grad_T[a] -= gdotn * n[a]
        g_diff = grad_new - grad_T
        beta = np.sum(grad_new * g_diff) * dV / (g_norm_sq + 1e-30)
        beta = max(beta, 0.0)  # PR+ (restart if beta < 0)

        if step % restart_every == 0:
            beta = 0.0

        # Transport old CG direction to new tangent plane
        ddotn = d[0] * n[0] + d[1] * n[1] + d[2] * n[2]
        for a in range(3):
            d[a] -= ddotn * n[a]

        # New CG direction
        d = -grad_new + beta * d

        grad = grad_new
        g_norm_sq = g_new_norm_sq

        # ---- Dual-tolerance topology trust region ----
        # Per-chunk: compare H against chunk start (allows cumulative drift)
        # Global: compare H against initial H0 (hard limit on total drift)
        if step % topo_check_every == 0:
            H_now = compute_hopf_charge(n, h)
            history['H'].append(H_now)
            dH_chunk = abs(H_now - H_chunk) / (abs(H_chunk) + 1e-30)
            dH_global = abs(H_now - H0) / (abs(H0) + 1e-30)

            if dH_global > global_H_tol:
                # Hard stop: total drift exceeded global limit
                print("\n  ** GLOBAL topology limit at step %d: H=%.4f "
                      "(dH_global=%.1f%% > %.1f%%), stopping"
                      % (step, H_now, dH_global * 100, global_H_tol * 100))
                print("     Reverting to best (E=%.4f)" % E_best)
                n = n_best
                break

            if dH_chunk > topo_tol:
                topo_reverts += 1
                print("\n  ** Per-chunk topology at step %d: H=%.4f "
                      "(dH_chunk=%.1f%% > %.1f%%, dH_global=%.1f%%), "
                      "reverting (revert #%d)"
                      % (step, H_now, dH_chunk * 100, topo_tol * 100,
                         dH_global * 100, topo_reverts))
                # Revert to chunk start
                n = n_chunk.copy()
                E, E2_r, E4_r, _, _, _ = compute_energy(n, h)
                E2, E4 = E2_r, E4_r
                virial = E2 / (E4 + 1e-30)
                # Halve dt and reduce dt ceiling to prevent regrowth
                dt = max(dt * 0.5, dt_min)
                dt_max = max(dt_max * 0.7, dt0 * 0.5)
                grad = compute_gradient(n, h)
                g_norm_sq = np.sum(grad**2) * dV
                d = -grad.copy()
                cg_resets += 1
                consec_accept = 0
                if topo_reverts >= max_topo_reverts:
                    print("  Too many topology reverts (%d), stopping"
                          % topo_reverts)
                    break
                print("  Continuing with dt=%.2e, dt_max=%.2e"
                      % (dt, dt_max))
                print(hdr)
                print("  " + "-" * (len(hdr) - 2))
            else:
                # Per-chunk OK -- update chunk checkpoint
                n_chunk = n.copy()
                E_chunk = E
                grad_chunk = grad.copy()
                H_chunk = H_now  # NEW: update chunk H for next comparison
                if E < E_best:
                    E_best = E
                    n_best = n.copy()
                print("\n  -- Topology OK at step %d: H=%.4f "
                      "(dH_chunk=%.1f%%, dH_global=%.1f%%), "
                      "E=%.4f (best=%.4f)"
                      % (step, H_now, dH_chunk * 100, dH_global * 100,
                         E, E_best))
                print(hdr)
                print("  " + "-" * (len(hdr) - 2))

        # Periodic checkpoint save
        if checkpoint_dir and step % checkpoint_every == 0:
            save_checkpoint(n, h, x1d, E, H0,
                            os.path.join(checkpoint_dir,
                                         'cg_N%d_step%d.npz' % (Ng, step)))

        # History logging
        if step % 10 == 0:
            history['E'].append(E)
            history['E2'].append(E2)
            history['E4'].append(E4)
            history['virial'].append(virial)
            history['step'].append(step)
            history['grad_norm'].append(np.sqrt(g_new_norm_sq))
            history['dt'].append(dt)
            history['beta'].append(beta)

        # Progress output
        if step % print_every == 0:
            elapsed = time.time() - t_start
            print("  %6d  %10.4f  %8.2f  %8.2f  %7.3f  %10.2e  %10.2e"
                  "  %6.3f  [%.0fs]"
                  % (step, E, E2, E4, virial,
                     np.sqrt(g_new_norm_sq), dt, beta, elapsed))

    else:
        elapsed = time.time() - t_start
        print("\n  Reached max steps (%d). E = %.4f, E2/E4 = %.4f"
              % (n_steps, E, virial))

    elapsed = time.time() - t_start
    H_final = compute_hopf_charge(n_best, h)
    Q_final = compute_hopf_charge_proxy(n_best, h)
    print()
    print("  CG flow complete in %.1fs (%.1f min)" % (elapsed, elapsed / 60))
    print("  Best E=%.4f, H=%.4f, Q=%.4f (init H=%.4f Q=%.4f)" %
          (E_best, H_final, Q_final, H0, Q0))
    print("  CG resets: %d, total rejects: %d" % (cg_resets, rejects_total))

    for key in history:
        history[key] = np.array(history[key])
    return n_best, history


# =========================================================================
#  8. POST-PROCESSING
# =========================================================================

def azimuthal_average(field_3d, x1d):
    """
    Compute azimuthally averaged field in (rho, z) from 3D Cartesian data.

    For each (rho, z) bin, average field values over the azimuthal angle phi.

    Parameters
    ----------
    field_3d : (N, N, N) array on the Cartesian grid.
    x1d      : (N,) coordinate array.

    Returns
    -------
    eps_rz   : (Nr, Nz) azimuthally averaged field.
    rho_1d   : (Nr,) radial coordinate.
    z_1d     : (N,) same as third coordinate axis.
    """
    N = len(x1d)
    h = x1d[1] - x1d[0]
    L = x1d[-1] + h / 2.0

    # rho bins
    Nr = N // 2
    rho_max = L
    drho = rho_max / Nr
    rho_1d = np.linspace(drho / 2.0, rho_max - drho / 2.0, Nr)

    # Build x, y grids for the first two axes
    X, Y = np.meshgrid(x1d, x1d, indexing='ij')
    rho_xy = np.sqrt(X**2 + Y**2)  # (N, N) radial distance

    # Bin index for each (ix, iy) pair
    bin_idx = np.floor(rho_xy / drho).astype(int)
    bin_idx = np.clip(bin_idx, 0, Nr - 1)

    # Azimuthal average: for each z-slice and rho-bin
    eps_rz = np.zeros((Nr, N), dtype=np.float64)
    counts = np.zeros((Nr, N), dtype=np.float64)

    for iz in range(N):
        slice_2d = field_3d[:, :, iz]  # (N, N)
        for ir in range(Nr):
            mask = (bin_idx == ir)
            if np.any(mask):
                eps_rz[ir, iz] = np.mean(slice_2d[mask])
                counts[ir, iz] = np.sum(mask)

    z_1d_out = x1d  # z axis is the same
    return eps_rz, rho_1d, z_1d_out


def extract_properties(n, h, x1d):
    """
    Extract soliton geometric properties from the converged field.

    Returns
    -------
    props : dict with keys:
        'E_total', 'E2', 'E4', 'virial', 'H',
        'rho_0', 'delta', 'aspect',
        'eps_rz', 'rho_1d', 'z_1d',
        'eps_profile_rho', 'eps_profile_z',
        'fwhm_rho', 'fwhm_z', 'C2'
    """
    E_total, E2, E4, eps2, eps4, dn = compute_energy(n, h)
    eps_total = eps2 + eps4

    H = compute_hopf_charge(n, h)
    Q = compute_hopf_charge_proxy(n, h)

    # Azimuthal average of total energy density
    eps_rz, rho_1d, z_1d = azimuthal_average(eps_total, x1d)

    # Core location: maximum of azimuthally-averaged energy density
    idx_rz = np.unravel_index(np.argmax(eps_rz), eps_rz.shape)
    rho_0 = rho_1d[idx_rz[0]]
    z_0 = z_1d[idx_rz[1]]

    # Radial profile at z = z_0
    eps_profile_rho = eps_rz[:, idx_rz[1]]

    # Axial profile at rho = rho_0
    eps_profile_z = eps_rz[idx_rz[0], :]

    # FWHM in rho
    half_max_rho = np.max(eps_profile_rho) / 2.0
    above_rho = rho_1d[eps_profile_rho >= half_max_rho]
    if len(above_rho) >= 2:
        fwhm_rho = above_rho[-1] - above_rho[0]
    else:
        fwhm_rho = rho_1d[1] - rho_1d[0]

    # FWHM in z
    half_max_z = np.max(eps_profile_z) / 2.0
    above_z = z_1d[eps_profile_z >= half_max_z]
    if len(above_z) >= 2:
        fwhm_z = above_z[-1] - above_z[0]
    else:
        fwhm_z = x1d[1] - x1d[0]

    # Shell thickness (geometric mean of FWHM)
    delta = np.sqrt(max(fwhm_rho, h) * max(fwhm_z, h))

    # Aspect ratio
    aspect = rho_0 / (delta + 1e-30)

    # C_2 computation: -pi^2 / (4 A^2)  from Eq. 13.20
    C2 = -np.pi**2 / (4.0 * aspect**2) if aspect > 0.1 else 0.0

    props = {
        'E_total': E_total,
        'E2': E2,
        'E4': E4,
        'virial': E2 / (E4 + 1e-30),
        'H': H,
        'Q_proxy': Q,
        'rho_0': rho_0,
        'z_0': z_0,
        'delta': delta,
        'fwhm_rho': fwhm_rho,
        'fwhm_z': fwhm_z,
        'aspect': aspect,
        'C2': C2,
        'eps_rz': eps_rz,
        'rho_1d': rho_1d,
        'z_1d': z_1d,
        'eps_profile_rho': eps_profile_rho,
        'eps_profile_z': eps_profile_z,
    }
    return props


# =========================================================================
#  9. PLOTTING
# =========================================================================

def plot_convergence(history, outdir):
    """Energy vs iteration with Battye-Sutcliffe target line."""
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    steps = history['step']
    if len(steps) < 2:
        print("  Skipping convergence plot (no iterations)")
        return
    mask = steps > 0  # skip step 0 for log scale

    ax.plot(steps[mask], history['E'][mask], color=PURPLE, linewidth=2,
            label='$E_{\\rm total}$')
    ax.plot(steps[mask], history['E2'][mask], color=TEAL, linewidth=1.5,
            linestyle='--', label='$E_2$ (sigma-model)')
    ax.plot(steps[mask], history['E4'][mask], color=CORAL, linewidth=1.5,
            linestyle='--', label='$E_4$ (Skyrme)')

    BS_E = 192.5
    ax.axhline(BS_E, color=GOLD, linestyle=':', linewidth=2,
               label='Battye-Sutcliffe: %.1f' % BS_E)

    ax.set_xscale('log')
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('Energy (soliton units)', fontsize=12)
    ax.set_title('3D FN Soliton Energy Convergence', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_3d_convergence.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


def plot_energy_density(eps_rz, rho_1d, z_1d, props, outdir):
    """2D contour plot of azimuthally-averaged energy density."""
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    RHO, Z = np.meshgrid(rho_1d, z_1d, indexing='ij')

    colors = ['white', TEAL, PURPLE]
    cmap = LinearSegmentedColormap.from_list('soliton', colors, N=256)

    eps_norm = eps_rz / (np.max(eps_rz) + 1e-30)
    levels = np.linspace(0.05, 1.0, 20)
    cs = ax.contourf(RHO, Z, eps_norm, levels=levels, cmap=cmap)
    ax.contour(RHO, Z, eps_norm, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
               colors='white', linewidths=0.5, alpha=0.5)

    # Mark core
    ax.plot(props['rho_0'], props['z_0'], 'o', color=CORAL, markersize=8,
            markeredgecolor='white', markeredgewidth=1.5, zorder=5)
    ax.annotate(
        "$\\rho_0 = %.2f$\n$\\delta = %.2f$\nA = %.1f"
        % (props['rho_0'], props['delta'], props['aspect']),
        xy=(props['rho_0'], props['z_0']),
        xytext=(props['rho_0'] + 1.2, props['z_0'] + 1.2),
        fontsize=10, color='white',
        bbox=dict(boxstyle='round,pad=0.3', facecolor=PURPLE, alpha=0.8),
        arrowprops=dict(arrowstyle='->', color='white', lw=1.5))

    plt.colorbar(cs, ax=ax,
                 label='$\\varepsilon / \\varepsilon_{\\max}$')
    ax.set_xlabel('$\\rho$ (soliton units)', fontsize=12)
    ax.set_ylabel('$z$ (soliton units)', fontsize=12)
    ax.set_title('FN $|H|=1$ Soliton -- Azimuthal Average', fontsize=13,
                 fontweight='bold')
    rmax = max(3.0 * props['rho_0'], 4.0)
    zmax = max(2.0 * props['rho_0'], 3.0)
    ax.set_xlim(0, rmax)
    ax.set_ylim(-zmax, zmax)
    ax.set_aspect('equal')

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_3d_energy_density.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


def plot_cross_section(props, outdir):
    """1D radial and axial profiles through the torus core with FWHM."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor(BG_COLOR)

    # -- Radial cross-section --
    ax = axes[0]
    ax.set_facecolor('white')
    rho_1d = props['rho_1d']
    eps_rho = props['eps_profile_rho']
    eps_rho_n = eps_rho / (np.max(eps_rho) + 1e-30)
    ax.plot(rho_1d, eps_rho_n, color=CORAL, linewidth=2.5)
    ax.axvline(props['rho_0'], color=PURPLE, linestyle='--', linewidth=1,
               label='$\\rho_0 = %.2f$' % props['rho_0'])

    half = 0.5
    ax.axhline(half, color=GOLD, linestyle=':', linewidth=1, alpha=0.7)
    rl = props['rho_0'] - props['fwhm_rho'] / 2.0
    rr = props['rho_0'] + props['fwhm_rho'] / 2.0
    ax.annotate('', xy=(max(rl, 0), half), xytext=(rr, half),
                arrowprops=dict(arrowstyle='<->', color=GOLD, lw=2))
    ax.text((max(rl, 0) + rr) / 2.0, half + 0.08,
            'FWHM = %.2f' % props['fwhm_rho'],
            ha='center', fontsize=10, color=GOLD, fontweight='bold')

    ax.set_xlabel('$\\rho$ (soliton units)', fontsize=12)
    ax.set_ylabel('$\\varepsilon / \\varepsilon_{\\max}$', fontsize=12)
    ax.set_title('Radial Cross-Section ($z = z_0$)', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, max(3.0 * props['rho_0'], 4.0))

    # -- Axial cross-section --
    ax = axes[1]
    ax.set_facecolor('white')
    z_1d = props['z_1d']
    eps_z = props['eps_profile_z']
    eps_z_n = eps_z / (np.max(eps_z) + 1e-30)
    ax.plot(z_1d, eps_z_n, color=TEAL, linewidth=2.5)
    ax.axvline(props['z_0'], color=PURPLE, linestyle='--', linewidth=1,
               label='$z_0 = %.2f$' % props['z_0'])

    ax.axhline(half, color=GOLD, linestyle=':', linewidth=1, alpha=0.7)
    zl = props['z_0'] - props['fwhm_z'] / 2.0
    zr = props['z_0'] + props['fwhm_z'] / 2.0
    ax.annotate('', xy=(zl, half), xytext=(zr, half),
                arrowprops=dict(arrowstyle='<->', color=GOLD, lw=2))
    ax.text((zl + zr) / 2.0, half + 0.08,
            'FWHM = %.2f' % props['fwhm_z'],
            ha='center', fontsize=10, color=GOLD, fontweight='bold')

    ax.set_xlabel('$z$ (soliton units)', fontsize=12)
    ax.set_ylabel('$\\varepsilon / \\varepsilon_{\\max}$', fontsize=12)
    ax.set_title('Axial Cross-Section ($\\rho = \\rho_0$)', fontsize=13,
                 fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    zmax = max(2.0 * props['rho_0'], 3.0)
    ax.set_xlim(-zmax, zmax)

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_3d_cross_section.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


def plot_virial(history, outdir):
    """Virial ratio E2/E4 vs iteration (should converge to 1)."""
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(BG_COLOR)
    ax.set_facecolor('white')

    steps = history['step']
    mask = steps > 0

    ax.plot(steps[mask], history['virial'][mask], color=PURPLE, linewidth=2)
    ax.axhline(1.0, color=GOLD, linestyle=':', linewidth=2,
               label='Virial equilibrium ($E_2/E_4 = 1$)')

    ax.set_xscale('log')
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_ylabel('$E_2 / E_4$', fontsize=12)
    ax.set_title('Virial Ratio Convergence', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(outdir, 'fn_soliton_3d_virial.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: %s" % path)


# =========================================================================
#  10. MAIN
# =========================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='FN |H|=1 Hopf Soliton 3D Solver')
    parser.add_argument('N', nargs='?', type=int, default=101,
                        help='Grid size (default: 101)')
    parser.add_argument('n_steps', nargs='?', type=int, default=6000,
                        help='Max flow steps (default: 6000)')
    parser.add_argument('--mode', choices=['single', 'basin', 'resume', 'topo', 'cg'],
                        default='single',
                        help='Run mode: single, basin, resume, topo, or cg (default: single)')
    parser.add_argument('--a', type=float, default=0.5,
                        help='Hopf map scale parameter (default: 0.5)')
    parser.add_argument('--dt', type=float, default=5e-4,
                        help='Initial time step (default: 5e-4)')
    parser.add_argument('--basins', type=int, default=15,
                        help='Number of basins for basin-hopping (default: 15)')
    parser.add_argument('--a-min', type=float, default=0.3,
                        help='Min scale param for basin scan (default: 0.3)')
    parser.add_argument('--a-max', type=float, default=3.0,
                        help='Max scale param for basin scan (default: 3.0)')
    parser.add_argument('--short-steps', type=int, default=400,
                        help='Steps per short relaxation in basin mode (default: 400)')
    parser.add_argument('--perturb-rounds', type=int, default=5,
                        help='Perturb-relax cycles per basin (default: 5)')
    parser.add_argument('--checkpoint', type=str, default=None,
                        help='Path to .npz checkpoint for resume mode')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for basin-hopping (default: 42)')
    return parser.parse_args()


def main():
    """Main computation pipeline."""
    args = parse_args()
    N = args.N
    L = 6.0
    dt0 = args.dt

    print("=" * 70)
    print("  Faddeev-Niemi |H|=1 Hopf Soliton -- 3D Cartesian Solver")
    print("  Mode: %s" % args.mode)
    print("=" * 70)
    print()

    # ---- Grid setup ------------------------------------------------------
    x1d, h = setup_grid(N, L)
    print("  Grid: %d x %d x %d = %d points" % (N, N, N, N**3))
    print("  Domain: [%.1f, %.1f]^3,  h = %.4f" % (-L, L, h))
    mem_mb = 10.0 * N**3 * 8.0 / 1e6
    print("  Estimated memory: ~%.0f MB" % mem_mb)
    print()

    # ================================================================
    #  MODE: BASIN-HOPPING
    # ================================================================
    if args.mode == 'basin':
        n_best, E_best, basin_log = basin_hopping(
            x1d, h,
            n_basins=args.basins,
            a_range=(args.a_min, args.a_max),
            short_steps=args.short_steps,
            perturb_rounds=args.perturb_rounds,
            dt0=dt0,
            checkpoint_dir=OUTDIR,
            seed=args.seed
        )

        # Post-processing on best configuration
        print()
        print("  Extracting soliton properties from best basin...")
        t0 = time.time()
        props = extract_properties(n_best, h, x1d)
        print("  Done in %.1fs." % (time.time() - t0))

        # Summary
        BS_E = 192.5
        print()
        print("=" * 60)
        print("  SOLITON PROPERTIES (Basin-Hopping Best)")
        print("=" * 60)
        print("  Energy:          E     = %.2f" % props['E_total'])
        print("                   E2    = %.2f" % props['E2'])
        print("                   E4    = %.2f" % props['E4'])
        print("  Virial ratio:    E2/E4 = %.4f  (target: 1.0)" % props['virial'])
        print("  Hopf charge:     H     = %.4f  (target: 1.0)" % props['H'])
        print("  BS reference:    E_BS  = %.1f" % BS_E)
        print("  Ratio E/E_BS:          = %.4f" % (props['E_total'] / BS_E))
        print()
        print("  Core radius:     rho_0 = %.3f" % props['rho_0'])
        print("  Shell thickness: delta = %.3f" % props['delta'])
        print("  Aspect ratio:    A     = %.3f" % props['aspect'])
        print("  C_2 (Eq. 13.20): C2    = %.4f" % props['C2'])
        print("  QED C_2:                 -0.3285")
        print("=" * 60)

        # Plots
        print()
        print("  Generating figures...")
        os.makedirs(OUTDIR, exist_ok=True)
        plot_basin_hopping(basin_log, OUTDIR)
        hist_final = basin_log['history_final']
        plot_convergence(hist_final, OUTDIR)
        plot_energy_density(props['eps_rz'], props['rho_1d'], props['z_1d'],
                            props, OUTDIR)
        plot_cross_section(props, OUTDIR)
        plot_virial(hist_final, OUTDIR)

        print()
        print("=" * 70)
        print("  COMPUTATION COMPLETE (basin-hopping)")
        print("=" * 70)
        return n_best, props, basin_log

    # ================================================================
    #  MODE: TOPOLOGY-SAFE GRADIENT + DERRICK
    # ================================================================
    elif args.mode == 'topo':
        a_hopf = args.a
        n_steps = args.n_steps

        print("  Building Hopf map (a = %.2f)..." % a_hopf)
        n = hopf_map_initial(x1d, a=a_hopf)
        E0, E2_0, E4_0, _, _, _ = compute_energy(n, h)
        H0 = compute_hopf_charge(n, h)
        print("  Initial: E=%.2f (E2=%.2f E4=%.2f E2/E4=%.2f) H=%.4f" %
              (E0, E2_0, E4_0, E2_0 / (E4_0 + 1e-30), H0))
        print()

        n, history = topology_safe_flow(
            n, h, x1d, n_steps=n_steps, dt0=dt0,
            derrick_weight=0.0,
            topo_check_every=5, topo_tol=0.03,
            print_every=max(1, n_steps // 30)
        )

        # Save checkpoint
        E_final = history['E'][-1]
        H_final = compute_hopf_charge(n, h)
        save_checkpoint(n, h, x1d, E_final, H_final,
                        os.path.join(OUTDIR, 'topo_N%d.npz' % N))

    # ================================================================
    #  MODE: CONJUGATE GRADIENT
    # ================================================================
    elif args.mode == 'cg':
        a_hopf = args.a
        n_steps = args.n_steps

        print("  Building Hopf map (a = %.2f)..." % a_hopf)
        n = hopf_map_initial(x1d, a=a_hopf)
        E0, E2_0, E4_0, _, _, _ = compute_energy(n, h)
        H0 = compute_hopf_charge(n, h)
        print("  Initial: E=%.2f (E2=%.2f E4=%.2f E2/E4=%.2f) H=%.4f" %
              (E0, E2_0, E4_0, E2_0 / (E4_0 + 1e-30), H0))
        print()

        n, history = cg_flow(
            n, h, x1d, n_steps=n_steps, dt0=dt0,
            restart_every=100,
            print_every=max(1, n_steps // 40),
            topo_check_every=50, topo_tol=0.03,
            global_H_tol=0.20,
            checkpoint_every=500,
            checkpoint_dir=OUTDIR
        )

        # Save final checkpoint
        E_final = history['E'][-1]
        H_final = compute_hopf_charge(n, h)
        save_checkpoint(n, h, x1d, E_final, H_final,
                        os.path.join(OUTDIR, 'cg_N%d.npz' % N))

    # ================================================================
    #  MODE: RESUME FROM CHECKPOINT
    # ================================================================
    elif args.mode == 'resume':
        if args.checkpoint is None:
            # Look for default checkpoint
            ckpt_path = os.path.join(OUTDIR, 'basin_final_N%d.npz' % N)
            if not os.path.exists(ckpt_path):
                ckpt_path = os.path.join(OUTDIR, 'basin_best_N%d.npz' % N)
            if not os.path.exists(ckpt_path):
                print("  ERROR: No checkpoint found. Run basin mode first.")
                sys.exit(1)
        else:
            ckpt_path = args.checkpoint

        n, h_ckpt, x1d_ckpt, E_ckpt, H_ckpt = load_checkpoint(ckpt_path)
        # Use checkpoint grid if sizes match
        if len(x1d_ckpt) == N:
            x1d = x1d_ckpt
            h = h_ckpt
        else:
            print("  WARNING: Checkpoint N=%d != requested N=%d, using new grid"
                  % (len(x1d_ckpt), N))

        print("  Resuming from E = %.4f, H = %.4f" % (E_ckpt, H_ckpt))
        print("  Running %d additional steps..." % args.n_steps)
        print()

        n, history = arrested_flow(
            n, h, n_steps=args.n_steps, dt0=dt0,
            print_every=max(1, args.n_steps // 30),
            converge_tol=1e-8, converge_window=50
        )

        # Save updated checkpoint
        E_final = history['E'][-1]
        H_final = compute_hopf_charge(n, h)
        save_checkpoint(n, h, x1d, E_final, H_final,
                        os.path.join(OUTDIR, 'resume_N%d.npz' % N))

        # Fall through to standard post-processing below
        # (reuse same code as single mode)

    # ================================================================
    #  MODE: SINGLE (default, backwards compatible)
    # ================================================================
    else:
        a_hopf = args.a
        n_steps = args.n_steps

        # ---- Initial condition: Hopf map ---------------------------------
        print("  Building Hopf map initial condition (a = %.2f)..." % a_hopf)
        t0 = time.time()
        n = hopf_map_initial(x1d, a=a_hopf)
        print("  Done in %.1fs." % (time.time() - t0))

        # Verify unit-vector constraint
        norm_check = np.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
        print("  |n| range: [%.6f, %.6f]"
              % (np.min(norm_check), np.max(norm_check)))

        # Initial energy
        E0, E2_0, E4_0, _, _, _ = compute_energy(n, h)
        print("  Initial energy: E = %.2f (E2 = %.2f, E4 = %.2f, ratio = %.3f)"
              % (E0, E2_0, E4_0, E2_0 / (E4_0 + 1e-30)))

        # Initial Hopf charge (FFT-based)
        print("  Computing initial Hopf charge (FFT)...")
        t0 = time.time()
        H0 = compute_hopf_charge(n, h)
        Q0 = compute_hopf_charge_proxy(n, h)
        print("  H (Whitehead) = %.4f,  Q (proxy) = %.4f  [%.1fs]"
              % (H0, Q0, time.time() - t0))
        print()

        # ---- Arrested damped Newton flow ---------------------------------
        print("  Starting normalised gradient descent...")
        print()
        n, history = arrested_flow(
            n, h, n_steps=n_steps, dt0=dt0,
            print_every=max(1, n_steps // 30),
            converge_tol=1e-7, converge_window=50
        )
        print()

    # ---- Post-processing (shared by single and resume modes) -------------
    print("  Extracting soliton properties...")
    t0 = time.time()
    props = extract_properties(n, h, x1d)
    print("  Done in %.1fs." % (time.time() - t0))
    print()

    # ---- Final summary ---------------------------------------------------
    BS_E = 192.5
    print("=" * 60)
    print("  SOLITON PROPERTIES (3D Cartesian)")
    print("=" * 60)
    print("  Energy:          E     = %.2f" % props['E_total'])
    print("                   E2    = %.2f" % props['E2'])
    print("                   E4    = %.2f" % props['E4'])
    print("  Virial ratio:    E2/E4 = %.4f  (target: 1.0)" % props['virial'])
    print("  Hopf charge:     H     = %.4f  (target: 1.0)" % props['H'])
    print("  Topology proxy:  Q     = %.4f" % props['Q_proxy'])
    print("  BS reference:    E_BS  = %.1f" % BS_E)
    print("  Ratio E/E_BS:          = %.4f" % (props['E_total'] / BS_E))
    print()
    print("  Core radius:     rho_0 = %.3f" % props['rho_0'])
    print("  Core z-position: z_0   = %.3f" % props['z_0'])
    print("  Shell thickness: delta = %.3f" % props['delta'])
    print("  FWHM(rho):             = %.3f" % props['fwhm_rho'])
    print("  FWHM(z):               = %.3f" % props['fwhm_z'])
    print("  Aspect ratio:    A     = %.3f" % props['aspect'])
    print()
    print("  C_2 (Eq. 13.20): C2    = %.4f" % props['C2'])
    print("  QED C_2:                 -0.3285")
    print("=" * 60)
    print()

    # ---- Plots -----------------------------------------------------------
    print("  Generating figures...")
    os.makedirs(OUTDIR, exist_ok=True)
    plot_convergence(history, OUTDIR)
    plot_energy_density(props['eps_rz'], props['rho_1d'], props['z_1d'],
                        props, OUTDIR)
    plot_cross_section(props, OUTDIR)
    plot_virial(history, OUTDIR)

    print()
    print("=" * 70)
    print("  COMPUTATION COMPLETE")
    print("=" * 70)

    return n, props, history


if __name__ == '__main__':
    n, props, history = main()
