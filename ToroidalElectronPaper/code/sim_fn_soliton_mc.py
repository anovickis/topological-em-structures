#!/usr/bin/env python3
"""
sim_fn_soliton_mc.py -- Monte Carlo minimisation of the Faddeev-Niemi soliton.

Vectorized stride-3 sublattice decomposition for correct parallel updates.
With stride-3, the minimum distance between sublattice sites is 3 in every
coordinate axis. Since the energy density eps(P) depends on n at P and its
6 face-neighbors (distance 1), and the stencil of a site extends to distance 1,
the total dependency radius is 2. With spacing 3, no two sublattice sites
share any dependency -> per-site dE is exact.

Usage:
    python sim_fn_soliton_mc.py [N] [n_sweeps] [--delta D] [--temp T] [--cool R]
"""

import numpy as np
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os, sys, time, argparse

CORAL  = '#e76f51'
TEAL   = '#2a9d8f'
GOLD   = '#e9c46a'
PURPLE = '#a855f7'
BG_COLOR = '#f8f9fa'

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(os.path.dirname(SCRIPT_DIR), 'images')


# =========================================================================
#  GRID + HOPF MAP
# =========================================================================

def setup_grid(N=61, L=6.0):
    h = 2.0 * L / N
    x1d = np.linspace(-L + h / 2.0, L - h / 2.0, N)
    return x1d, h

def hopf_map_initial(x1d, a=1.5):
    N = len(x1d)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing='ij')
    r2 = X**2 + Y**2 + Z**2
    denom = r2 + a**2
    n = np.empty((3, N, N, N), dtype=np.float64)
    n[0] = 2.0 * (X * Z + Y * a) / denom
    n[1] = 2.0 * (Y * Z - X * a) / denom
    n[2] = (X**2 + Y**2 - Z**2 - a**2) / denom
    norm = np.sqrt(n[0]**2 + n[1]**2 + n[2]**2)
    norm = np.maximum(norm, 1e-30)
    for c in range(3):
        n[c] /= norm
    return n


# =========================================================================
#  DERIVATIVES (vectorised)
# =========================================================================

def partial_x(f, h):
    d = np.empty_like(f)
    d[1:-1] = (f[2:] - f[:-2]) / (2.0 * h)
    d[0] = (f[1] - f[0]) / h
    d[-1] = (f[-1] - f[-2]) / h
    return d

def partial_y(f, h):
    d = np.empty_like(f)
    d[:, 1:-1] = (f[:, 2:] - f[:, :-2]) / (2.0 * h)
    d[:, 0] = (f[:, 1] - f[:, 0]) / h
    d[:, -1] = (f[:, -1] - f[:, -2]) / h
    return d

def partial_z(f, h):
    d = np.empty_like(f)
    d[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2.0 * h)
    d[:, :, 0] = (f[:, :, 1] - f[:, :, 0]) / h
    d[:, :, -1] = (f[:, :, -1] - f[:, :, -2]) / h
    return d


# =========================================================================
#  ENERGY: full grid + eps grid
# =========================================================================

def _compute_dn(n, h):
    """Compute all 9 partial derivatives dn[a,i] = d(n_a)/d(x_i)."""
    partials = [partial_x, partial_y, partial_z]
    dn = np.empty((3, 3, n.shape[1], n.shape[2], n.shape[3]))
    for a in range(3):
        for i in range(3):
            dn[a, i] = partials[i](n[a], h)
    return dn

def _cij(n, dn, i, j):
    """Topological current density C_ij = n . (d_i n x d_j n)."""
    return (n[0] * (dn[1, i] * dn[2, j] - dn[2, i] * dn[1, j])
          + n[1] * (dn[2, i] * dn[0, j] - dn[0, i] * dn[2, j])
          + n[2] * (dn[0, i] * dn[1, j] - dn[1, i] * dn[0, j]))

def compute_eps_grid(n, h):
    """Return total energy density eps = eps2 + eps4 at every grid point."""
    dn = _compute_dn(n, h)
    eps2 = np.zeros_like(n[0])
    for a in range(3):
        for i in range(3):
            eps2 += dn[a, i]**2
    eps2 *= 0.5
    C01 = _cij(n, dn, 0, 1)
    C02 = _cij(n, dn, 0, 2)
    C12 = _cij(n, dn, 1, 2)
    eps4 = 0.5 * (C01**2 + C02**2 + C12**2)
    return eps2, eps4

def compute_energy(n, h):
    """Return (E, E2, E4) integrated over the grid."""
    eps2, eps4 = compute_eps_grid(n, h)
    dV = h**3
    E2 = np.sum(eps2) * dV
    E4 = np.sum(eps4) * dV
    return E2 + E4, E2, E4


# =========================================================================
#  DERRICK RESCALING  n(x) -> n(lambda*x)
# =========================================================================

def derrick_rescale(n, h, x1d):
    """
    Apply optimal Derrick rescaling to bring virial ratio E2/E4 toward 1.
    Under n(x)->n(lam*x):  E2 -> E2/lam,  E4 -> lam*E4
    Optimal: lam = sqrt(E2/E4).
    Uses cubic interpolation via scipy.
    Returns (n_new, lam, E_new, E2_new, E4_new).
    """
    E, E2, E4 = compute_energy(n, h)
    lam = np.sqrt(E2 / (E4 + 1e-30))
    if abs(lam - 1.0) < 0.02:
        return n, 1.0, E, E2, E4

    N = n.shape[1]
    c = (N - 1) / 2.0   # grid center in index space
    ix = np.arange(N, dtype=np.float64)
    IX, IY, IZ = np.meshgrid(ix, ix, ix, indexing='ij')

    # Rescaled index coordinates: x_new = center + (x_old - center) * lam
    IX_s = c + (IX - c) * lam
    IY_s = c + (IY - c) * lam
    IZ_s = c + (IZ - c) * lam

    coords = np.array([IX_s.ravel(), IY_s.ravel(), IZ_s.ravel()])

    n_new = np.empty_like(n)
    for a in range(3):
        n_new[a] = map_coordinates(n[a], coords, order=3,
                                   mode='constant', cval=0.0).reshape(N, N, N)

    # Set boundary to vacuum (n3 = -1 for south pole / vacuum)
    # Actually, preserve whatever the Hopf map gives at the boundary
    # Just renormalize to S^2
    nrm = np.sqrt(np.sum(n_new**2, axis=0, keepdims=True))
    # Where norm is too small (outside original domain), set to vacuum
    vacuum_mask = nrm[0] < 0.1
    n_new[0][vacuum_mask] = 0.0
    n_new[1][vacuum_mask] = 0.0
    n_new[2][vacuum_mask] = -1.0
    nrm = np.sqrt(np.sum(n_new**2, axis=0, keepdims=True))
    n_new /= np.maximum(nrm, 1e-30)

    E_n, E2_n, E4_n = compute_energy(n_new, h)
    return n_new, lam, E_n, E2_n, E4_n


# =========================================================================
#  HOPF CHARGE (FFT-based Whitehead integral)
# =========================================================================

def compute_hopf_charge(n, h):
    dn = _compute_dn(n, h)
    Bx = _cij(n, dn, 1, 2)
    By = -_cij(n, dn, 0, 2)
    Bz = _cij(n, dn, 0, 1)
    N = n.shape[1]
    kx = np.fft.fftfreq(N, d=h) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
    k2 = KX**2 + KY**2 + KZ**2
    k2[0, 0, 0] = 1.0
    Bxh, Byh, Bzh = np.fft.fftn(Bx), np.fft.fftn(By), np.fft.fftn(Bz)
    Axh = 1j * (KY * Bzh - KZ * Byh) / k2
    Ayh = 1j * (KZ * Bxh - KX * Bzh) / k2
    Azh = 1j * (KX * Byh - KY * Bxh) / k2
    Axh[0, 0, 0] = Ayh[0, 0, 0] = Azh[0, 0, 0] = 0.0
    Ax = np.real(np.fft.ifftn(Axh))
    Ay = np.real(np.fft.ifftn(Ayh))
    Az = np.real(np.fft.ifftn(Azh))
    return np.sum(Ax * Bx + Ay * By + Az * Bz) * h**3 / (4 * np.pi**2)


# =========================================================================
#  STENCIL ENERGY at sublattice sites
# =========================================================================

# 7-point stencil offsets (center + 6 face neighbors)
_STENCIL_OFFSETS = [(0,0,0), (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)]

def _stencil_energy(eps, si, sj, sk):
    """Sum eps at 7-point stencil for all sublattice sites.
    Returns array of shape (len(si), len(sj), len(sk))."""
    result = eps[np.ix_(si, sj, sk)].copy()
    for di, dj, dk in _STENCIL_OFFSETS[1:]:
        result += eps[np.ix_(si + di, sj + dj, sk + dk)]
    return result


# =========================================================================
#  VECTORISED MC SWEEP (stride-3 sublattice)
# =========================================================================

def mc_sweep(n, h, delta, rng, T=0.0, margin=2):
    """
    One MC sweep over all 27 stride-3 sublattices.

    For each sublattice:
      1. Compute old stencil energy (from current eps grid)
      2. Apply Rodrigues rotation proposals to all sublattice sites
      3. Compute new stencil energy (from updated eps grid)
      4. Accept if dE < 0 (or Metropolis at temperature T)
      5. Restore rejected sites

    Returns (n_accepted, n_proposed).
    """
    Ng = n.shape[1]
    accepted_total = 0
    proposed_total = 0

    for a in range(3):
        for b in range(3):
            for c in range(3):
                # Sublattice site indices (skip margin from boundaries)
                si = np.arange(a, Ng, 3)
                sj = np.arange(b, Ng, 3)
                sk = np.arange(c, Ng, 3)
                # Clip to interior (margin from each boundary)
                si = si[(si >= margin) & (si < Ng - margin)]
                sj = sj[(sj >= margin) & (sj < Ng - margin)]
                sk = sk[(sk >= margin) & (sk < Ng - margin)]
                if len(si) == 0 or len(sj) == 0 or len(sk) == 0:
                    continue

                M = len(si) * len(sj) * len(sk)
                proposed_total += M

                # --- Old eps and stencil energy ---
                eps2, eps4 = compute_eps_grid(n, h)
                eps_old = eps2 + eps4
                old_se = _stencil_energy(eps_old, si, sj, sk)

                # --- Save old field and generate proposals ---
                ix4 = np.ix_(np.arange(3), si, sj, sk)
                n_old = n[ix4].copy()           # (3, Mi, Mj, Mk)

                shape3d = (len(si), len(sj), len(sk))
                # Random unit rotation axes
                ax_raw = rng.standard_normal((3, *shape3d))
                ax_norm = np.sqrt(np.sum(ax_raw**2, axis=0, keepdims=True))
                axes = ax_raw / np.maximum(ax_norm, 1e-30)
                # Random rotation angles
                angles = delta * rng.standard_normal(shape3d)

                # Rodrigues rotation: n' = n cos a + (k x n) sin a + k(k.n)(1-cos a)
                ca = np.cos(angles)[np.newaxis]     # (1, Mi, Mj, Mk)
                sa = np.sin(angles)[np.newaxis]
                kdn = np.sum(axes * n_old, axis=0, keepdims=True)   # (1,Mi,Mj,Mk)
                kcn = np.cross(axes, n_old, axisa=0, axisb=0, axisc=0)  # (3,Mi,Mj,Mk)

                n_prop = n_old * ca + kcn * sa + axes * kdn * (1.0 - ca)
                # Renormalise to S^2
                pnorm = np.sqrt(np.sum(n_prop**2, axis=0, keepdims=True))
                n_prop /= np.maximum(pnorm, 1e-30)

                # --- Apply proposals ---
                n[ix4] = n_prop

                # --- New eps and stencil energy ---
                eps2_new, eps4_new = compute_eps_grid(n, h)
                eps_new = eps2_new + eps4_new
                new_se = _stencil_energy(eps_new, si, sj, sk)

                # --- Accept / reject ---
                dE = new_se - old_se

                if T > 0:
                    accept = (dE < 0) | (rng.random(shape3d) < np.exp(
                        np.minimum(-dE / T, 0.0)))
                else:
                    accept = dE < 0

                reject = ~accept
                nacc = int(np.sum(accept))
                accepted_total += nacc

                # Restore rejected sites
                if nacc < M:
                    n_cur = n[ix4].copy()
                    n_cur[:, reject] = n_old[:, reject]
                    n[ix4] = n_cur

    return accepted_total, proposed_total


# =========================================================================
#  POST-PROCESSING + PLOTTING
# =========================================================================

def azimuthal_average(field_3d, x1d):
    N = len(x1d)
    h = x1d[1] - x1d[0]
    L = x1d[-1] + h / 2.0
    Nr = N // 2
    drho = L / Nr
    rho_1d = np.linspace(drho / 2, L - drho / 2, Nr)
    X, Y = np.meshgrid(x1d, x1d, indexing='ij')
    rho_xy = np.sqrt(X**2 + Y**2)
    bin_idx = np.clip(np.floor(rho_xy / drho).astype(int), 0, Nr - 1)
    eps_rz = np.zeros((Nr, N))
    for iz in range(N):
        s = field_3d[:, :, iz]
        for ir in range(Nr):
            m = bin_idx == ir
            if np.any(m):
                eps_rz[ir, iz] = np.mean(s[m])
    return eps_rz, rho_1d, x1d


def plot_results(history, n, h, x1d, outdir):
    os.makedirs(outdir, exist_ok=True)

    # --- Convergence ---
    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(BG_COLOR); ax.set_facecolor('white')
    sw = np.array(history['sweep'])
    ax.plot(sw, history['E'], color=PURPLE, lw=2, label='$E$')
    ax.plot(sw, history['E2'], color=TEAL, lw=1.5, ls='--', label='$E_2$')
    ax.plot(sw, history['E4'], color=CORAL, lw=1.5, ls='--', label='$E_4$')
    ax.axhline(192.5, color=GOLD, ls=':', lw=2, label='BS: 192.5')
    ax.set_xlabel('MC Sweep'); ax.set_ylabel('Energy')
    ax.set_title('FN Soliton -- MC Convergence', fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fn_soliton_mc_convergence.png'),
                dpi=150, bbox_inches='tight')
    plt.close()

    # --- Energy density (azimuthal average) ---
    eps2, eps4 = compute_eps_grid(n, h)
    eps_tot = eps2 + eps4
    eps_rz, rho_1d, z_1d = azimuthal_average(eps_tot, x1d)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(BG_COLOR); ax.set_facecolor('white')
    RHO, Z = np.meshgrid(rho_1d, z_1d, indexing='ij')
    cmap = LinearSegmentedColormap.from_list('s', ['white', TEAL, PURPLE], 256)
    eps_n = eps_rz / (np.max(eps_rz) + 1e-30)
    cs = ax.contourf(RHO, Z, eps_n, levels=np.linspace(0.05, 1, 20), cmap=cmap)
    ax.contour(RHO, Z, eps_n, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
               colors='white', linewidths=0.5, alpha=0.5)
    plt.colorbar(cs, ax=ax, label='$\\varepsilon/\\varepsilon_{max}$')
    ax.set_xlabel('$\\rho$'); ax.set_ylabel('$z$')
    ax.set_title('FN Soliton Energy Density (MC)', fontweight='bold')
    ax.set_xlim(0, 5); ax.set_ylim(-4, 4); ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fn_soliton_mc_energy_density.png'),
                dpi=150, bbox_inches='tight')
    plt.close()

    # --- Acceptance + virial ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor(BG_COLOR)
    ax = axes[0]; ax.set_facecolor('white')
    ax.plot(sw, np.array(history['acc']) * 100, color=TEAL, lw=1.5)
    ax.set_xlabel('Sweep'); ax.set_ylabel('Acceptance %')
    ax.set_title('MC Acceptance Rate', fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax = axes[1]; ax.set_facecolor('white')
    ax.plot(sw, history['virial'], color=PURPLE, lw=1.5)
    ax.axhline(1.0, color=GOLD, ls=':', lw=2, label='Virial = 1')
    ax.set_xlabel('Sweep'); ax.set_ylabel('$E_2/E_4$')
    ax.set_title('Virial Ratio', fontweight='bold')
    ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fn_soliton_mc_acceptance.png'),
                dpi=150, bbox_inches='tight')
    plt.close()

    # --- Cross-section through energy density peak ---
    mid = n.shape[1] // 2
    eps_xz = eps_tot[:, mid, :]
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(BG_COLOR); ax.set_facecolor('white')
    XX, ZZ = np.meshgrid(x1d, x1d, indexing='ij')
    eps_xz_n = eps_xz / (np.max(eps_xz) + 1e-30)
    cs = ax.contourf(XX, ZZ, eps_xz_n, levels=np.linspace(0.05, 1, 20), cmap=cmap)
    ax.contour(XX, ZZ, eps_xz_n, levels=[0.1, 0.3, 0.5, 0.7, 0.9],
               colors='white', linewidths=0.5, alpha=0.5)
    plt.colorbar(cs, ax=ax, label='$\\varepsilon/\\varepsilon_{max}$')
    ax.set_xlabel('$x$'); ax.set_ylabel('$z$')
    ax.set_title('FN Soliton -- Cross Section y=0 (MC)', fontweight='bold')
    ax.set_xlim(-5, 5); ax.set_ylim(-5, 5); ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'fn_soliton_mc_cross_section.png'),
                dpi=150, bbox_inches='tight')
    plt.close()

    print("  Saved plots to %s" % outdir)


# =========================================================================
#  CHECKPOINT
# =========================================================================

def save_checkpoint(n, history, fname):
    np.savez_compressed(fname, n=n,
                        sweep=history['sweep'], E=history['E'],
                        E2=history['E2'], E4=history['E4'],
                        virial=history['virial'], acc=history['acc'])
    print("  Checkpoint saved: %s" % fname)


# =========================================================================
#  MAIN
# =========================================================================

def main():
    parser = argparse.ArgumentParser(description='FN Soliton MC Solver')
    parser.add_argument('N', type=int, nargs='?', default=41,
                        help='Grid points per axis (default: 41)')
    parser.add_argument('n_sweeps', type=int, nargs='?', default=200,
                        help='Number of MC sweeps (default: 200)')
    parser.add_argument('--delta', type=float, default=0.15,
                        help='Initial rotation step size in radians (default: 0.15)')
    parser.add_argument('--temp', type=float, default=0.0,
                        help='Initial temperature for SA (0 = greedy, default: 0)')
    parser.add_argument('--cool', type=float, default=0.99,
                        help='SA cooling factor per sweep (default: 0.99)')
    parser.add_argument('--a', type=float, default=1.5,
                        help='Hopf map scale parameter (default: 1.5)')
    parser.add_argument('--L', type=float, default=6.0,
                        help='Half-domain size (default: 6.0)')
    parser.add_argument('--resume', type=str, default=None,
                        help='Resume from checkpoint .npz file')
    parser.add_argument('--derrick', type=int, default=50,
                        help='Apply Derrick rescaling every N sweeps (0=off, default: 50)')
    parser.add_argument('--htol', type=float, default=0.0,
                        help='Topology guard: stop if |H| drops below this (0=off)')
    args = parser.parse_args()

    print("=" * 70)
    print("  Faddeev-Niemi Hopf Soliton -- MC Solver (stride-3 vectorised)")
    print("=" * 70)
    print("  N = %d,  sweeps = %d,  delta = %.4f" % (args.N, args.n_sweeps, args.delta))
    print("  L = %.1f,  a_hopf = %.2f" % (args.L, args.a))
    if args.temp > 0:
        print("  SA: T0 = %.2e, cool = %.4f" % (args.temp, args.cool))
    print()

    x1d, h = setup_grid(args.N, args.L)
    print("  Grid: %d^3 = %d,  h = %.4f" % (args.N, args.N**3, h))

    if args.resume and os.path.exists(args.resume):
        print("  Resuming from %s" % args.resume)
        data = np.load(args.resume)
        n = data['n']
        history = {k: list(data[k]) for k in ['sweep', 'E', 'E2', 'E4', 'virial', 'acc']}
        sweep_offset = int(history['sweep'][-1])
        print("  Resuming from sweep %d" % sweep_offset)
    else:
        n = hopf_map_initial(x1d, a=args.a)
        history = {'sweep': [], 'E': [], 'E2': [], 'E4': [], 'virial': [], 'acc': []}
        sweep_offset = 0

    # --- Derrick rescaling to improve virial ratio ---
    E_pre, E2_pre, E4_pre = compute_energy(n, h)
    v_pre = E2_pre / (E4_pre + 1e-30)
    print("  Pre-Derrick: E=%.2f (E2=%.2f, E4=%.2f, E2/E4=%.2f)" %
          (E_pre, E2_pre, E4_pre, v_pre))
    if v_pre > 1.5 or v_pre < 0.67:
        n, lam, E_d, E2_d, E4_d = derrick_rescale(n, h, x1d)
        v_d = E2_d / (E4_d + 1e-30)
        print("  Derrick lambda=%.3f: E=%.2f (E2=%.2f, E4=%.2f, E2/E4=%.2f)" %
              (lam, E_d, E2_d, E4_d, v_d))

    E0, E2_0, E4_0 = compute_energy(n, h)
    v0 = E2_0 / (E4_0 + 1e-30)
    print("  Initial: E = %.2f (E2 = %.2f, E4 = %.2f, E2/E4 = %.2f)" %
          (E0, E2_0, E4_0, v0))

    H0 = compute_hopf_charge(n, h)
    print("  Hopf charge: H = %.4f" % H0)
    print()

    if not history['sweep']:
        history['sweep'].append(0)
        history['E'].append(E0)
        history['E2'].append(E2_0)
        history['E4'].append(E4_0)
        history['virial'].append(v0)
        history['acc'].append(1.0)

    rng = np.random.default_rng(42)
    delta = args.delta
    T = args.temp
    best_E = E0
    best_n = n.copy()

    print_every = max(1, args.n_sweeps // 40)
    hopf_every = max(10, args.n_sweeps // 10)
    ckpt_every = max(50, args.n_sweeps // 4)

    print("  %6s  %10s  %8s  %8s  %7s  %6s  %8s  %6s" %
          ("Sweep", "E", "E2", "E4", "E2/E4", "Acc%", "delta", "T"))
    print("  " + "-" * 72)

    t_start = time.time()

    for sweep in range(1, args.n_sweeps + 1):
        nacc, nprop = mc_sweep(n, h, delta, rng, T=T)
        rate = nacc / max(nprop, 1)

        # Adaptive delta: target 20-50% acceptance
        if rate > 0.50:
            delta = min(delta * 1.05, 1.0)
        elif rate < 0.20:
            delta = max(delta * 0.92, 0.005)

        # SA cooling
        if T > 0:
            T *= args.cool

        # Logging
        sw_abs = sweep_offset + sweep
        if sweep % print_every == 0 or sweep <= 5 or sweep == args.n_sweeps:
            E, E2, E4 = compute_energy(n, h)
            v = E2 / (E4 + 1e-30)
            elapsed = time.time() - t_start
            print("  %6d  %10.2f  %8.2f  %8.2f  %7.3f  %5.1f%%  %8.5f  %.1e  [%.0fs]" %
                  (sw_abs, E, E2, E4, v, rate * 100, delta, T, elapsed))
            history['sweep'].append(sw_abs)
            history['E'].append(E)
            history['E2'].append(E2)
            history['E4'].append(E4)
            history['virial'].append(v)
            history['acc'].append(rate)

            if E < best_E:
                best_E = E
                best_n = n.copy()

        # Hopf charge check + topology guard
        if sweep % hopf_every == 0:
            H = compute_hopf_charge(n, h)
            print("    -- Hopf charge: H = %.4f" % H)
            if args.htol > 0 and abs(H) < args.htol:
                print("    ** Topology guard: |H|=%.3f < %.3f, stopping." %
                      (abs(H), args.htol))
                break
            if abs(H) < 0.5:
                print("    ** WARNING: topology may be destroyed (|H| < 0.5)")

        # Periodic Derrick rescaling
        if args.derrick > 0 and sweep % args.derrick == 0:
            E_pre, E2_pre, E4_pre = compute_energy(n, h)
            v_pre = E2_pre / (E4_pre + 1e-30)
            if v_pre > 1.5 or v_pre < 0.67:
                n, lam, E_d, E2_d, E4_d = derrick_rescale(n, h, x1d)
                v_d = E2_d / (E4_d + 1e-30)
                print("    -- Derrick lam=%.3f: E %.2f->%.2f, virial %.2f->%.2f" %
                      (lam, E_pre, E_d, v_pre, v_d))

        # Checkpoint
        if sweep % ckpt_every == 0:
            ckpt_path = os.path.join(OUTDIR, 'fn_mc_checkpoint.npz')
            save_checkpoint(n, history, ckpt_path)

    elapsed = time.time() - t_start
    print()
    print("  Done in %.1f s (%.1f min)" % (elapsed, elapsed / 60))

    # Final diagnostics
    E_f, E2_f, E4_f = compute_energy(n, h)
    H_f = compute_hopf_charge(n, h)
    print()
    print("=" * 60)
    print("  FINAL: E = %.2f (E2=%.2f, E4=%.2f, E2/E4=%.3f)" %
          (E_f, E2_f, E4_f, E2_f / (E4_f + 1e-30)))
    print("  Hopf charge: H = %.4f (initial: %.4f)" % (H_f, H0))
    print("  BS target: 192.5,  ratio E/E_BS = %.4f" % (E_f / 192.5))
    print("  Best E encountered: %.2f" % best_E)
    print("=" * 60)

    print()
    print("  Generating plots...")
    plot_results(history, n, h, x1d, OUTDIR)

    # Save final state
    final_path = os.path.join(OUTDIR, 'fn_mc_final.npz')
    save_checkpoint(n, history, final_path)
    print("  Done.")

    return n, history


if __name__ == '__main__':
    main()
