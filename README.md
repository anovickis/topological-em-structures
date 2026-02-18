# Topological Electromagnetic Structures in Particle Physics

Research papers and computational tools exploring topological soliton models of fundamental particles.

## Papers

### 1. The Toroidal Electron
**"The Toroidal Electron: A Unified Geometric Theory of Electromagnetic Structure, Mass, and the Fine Structure Constant"**

The electron is proposed as a topological Hopf soliton (linking number $H = \pm 1$) in the Faddeev-Niemi nonlinear sigma model, derived from SU(2) gauge theory via the Cho-Faddeev-Niemi decomposition. The framework produces: charge quantization from topology, spin-1/2 from the Finkelstein-Rubinstein mechanism, the de Broglie wavelength from Lorentz-boosted internal circulation, and a point-like electric form factor from a topological Ward identity.

- **Paper:** [`ToroidalElectronPaper/Toroidal_Electron_Full_Paper.md`](ToroidalElectronPaper/Toroidal_Electron_Full_Paper.md)
- **FAQ:** [`ToroidalElectronPaper/FAQ.md`](ToroidalElectronPaper/FAQ.md)
- **Multi-linking extension (quarks):** [`ToroidalElectronPaper/Multi-Linking_Soliton_Spectrum.md`](ToroidalElectronPaper/Multi-Linking_Soliton_Spectrum.md)

### 2. Dark Matter as Topological EM Structures
**"Dark Matter as Topological Electromagnetic Structures: Mathematical Framework for H = 0 Stable Configurations"**

Extends the toroidal electron model to propose that dark matter consists of stable, uncharged ($H = 0$) knotted solitons (trefoils, figure-8 knots) in the Faddeev-Niemi model. Predicts MeV-scale dark matter with specific gamma-ray signatures testable by the COSI mission (2027).

- **Paper:** [`DarkMatterPaper/Dark_Matter_Topological_EM_Structures.md`](DarkMatterPaper/Dark_Matter_Topological_EM_Structures.md)
- **Critical Review:** [`DarkMatterPaper/Dark_Matter_Critical_Review.md`](DarkMatterPaper/Dark_Matter_Critical_Review.md)

## Computational Tools

All Python scripts use `matplotlib` with headless (`Agg`) backend.

### Toroidal Electron
| Script | Purpose |
|--------|---------|
| `sim_fn_soliton_c2.py` | FN soliton gradient flow on cylindrical grid, C2 computation |
| `sim_fn_soliton_3d.py` | 3D FN soliton solver with basin-hopping Monte Carlo |
| `sim_mass_formula_search.py` | Exhaustive search confirming {21,15} uniqueness |
| `sim_coupling_threshold.py` | One-loop threshold matching for FN couplings |
| `sim_lattice_su2_matching.py` | Lattice SU(2)+Higgs Monte Carlo |
| `sim_higgs_backreaction.py` | Higgs VEV suppression and self-consistency |

### Dark Matter
| Script | Purpose |
|--------|---------|
| `sim_faddeev_energy.py` | FN energy minimization for unknot/trefoil |
| `sim_knot_ropelength.py` | Mass predictions from ropelength scaling |
| `sim_relic_abundance.py` | Boltzmann freeze-out calculation |
| `sim_cross_sections.py` | Polarizability and annihilation cross-sections |
| `sim_gamma_ray_flux.py` | Gamma-ray flux from DM annihilation |
| `sim_neff_bbn.py` | BBN $N_{\text{eff}}$ constraint analysis |
| `sim_sensitivity_comparison.py` | Experimental sensitivity comparison |

## Key Results

| Result | Value | Paper |
|--------|-------|-------|
| Soliton energy match | $\kappa_2 \sim \alpha\hbar c$ (15%) | Toroidal Electron |
| $C_2$ anomalous magnetic moment | $\approx -0.33$ vs QED's $-0.3285$ | Toroidal Electron |
| Mass formula accuracy | 0.008% with {21, 15} | Toroidal Electron |
| Best dark matter mass (trefoil) | $2.0^{+0.7}_{-0.8}$ MeV | Dark Matter |
| Dark matter free-streaming | $< 10^{-2}$ Mpc (cold) | Dark Matter |

## Author

**Alexander Novickis** - alex.novickis@gmail.com

## License

These papers and code are provided for academic and research purposes. Please cite the relevant paper if you use this work.
