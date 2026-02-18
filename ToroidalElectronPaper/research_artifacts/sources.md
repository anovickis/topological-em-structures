# Research Artifacts — Source Log

This file logs all fetched URLs, papers, and external sources consulted during research for the Toroidal Electron project.

---

## 2026-02-17: Multi-Linking / SU(3) / Quark Substructure Research

### Topic: SU(3) Cho-Faddeev-Niemi Decomposition

| Source | Key Finding |
|--------|-------------|
| Cho, Y.M. (1980). "Restricted gauge theory." Phys. Rev. D 21, 1080. | Original CFN decomposition for SU(2). |
| Cho, Y.M. & Pak, D.G. (2002). "Monopole condensation in SU(3)." Phys. Rev. D 65, 074027. | Extended CFN to SU(3): target space is flag manifold $F_2 = \text{SU}(3)/[\text{U}(1) \times \text{U}(1)]$, two independent topological charges from $\pi_2(F_2) = \mathbb{Z} \oplus \mathbb{Z}$. |
| Faddeev, L. & Niemi, A.J. (1997). "Stable knot-like structures in classical field theory." Nature 387, 58. | Original FN soliton paper. Target space $S^2 = \text{SU}(2)/\text{U}(1)$, admits Hopf solitons. |
| Shabanov, S. (2000). "An effective action for monopoles and knot solitons in Yang-Mills theory." Phys. Lett. B 458, 322. | Effective FN Lagrangian from integrating out off-diagonal modes; discusses non-perturbative matching. |

### Topic: CP2 Homotopy and Higher Target Spaces

| Source | Key Finding |
|--------|-------------|
| Hatcher, A. (2002). *Algebraic Topology*. Cambridge. (Ch. 4, homotopy groups table) | $\pi_3(\mathbb{CP}^2) = 0$ --- NO Hopf-like 3D solitons for CP2 target! Critical negative result. |
| Eichenherr, H. & Forger, M. (1981). "More about nonlinear sigma models on symmetric spaces." Nucl. Phys. B 164, 528. | Classification of sigma models on coset spaces; flag manifold sigma models have richer topology than CP^N. |

### Topic: Skyrmion Quark Substructure

| Source | Key Finding |
|--------|-------------|
| Atiyah, M.F. & Manton, N.S. (1989). "Skyrmions from Instantons." Phys. Lett. B 222, 438. | B=1 skyrmion from instanton holonomy reveals 3-quark internal structure at vertices of equilateral triangle. |
| Manton, N.S. (1987). "Geometry of Skyrmions." Comm. Math. Phys. 111, 469. | Skyrmion moduli space geometry; B=1 has natural Z3 decomposition. |
| Houghton, C.J., Manton, N.S. & Sutcliffe, P.M. (1998). "Rational maps, monopoles and skyrmions." Nucl. Phys. B 510, 507. | Rational map ansatz for multi-skyrmions; reveals polyhedral (including triangular) symmetry of B=1. |
| Battye, R.A. & Sutcliffe, P.M. (1998). "Knots as Stable Soliton Solutions in a Three-Dimensional Classical Field Theory." PRL 81, 4798. | Numerical Hopf soliton profiles. |H|=1 is fat torus, |H|=2 is bent torus, higher charges form torus knots. |

### Topic: Fractional Topology and Confinement

| Source | Key Finding |
|--------|-------------|
| 't Hooft, G. (1981). "Topology of the gauge condition and new confinement phases in non-Abelian gauge theories." Nucl. Phys. B 190, 455. | Center vortex model of confinement: Z_N center of SU(N), Wilson loop area law from vortex linking numbers. |
| González-Arroyo, A. & Montero, A. (1998). "Selfdual vortex-like configurations in SU(2) Yang-Mills theory." Phys. Lett. B 442, 273. | Fractional instantons on twisted torus carry topological charge Q=1/N. For SU(3), fractional charge = 1/3! |
| Duan, Y.S., Liu, X. & Zhang, P.M. (2003). "Decomposition of the Hopf invariant." J. Phys. A 36, 563. | Hopf invariant decomposes as: $H = \sum_i \text{SL}(C_i) + \sum_{i<j} \text{Lk}(C_i, C_j)$ (self-linking + pairwise linking of preimage knots). Key for sector decomposition! |
| Greensite, J. (2011). *An Introduction to the Confinement Problem*. Springer. | Review of center vortex confinement, dual superconductor model, and Abelian projection. |

### Topic: Weyl Group and Z3 Symmetry

| Source | Key Finding |
|--------|-------------|
| Balachandran, A.P. et al. (1983). "Monopole topology and the problem of color." PRL 50, 1553. | Topological origin of color from homotopy of SU(3)/Z3; explains why quarks are confined. |
| Preskill, J. (1984). "Magnetic monopoles." Ann. Rev. Nucl. Part. Sci. 34, 461. | Review of monopoles in non-Abelian gauge theories; Z_N structure and fractional charges. |

---

## Instructions

When fetching new web content or consulting external sources during research:
1. Add a new dated section to this file
2. Log the URL or paper reference
3. Note the key finding relevant to the project
4. Save any substantial fetched content as a separate `.md` file in this folder
