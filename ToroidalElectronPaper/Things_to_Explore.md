# Things to Explore --- Toroidal Electron Project

Tags: #physics #topology #electron #future-work
Created: 2026-02-17
Updated: 2026-02-17

> [!summary] TL;DR
> Running list of ideas, directions, and questions worth investigating that don't yet fit into the paper or companion documents. Organized by topic with rough priority/feasibility assessment.

---

## Vacuum Structure and Topological Fluctuations

### The Vacuum as Topological Ground State

**Idea:** In the Faddeev-Niemi model, the vacuum is the $H = 0$ sector where $\mathbf{n} = \text{const}$ everywhere. But quantum fluctuations of the $\mathbf{n}$-field include transient topological excitations. The question: does the topological structure of the $\mathbf{n}$-field vacuum differ from the standard QFT vacuum?

**Key physics:**
- The FN model has a one-loop effective action around $\mathbf{n} = \text{const}$. The spectrum includes massless Goldstone-like modes (spin waves of $\mathbf{n}$ on $S^2$) and the soliton sector
- Virtual $e^+e^-$ pairs correspond to transient $H = +1, -1$ fluctuations of the $\mathbf{n}$-field
- The Vakulenko-Kapitanski bound $E \geq C_{VK}|H|^{3/4}$ sets the minimum energy for any $|H| = 1$ fluctuation --- matching the threshold $2m_e c^2$
- The vacuum polarization (Uehling potential, running of $\alpha$) should emerge from integrating over topological fluctuations

**What to calculate:**
1. One-loop effective action of the FN model around the trivial vacuum $\mathbf{n} = \hat{z}$
2. Whether the topological constraint (integer $H$) suppresses the vacuum energy relative to standard QFT
3. The vacuum polarization contribution from $H = \pm 1$ fluctuations vs standard QED virtual pairs

**Feasibility:** Medium. The one-loop calculation is well-defined but technically demanding (functional determinant on $S^2$ target space).

**Priority:** High --- this connects directly to the bootstrap hypothesis in Appendix A.4.

### Casimir Effect from Topological Perspective

**Idea:** The Casimir effect arises because conducting plates restrict which field modes can exist between them. In the FN model, this restriction applies to the $\mathbf{n}$-field: fewer allowed topological fluctuations between the plates → lower vacuum energy → attractive force.

**Key physics:**
- Standard Casimir calculation: sum over modes of the EM field with boundary conditions → $F/A = -\pi^2 \hbar c / (240 d^4)$
- In the FN model, the $\mathbf{n}$-field has the same mode spectrum as the EM field (they're related by the CFN decomposition) plus additional nonlinear corrections
- The Skyrme-like term in the FN Lagrangian contributes a correction to the Casimir energy at fourth order in $1/d$
- This correction would be suppressed by $(\bar{\lambda}_C / d)^2$ relative to the standard result

**What to calculate:**
1. The FN model's Casimir energy between parallel plates with Dirichlet b.c. for $\mathbf{n}$
2. Whether the nonlinear term gives a measurable correction at experimentally accessible plate separations

**Feasibility:** Medium-high. The linear (sigma model) part reproduces standard Casimir; the Skyrme correction is perturbative and calculable.

**Priority:** Medium --- would be a nice consistency check but unlikely to produce measurable deviations.

### Cosmological Constant Problem

**Idea:** The QFT vacuum energy is $\sim 10^{120}$ times larger than the observed cosmological constant. Could the topological constraint on the $\mathbf{n}$-field (Hopf invariant must be integer) suppress the vacuum energy?

**Speculative argument:**
- In standard QFT, ALL field modes contribute to vacuum energy
- In the FN model, the $\mathbf{n}$-field lives on $S^2$, which is compact --- this automatically regulates the UV divergence differently from a flat target space
- The Vakulenko-Kapitanski bound means that topological excitations have a minimum energy gap --- below this, only "smooth" fluctuations contribute
- If the smooth fluctuations partially cancel (as in supersymmetry but from topology), the vacuum energy could be suppressed

**What to calculate:**
1. The vacuum energy of the FN sigma model on $S^2$ vs flat target space $\mathbb{R}^2$
2. Whether compactness of the target space provides a natural UV cutoff
3. Comparison with the observed $\Lambda \sim (2.3 \text{ meV})^4$

**Feasibility:** Low --- this is essentially the cosmological constant problem, one of the deepest unsolved problems in physics. But the FN model provides a new angle.

**Priority:** Low (fascinating but extremely difficult).

---

## Multi-Linking Extensions

### Flag Manifold Sigma Model Numerics

**Idea:** Compute the minimum-energy soliton of the flag manifold sigma model $F_2 = \text{SU}(3)/[\text{U}(1) \times \text{U}(1)]$ with a Skyrme-like stabilizer.

**What to compute:**
1. Discretize $F_2$ on a 3D grid (6 real fields with constraints)
2. Initialize with a generalized Hopf map adapted to the flag manifold
3. Gradient flow with topology monitoring (two charges $(n_1, n_2)$)
4. Check for three-fold internal structure in the energy density

**Feasibility:** High (given existing sim_fn_soliton_3d.py infrastructure). Main challenge: parameterizing the flag manifold and defining the Skyrme term.

**Priority:** Very high --- this is THE decisive numerical test for the multi-linking quark model.

### Meson Sector and String Tension

**Idea:** Compute the energy of a "stretched" soliton where two sectors are pulled apart (quark-antiquark separation) and verify linear confinement.

**Feasibility:** Medium. Requires constrained energy minimization with sector positions fixed.

**Priority:** High --- would directly test the confinement mechanism.

### Neutrino in the Multi-Linking Framework

**Idea:** The neutrino has $H = 0$ in the multi-linking picture. How is it different from the vacuum? Possibilities:
- Neutrino is a topological excitation of the $\mathbf{n}$-field with $H = 0$ but non-trivial knot type (like the DM paper's trefoil, but lighter)
- Neutrino is a soliton in a different sector (e.g., the "angular momentum" sector rather than the "linking" sector)
- Neutrino mass arises from non-zero $\pi_1$ of the configuration space (analogous to FR mechanism giving spin-1/2)

**Feasibility:** Speculative. No clear mathematical framework yet.

**Priority:** Medium --- important for completeness but may require new ideas.

---

## Electron Structure Calculations

### Topology-Preserving Energy Minimization

**Status:** Basin-hopping MC implemented (sim_fn_soliton_3d.py), launched on remote server. Early basins struggling with virial ratio.

**Next steps:**
1. Analyze basin-hopping results when remote computation completes
2. If best basin reaches $E < 250$: refine with arrested Newton flow
3. If all basins fail: try different initial conditions (not just Hopf maps)
4. Target: $E = 192.5$ (Battye-Sutcliffe minimum)

### C2 from Equilibrium Soliton

Once the soliton converges to the BS minimum:
1. Extract aspect ratio $A$ from the equilibrium profile
2. Compute $C_2 = -\pi^2/(4A^2)$ from Eq. 13.20
3. Compare with QED's $C_2 = -0.3285$
4. If agreement: this is a parameter-free prediction of a quantity known to 4 significant figures

### Deriving QED from Soliton Dynamics

**The big prize:** Show that QED's perturbative expansion emerges from the collective coordinate quantization of the Hopf soliton. This would require:
1. Compute the soliton's propagator in the collective coordinate framework
2. Show it reproduces the Dirac propagator at low energies
3. Derive the photon-electron vertex from soliton-photon scattering
4. Recover Feynman diagrams as a systematic expansion

**Feasibility:** Very low (this is essentially solving the theory). But partial results (like deriving the vertex function to lowest order) might be achievable.

**Priority:** Highest possible if achievable --- would prove the theory.

---

## Connections to Other Physics

### Unruh Effect Inside the Soliton

Already noted in Appendix A.7. The circular acceleration $a = c^2/R$ gives Unruh temperature $T_U = m_e c^2/(2\pi k_B)$. Open question: does this modify the soliton's effective energy?

### Gravitational Properties

Does the soliton's stress-energy tensor produce the correct gravitational mass? The FN stress-energy tensor is known; computing the ADM mass from the soliton's metric perturbation would verify $m_{\text{grav}} = m_{\text{inertial}}$ (equivalence principle).

### Black Hole Information Paradox

If particles are topological solitons, what happens to the topology when they fall into a black hole? The Hopf invariant is conserved under smooth evolution --- but the black hole singularity is not smooth. Does the topology survive Hawking evaporation?

---

## Quick Priority Summary

| Topic | Priority | Feasibility | Placement if successful |
|-------|----------|-------------|------------------------|
| Flag manifold soliton numerics | Very high | High | Paper + Multi-Linking doc |
| C2 from equilibrium soliton | Very high | Medium | Paper §13.5, §14.2 |
| Vacuum one-loop effective action | High | Medium | Paper Appendix A |
| Meson string tension | High | Medium | Multi-Linking doc |
| Casimir correction | Medium | Medium-high | Paper Appendix A |
| Neutrino mechanism | Medium | Low | Multi-Linking doc |
| Cosmological constant | Low | Very low | Paper Appendix A |
| Deriving QED | Highest | Very low | New paper |
