# Skyrmion Quark Substructure — Research Notes

Tags: #physics #topology #skyrmions #quarks
Created: 2026-02-17

---

## Summary

The Skyrme model provides a well-established precedent for topological solitons containing quark-like sub-structure. This is the closest analogue to the multi-linking Hopf soliton model.

---

## 1. The Skyrme Model

The Skyrme Lagrangian for pion fields:

$$\mathcal{L} = \frac{f_\pi^2}{16} \text{Tr}(\partial_\mu U \partial^\mu U^\dagger) + \frac{1}{32e^2} \text{Tr}([U^\dagger \partial_\mu U, U^\dagger \partial_\nu U]^2) \tag{1}$$

where $U: \mathbb{R}^3 \to \text{SU}(2)$ with boundary condition $U \to 1$ at infinity.

- Topological charge: Baryon number $B = \frac{1}{24\pi^2} \int \text{Tr}(U^\dagger dU \wedge U^\dagger dU \wedge U^\dagger dU)$
- $B \in \mathbb{Z}$ (winding number of $S^3 \to S^3$)
- Fermionic quantization via Finkelstein-Rubinstein (same as Hopf soliton!)

---

## 2. Atiyah-Manton Construction: Instantons Reveal Quarks

**Key insight (Atiyah & Manton, 1989):**

The B=1 skyrmion can be constructed from a Yang-Mills instanton via holonomy:

$$U(\mathbf{x}) = \mathcal{P} \exp\left(-\int_{-\infty}^{\infty} A_4(\mathbf{x}, x_4)\, dx_4\right) \tag{2}$$

When the instanton is decomposed into constituent instantons (monopole-instanton constituents on $\mathbb{R}^3 \times S^1$), the B=1 skyrmion decomposes into **three "quarks"**:

- The three quarks sit at the vertices of an equilateral triangle
- Each carries baryon number $B = 1/3$
- The triangle has $\mathbb{Z}_3$ symmetry
- The quarks are NOT independent topological objects — they are features of a single $B=1$ configuration

### Parallels to Hopf Soliton

| Skyrme model | Hopf soliton model |
|---|---|
| $B = 1$ skyrmion | $H = 1$ Hopf soliton |
| SU(2) target, $\pi_3(\text{SU}(2)) = \mathbb{Z}$ | $S^2$ target, $\pi_3(S^2) = \mathbb{Z}$ |
| 3 quarks at triangle vertices | 3 sectors of torus |
| Each quark: $B = 1/3$ | Each sector: $H_{\text{eff}} = 1/3$ or $2/3$ |
| Color = triangle vertex label | Color = torus sector label |
| Confinement = can't separate vertices | Confinement = can't cut torus |
| FR fermionic quantization | FR fermionic quantization |

---

## 3. Rational Map Ansatz

**Houghton, Manton & Sutcliffe (1998):**

For $B = N$ skyrmions, the angular dependence is captured by a rational map $R: S^2 \to S^2$ of degree $N$:

$$U(r, z) = \exp(if(r) \hat{\mathbf{n}}(z) \cdot \boldsymbol{\tau}) \tag{3}$$

where $z$ is the Riemann sphere coordinate and $\hat{\mathbf{n}}$ is determined by the rational map.

For $B = 1$: $R(z) = z$ (identity map)
- The energy density is spherically symmetric ("hedgehog")
- But the **color orientation** varies: different points on the sphere have different "color" directions

For $B = 2$: $R(z) = z^2$ (axially symmetric "torus")

For $B = 3$: $R(z) = \frac{z^3 - \sqrt{3}iz}{iz^3 - \sqrt{3}}$ (tetrahedral symmetry)

### Key Point for Hopf Solitons

The rational map reveals that even the "simple" $B=1$ skyrmion has internal structure — the color direction rotates as you go around the soliton. For Hopf solitons, the analogous statement is that the $\mathbf{n}$-field direction rotates around the torus, and this rotation can be decomposed into sectors.

---

## 4. SU(3) Skyrmions and Strangeness

When the Skyrme model is extended from SU(2) to SU(3) flavor:

$$U: \mathbb{R}^3 \to \text{SU}(3) \tag{4}$$

The embedding $\text{SU}(2) \hookrightarrow \text{SU}(3)$ gives:

- The $B=1$ SU(3) skyrmion is an SU(2) hedgehog embedded in the upper-left $2 \times 2$ block
- Quantization over the SU(3)/SU(2) coset (kaon zero-mode space) gives the baryon octet and decuplet
- Strangeness arises from excitations in the "strange" direction of flavor space

### Relevance

This shows that the quark flavor structure (u, d, s) emerges from the quantization of a single soliton's orientation in a larger symmetry space. The same mechanism could apply to Hopf solitons in the flag manifold: different "flavors" correspond to different orientations in the SU(3)/[U(1)×U(1)] target space.

---

## 5. Witten's Identification: Skyrmion = Baryon at Large $N_c$

**Witten (1979):** In the large-$N_c$ limit of QCD:

- Mesons are weakly coupled ($g_{\text{eff}} \sim 1/\sqrt{N_c}$)
- Baryons are heavy ($m_B \sim N_c$)
- Baryons are solitons in the meson field

The Skyrme model is the $N_c \to \infty$ effective theory of QCD, and skyrmions ARE baryons. This is not an analogy — it's a derivation.

### For Hopf Solitons

If the Faddeev-Niemi model can be derived from SU(3) gauge theory in a similar limit:
- Hopf solitons in the SU(3) FN model would BE hadrons
- The sector decomposition would BE the quark structure
- Confinement would be topological by construction

---

## 6. Open Questions

1. **Instanton construction for Hopf solitons?** Is there an analogue of the Atiyah-Manton construction that relates Yang-Mills instantons to Hopf solitons? This would provide the formal link between the gauge theory and the soliton.

2. **Large-$N_c$ limit for FN model?** Does the FN model emerge as the effective theory of SU(N) Yang-Mills at large $N_c$? Faddeev & Niemi argued "yes" but a rigorous derivation is lacking.

3. **Numerical Z₃ structure?** Can the $H=1$ Hopf soliton be shown numerically to have three-fold internal structure (e.g., by computing the energy density in $\mathbf{n}$-field sectors)?

---

## Key References

1. Atiyah, M.F. & Manton, N.S. (1989). Phys. Lett. B 222, 438.
2. Manton, N.S. (1987). Comm. Math. Phys. 111, 469.
3. Houghton, C.J., Manton, N.S. & Sutcliffe, P.M. (1998). Nucl. Phys. B 510, 507.
4. Witten, E. (1979). "Baryons in the 1/N expansion." Nucl. Phys. B 160, 57.
5. Guadagnini, E. (1984). "Baryons as solitons and mass formulae." Nucl. Phys. B 236, 35.
6. Manton, N.S. & Sutcliffe, P.M. (2004). *Topological Solitons*. Cambridge. Ch. 9.
