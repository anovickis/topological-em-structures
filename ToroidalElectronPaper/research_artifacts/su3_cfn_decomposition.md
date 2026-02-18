# SU(3) Cho-Faddeev-Niemi Decomposition — Research Notes

Tags: #physics #topology #SU3 #CFN
Created: 2026-02-17

---

## Summary

The Cho-Faddeev-Niemi (CFN) decomposition extends from SU(2) to SU(3), yielding a richer topological structure relevant to the multi-linking quark model.

## SU(2) CFN (Review)

For SU(2) gauge field $A_\mu^a$:

$$A_\mu = C_\mu \hat{n} + \frac{1}{g} \hat{n} \times \partial_\mu \hat{n} + W_\mu \tag{1}$$

- $C_\mu$: Abelian (photon-like) component
- $\hat{n} \times \partial_\mu \hat{n}$: topological (valence) potential — carries the soliton
- $W_\mu$: off-diagonal massive modes (W bosons)
- Target space: $S^2 = \text{SU}(2)/\text{U}(1)$
- Topology: $\pi_3(S^2) = \mathbb{Z}$ (Hopf invariant)

## SU(3) CFN Extension (Cho-Pak 2002)

For SU(3), the decomposition involves **two** unit vector fields (color directions):

$$A_\mu = \sum_{i=1}^{2} C_\mu^{(i)} \hat{n}_i + \text{valence terms} + W_\mu \tag{2}$$

### Target Space

The relevant coset space is the **flag manifold**:

$$F_2 = \frac{\text{SU}(3)}{\text{U}(1) \times \text{U}(1)}$$

This is a 6-dimensional manifold (dim SU(3) = 8, dim U(1)xU(1) = 2).

**Crucially different from $\mathbb{CP}^2$:** While $\mathbb{CP}^2 = \text{SU}(3)/[\text{U}(2)]$ has $\pi_3(\mathbb{CP}^2) = 0$ (no 3D Hopf-like solitons!), the flag manifold has:

$$\pi_2(F_2) = \mathbb{Z} \oplus \mathbb{Z}$$

This means there are **two independent topological charges** — exactly what we need for the two types of quarks (up-type and down-type charges).

### Why $\mathbb{CP}^2$ Fails

- $\mathbb{CP}^2 = \text{SU}(3)/\text{U}(2)$ is a 4-dimensional manifold
- $\pi_3(\mathbb{CP}^2) = 0$ — there are NO Hopf-like solitons in 3D with $\mathbb{CP}^2$ target
- $\pi_2(\mathbb{CP}^2) = \mathbb{Z}$ — only one topological charge (not enough for quarks)
- The naive generalization "replace $S^2$ with $\mathbb{CP}^2$" does NOT work

### Why the Flag Manifold Works

- $F_2 = \text{SU}(3)/[\text{U}(1) \times \text{U}(1)]$ is 6-dimensional
- $\pi_2(F_2) = \mathbb{Z} \oplus \mathbb{Z}$ — two independent charges $(n_1, n_2)$
- The maximal torus $\text{U}(1) \times \text{U}(1)$ corresponds to the two Cartan generators of SU(3)
- The Weyl group $S_3$ (permutation group of 3 elements) acts on the charges — this is the origin of the $\mathbb{Z}_3$ color symmetry!

## Connection to Multi-Linking Model

In the flag manifold sigma model:
- The two topological charges $(n_1, n_2)$ map to the two independent quark charges
- The Weyl group $S_3$ permutes color labels (R, G, B)
- A "baryon" configuration has charges summing to give the correct total
- The flag manifold naturally decomposes into three sectors related by the $\mathbb{Z}_3$ center of SU(3)

### Charge Identification

If we identify:
- Up-type quark: $(n_1, n_2) = (1, 0)$ → charge $+2/3$
- Down-type quark: $(n_1, n_2) = (0, -1)$ → charge $-1/3$

Then:
- Proton (uud): $(1,0) + (1,0) + (0,-1) = (2,-1)$ → total charge $+1$
- Neutron (udd): $(1,0) + (0,-1) + (0,-1) = (1,-2)$ → total charge $0$

The integer charges of the flag manifold sigma model give fractional physical charges through a rescaling by $1/3$.

## Key References

1. Cho, Y.M. & Pak, D.G. (2002). Phys. Rev. D 65, 074027.
2. Faddeev, L. & Niemi, A.J. (1999). "Partially dual variables in SU(2) Yang-Mills theory." Phys. Rev. Lett. 82, 1624.
3. Shabanov, S. (2000). Phys. Lett. B 458, 322.
4. Eichenherr, H. & Forger, M. (1981). Nucl. Phys. B 164, 528.
