# Exciton transition matrix elements in BerkeleyGW `absorption`

What is written in columns 3–4 of `eigenvalues.dat` (velocity operator) and
`eigenvalues_b{1,2,3}.dat` (momentum operator) for a TDA, complex-build
BSE diagonalization run.

Example runs:
- velocity: `/p/lustre5/simmerman1/bulk_silicon/test_ideas/14-pbe0_scf_scfq/01-nk8_shifted/06-abs-eh-vel`
- momentum: `/p/lustre5/simmerman1/bulk_silicon/test_ideas/14-pbe0_scf_scfq/01-nk8_shifted/05-abs-eh-mom`

## File columns

`write_eigenvalues` (`BSE/absp_io.f90:120`) writes, for each exciton $S$:

| column | quantity |
|---|---|
| 1 | $\Omega_S$ (eV) |
| 2 | $\lvert d_S\rvert^2$ |
| 3 | $\mathrm{Re}\, d_S$ |
| 4 | $\mathrm{Im}\, d_S$ |

with $d_S$ computed in `BSE/diag.f90:710`:

$$
d_S \;=\; \frac{1}{\sqrt{N_s}}\sum_{c v \mathbf{k}\sigma} A^S_{cv\mathbf{k}\sigma}\;\big[s_{cv\mathbf{k}\sigma}\big]^*
$$

- $A^S_{cv\mathbf k}$ is the exciton envelope function, normalized as
  $\sum_{cv\mathbf k}\lvert A^S_{cv\mathbf k}\rvert^2 = 1$ (no $1/\sqrt{N_k}$).
- $N_s$ = `nspin` (the factor makes singlet results agree between `nspin=1`
  and `nspin=2`; see Rohlfing & Louie, PRB 62, 4927 (2000)).
- $s_{cv\mathbf k}$ is the single-particle transition matrix element from
  `Common/mtxel_optical.f90`; this is the only thing that differs between the
  two approximations.

The spectrum is (`BSE/absp.f90:65`)

$$
\varepsilon_2(\omega) = \frac{16\pi^2}{V\, n_{\rm spinor}} \sum_S \lvert d_S\rvert^2\, \delta(\omega-\Omega_S),
$$

where $V$ = header `vol` $= N_k\,\Omega_{\rm cell}$, which cancels the
$N_k$ scaling of $\lvert d_S\rvert^2$ implied by the normalization of $A^S$.

Units: Rydberg atomic units ($\hbar = 1$, $m = 1/2$), so $\hat{\mathbf p} = -i\nabla$
and the local velocity operator is $\hat{\mathbf v} = 2\hat{\mathbf p}$. $d_S$ has
units of length (bohr).

## Momentum operator (`use_momentum`)

`mtxel_m` with `divide_energy=.true.`:

$$
s^{\rm mom}_{cv\mathbf{k}} = \frac{\langle \psi_{c\mathbf{k}}|\,\hat{\mathbf e}\cdot 2\hat{\mathbf p}\,|\psi_{v\mathbf{k}}\rangle}{\varepsilon^{\rm MF}_{c\mathbf k}-\varepsilon^{\rm MF}_{v\mathbf k}}
= \frac{\sum_{\mathbf G} c^*_{c\mathbf k}(\mathbf G)\,c_{v\mathbf k}(\mathbf G)\;2\,\hat{\mathbf e}\cdot\mathbf G}{\varepsilon^{\rm MF}_{c\mathbf k}-\varepsilon^{\rm MF}_{v\mathbf k}}
$$

(the $\mathbf k$ part of $\mathbf k+\mathbf G$ drops out by orthogonality for $c\neq v$), so

$$
d_S^{\rm mom} = \frac{1}{\sqrt{N_s}}\sum_{cv\mathbf k} A^S_{cv\mathbf k}\,
\frac{\langle v\mathbf k|\hat{\mathbf e}\cdot 2\hat{\mathbf p}|c\mathbf k\rangle}{\varepsilon^{\rm MF}_{c\mathbf k}-\varepsilon^{\rm MF}_{v\mathbf k}}
\;\approx\; \frac{-i}{\sqrt{N_s}}\sum_{cv\mathbf k} A^S_{cv\mathbf k}\,\langle v\mathbf k|\hat{\mathbf e}\cdot\mathbf r|c\mathbf k\rangle .
$$

Notes:
- **Direction.** $\hat{\mathbf e} = \hat{\mathbf b}_j$ for file `b`$j$: the
  polarization is a unit vector in reciprocal-lattice crystal coordinates,
  normalized with the `bdot` metric. For the cubic (8-atom, $\Omega \approx 1081$
  bohr³) Si cell used here, $b_1, b_2, b_3 = \hat x, \hat y, \hat z$.
- **Energy denominator.** Mean-field energies (`eqp%eclda`, `eqp%evlda`) from
  `WFN_fi`, *not* quasiparticle energies and *not* $\Omega_S$. If the gap is below
  `TOL_Degeneracy`, $s$ is set to zero.
- **Missing commutator.** The "$\approx \mathbf r$" step uses
  $\hat{\mathbf v} = i[H,\mathbf r] = 2\hat{\mathbf p}$, valid only for a local potential. The
  commutator $[V_{\rm NL},\mathbf r]$ (nonlocal pseudopotential projectors and, for
  PBE0, the nonlocal Fock exchange) is omitted.

## Velocity operator (`use_velocity`)

`mtxel_v`:

$$
s^{\rm vel}_{cv\mathbf{k}} = \frac{\langle \psi_{c\mathbf{k}}|e^{-i\mathbf q\cdot\mathbf r}|\psi_{v\mathbf{k+q}}\rangle}{|\mathbf q|}
= \frac{\langle u_{c\mathbf k}|u_{v\mathbf k+\mathbf q}\rangle}{|\mathbf q|}
= \frac{1}{|\mathbf q|}\sum_{\mathbf G}c^*_{c\mathbf k}(\mathbf G)\,c_{v\mathbf k+\mathbf q}(\mathbf G),
$$

so

$$
d_S^{\rm vel} = \frac{1}{\sqrt{N_s}}\sum_{cv\mathbf k} A^S_{cv\mathbf k}\,
\frac{\langle \psi_{v\mathbf k+\mathbf q}|e^{i\mathbf q\cdot\mathbf r}|\psi_{c\mathbf k}\rangle}{|\mathbf q|}.
$$

- $\mathbf q$ is the shift between `WFN_fi` and `WFNq_fi`; here $(0,0,10^{-5})$ in
  crystal coordinates, so $\hat{\mathbf e} = \hat{\mathbf q} = \hat{\mathbf b}_3$, and
  $|\mathbf q| = \sqrt{\mathbf q\cdot B\cdot\mathbf q}$ (bohr⁻¹). Only one polarization is
  available.
- $A^S$ is expressed in the basis $\{|c\mathbf k\rangle, |v\mathbf k+\mathbf q\rangle\}$, since the
  valence states come from `WFNq_fi`.

First-order $\mathbf k\cdot\mathbf p$ expansion, with $\partial H_{\mathbf k}/\partial\mathbf k = \hat{\mathbf v}_{\rm full}$:

$$
s^{\rm vel}_{cv\mathbf k} \approx \frac{\langle c\mathbf k|\hat{\mathbf q}\cdot\hat{\mathbf v}_{\rm full}|v\mathbf k\rangle}{\varepsilon^{\rm MF}_{v\mathbf k}-\varepsilon^{\rm MF}_{c\mathbf k}}
= -\,\hat{\mathbf q}\cdot\langle c\mathbf k|\,i\mathbf r\,|v\mathbf k\rangle + O(q).
$$

- $\hat{\mathbf v}_{\rm full} = i[H,\mathbf r]$ with the full mean-field Hamiltonian,
  **including** $V_{\rm NL}$ and the PBE0 exchange, because `WFNq_fi` was generated
  with the same $H$. This is the main physical difference from the momentum
  form, and is consistent with the larger f-sum rule in the example runs
  (≈1.04 velocity vs. ≈0.77 momentum, non-interacting).
- The sign is opposite to the momentum form (irrelevant for $\lvert d_S\rvert^2$);
  there is an $O(q)$ finite-difference error.

## Caveats when reading columns 3–4

1. **Re and Im separately are not physical.** The eigensolver fixes the overall
   phase of $A^S$ arbitrarily, and Bloch states carry arbitrary gauge phases.
   Only $\lvert d_S\rvert^2$ is meaningful.
2. **Degenerate excitons.** For a degenerate multiplet (e.g. the first four
   states at 1.87436 eV in the example), only
   $\sum_{S\in\text{multiplet}}\lvert d_S\rvert^2$ is invariant; individual values depend
   on the basis the diagonalizer picked within the degenerate subspace.
3. **Scaling.** $d_S$ is a length (bohr) and scales as $\sqrt{N_k}$.
4. **`_noeh` files.** `eigenvalues*_noeh.dat` list the bare $s_{cv\mathbf k}$
   (`write_eigenvalues_noeh`, `BSE/absp_io.f90:189`). Velocity and momentum
   files use the same definitions as above but with different $\hat{\mathbf e}$, so rows
   need not match one-to-one.

## Source references

- `BSE/diag.f90:706-735` — contraction of $A^S$ with $s$
- `BSE/absp_io.f90:67` — `write_eigenvalues`
- `BSE/vmtxel.f90:380` — `compute_ik_vmtxel`, choice of operator / polarization
- `Common/mtxel_optical.f90` — `mtxel_m` (momentum), `mtxel_v` (velocity)
- `BSE/input_fi_q.f90:257` — definition of `xct%qshift`
- `BSE/absp.f90:65` — $\varepsilon_2$ prefactor
