# DFT Exact Exchange Singularity Subtraction

## Quantum Espresso

Source tree inspected: `/global/cfs/cdirs/m4480/emma/code/qe/qe-cpu-7.3.1`.

### Summary

Quantum Espresso implements the Coulomb-singularity treatment in the `pw.x` exact-exchange infrastructure. It is not part of the ordinary local or semi-local exchange-correlation evaluation: LDA, GGA, and meta-GGA pieces are handled through the usual XC/potential path, while the singularity machinery is entered when `xclib_dft_is('hybrid')` is true and the Fock operator is active.

For a hybrid or Hartree-Fock calculation, setup calls `setup_exx()`, which initializes the exchange q-grid and calls `exx_div_check()`. During the hybrid SCF loop, `exxinit()` computes the scalar `exxdiv = exx_divergence()`, and the Fock application, exchange-energy, stress, and ACE-construction routines use the regularized reciprocal-space kernel from `g2_convolution()`/`g2_convolution_all()`.

### Main Source Files

The files most directly involved are:

- `PW/src/exx_base.f90`: central implementation. It defines `exxdiv_treatment`, `x_gamma_extrapolation`, `exxdiv`, `nq1/nq2/nq3`, `coulomb_fac`, `exx_div_check()`, `g2_convolution()`, `g2_convolution_all()`, and `exx_divergence()`.
- `PW/src/exx.f90`: main Fock/EXX implementation. `exxinit()` obtains screening parameters and computes `exxdiv`; `vexx_gamma()`, `vexx_k()`, `exxenergy2_gamma()`, `exxenergy2_k()`, and `aceinit()` consume the regularized `coulomb_fac`.
- `PW/src/setup.f90` and `PW/src/run_pwscf.f90`: call `setup_exx()`/`reset_exx()` for hybrid runs and after G-vector/cell changes.
- `PW/src/electrons.f90`: hybrid SCF loop. It performs an initial non-hybrid SCF, activates exact exchange, then iterates with a fixed/updateable Fock operator; its `dexx` diagnostic explicitly warns that a negative value can indicate failure of the EXX divergence treatment.
- `PW/src/wfcinit.f90` and `PW/src/h_psi.f90`: non-SCF/bands-side usage. `aceinit0()` reads an ACE potential from restart files for non-SCF calculations, and `h_psi()` applies either `vexxace_*` or the full `vexx` operator when `exx_is_active()`.
- `Modules/coulomb_vcut.f90`: Wigner-Seitz and spherical Coulomb-cutoff alternatives used when `exxdiv_treatment='vcut_ws'` or `'vcut_spherical'`.
- `PW/Doc/INPUT_PW.txt` and `PW/examples/EXX_example/README`: input documentation and the clearest prose description of the q -> 0 treatment.

There is not a separate singularity-subtraction implementation for SCF versus NSCF/bands. SCF hybrid calculations actively build the regularized Fock operator and, by default, its ACE representation. Non-SCF hybrid calculations normally read and apply the ACE operator from the preceding SCF restart, so the singularity treatment is inherited from the SCF calculation. The full `vexx` path still calls the same `g2_convolution()` machinery if ACE is disabled, but the EXX example README notes a practical limitation for band structures at new k-points: the code would need all bands at k+q points that were not generated in the SCF run.

### When It Is Used

The input variable `exxdiv_treatment` is documented as "Specific for EXX" and is checked only inside the hybrid setup path:

```text
setup()
  if xclib_dft_is('hybrid'):
    setup_exx()
      exx_grid_init()
      exx_mp_init()
      exx_div_check()
```

The local/semi-local part of a hybrid functional is not singularity-subtracted; only the nonlocal exact-exchange Coulomb kernel is. Pure local/semi-local functionals do not enter this path. Within exact exchange:

- Bare HF/PBE0-like exchange uses the 4 pi e^2 / |q+G|^2 kernel.
- HSE-like screened exchange uses `erfc_scrlen = screening_parameter` and the kernel factor `1 - exp(-|q+G|^2 / (4 erfc_scrlen^2))`.
- GAU-PBE uses a Gaussian kernel that is finite at q -> 0; the README and input docs say to use `exxdiv_treatment='none'` and `x_gamma_extrapolation=.false.`.

### Implemented Scheme

The default scheme is `exxdiv_treatment='gygi-baldereschi'`, with comments and documentation pointing to Gygi and Baldereschi, Phys. Rev. B 34, 4405 (1986), and Appendix A.5 of the Quantum Espresso paper. The implementation follows the standard add-and-subtract idea: replace the divergent Brillouin-zone integrand by a regular residue plus an analytically/numerically evaluated compensating term.

For each Fock pair, `g2_convolution()` forms

```text
Q = q + G = xk - xkq + G
qq = |Q|^2
```

in Cartesian units (`tpiba2` is included in the code). Away from the singular point, the reciprocal-space kernel is

```text
v(Q) = 4 pi e^2 / |Q|^2
```

or, for HSE-like screened exchange,

```text
v_erfc(Q) = 4 pi e^2 / |Q|^2 * [1 - exp(-|Q|^2 / (4 omega^2))]
```

where `omega` is `erfc_scrlen`. The code stores this in `fac(ig)` and multiplies by `grid_factor_track(ig)`. If `|Q|^2 <= eps_qdiv`, the special value is

```text
fac(Q = 0) = -exxdiv
```

with additional finite screened/Yukawa terms added when `x_gamma_extrapolation` is disabled. Thus `exxdiv` is the finite replacement generated by the singularity-subtraction calculation.

`exx_divergence()` computes `exxdiv` using an Ewald-damped reciprocal-space sum over the EXX q-grid and G-vectors, with

```text
alpha = 10 / gcutw
S = sum_q sum_G' exp(-alpha |q+G|^2) / |q+G|^2
```

where the prime omits points on the doubled grid when `x_gamma_extrapolation` is enabled and omits the singular `q+G=0` point. In the bare-Coulomb case, QE then applies the normalization and analytic subtraction:

```text
div = (4 pi e^2 / tpiba^2 / nqs) * S
aa  = 8/(4 pi) * integral_0^{5/sqrt(alpha)} dq [ ... ] + 1/sqrt(alpha*pi)
div = div - e^2 * Omega * aa
exxdiv = nqs * div
```

The code evaluates the one-dimensional integral numerically with `nqq = 100000` points. For HSE-like screened exchange, the sum and integral include the screened factors from `erfc_scrlen`; for `exxdiv_treatment='none'`, `exx_divergence()` returns zero.

By default, `x_gamma_extrapolation=.true.` also implements QE's double-grid extrapolation of the non-analytic `Q -> 0` limit. In `g2_convolution()`, points that lie on the doubled q-grid get `grid_factor_track=0`, while the remaining points get `grid_factor=8/7`. The EXX README describes this as internally extracting the missing `Q=0` contribution from a calculation on the requested q-grid and one twice as coarse in each direction. This is incompatible with the `vcut_ws` and `vcut_spherical` options.

For strongly anisotropic cells, QE provides Coulomb-cutoff alternatives:

```text
exxdiv_treatment = 'vcut_ws'
exxdiv_treatment = 'vcut_spherical'
```

These bypass the Gygi-Baldereschi branch in `g2_convolution()` and call `vcut_get()` or `vcut_spheric_get()` from `Modules/coulomb_vcut.f90`. The Wigner-Seitz cutoff builds a superperiodic lattice with columns scaled by `nq1`, `nq2`, and `nq3`; `ecutvcut` controls the reciprocal-space cutoff for the correction.

### References Found in Comments or Documentation

- Gygi and Baldereschi, Phys. Rev. B 34, 4405 (1986): cited in `PW/src/exx_base.f90` and `PW/examples/EXX_example/README` for the default singularity treatment.
- Quantum Espresso paper, Appendix A.5: referenced in the EXX example README for the q -> 0 treatment.
- Chawla and Voth, J. Chem. Phys. 108, 4697 (1998); Sorouri, Foulkes, and Hine, J. Chem. Phys. 124, 064105 (2006); Spencer and Alavi, Phys. Rev. B 77, 193110 (2008): cited in the EXX example README as algorithmic references.

## INQ

Source tree inspected: `/global/homes/e/emsi/code/inq/inq-latest_main`.

### Summary

INQ implements singularity subtraction for the non-local exact-exchange/Fock operator used by Hartree-Fock and hybrid functionals. The implementation is not part of the ordinary local or semi-local exchange-correlation path. Local, semi-local, and meta-GGA functionals are evaluated through libxc in `hamiltonian::xc_term`/`hamiltonian::xc_functional`, and the ordinary Hartree potential uses the Poisson solver with its default zero `G = 0` term.

The singularity-subtraction implementation is based on the method of Carrier et al., Phys. Rev. B 75, 205126 (2007), explicitly referenced in `src/ionic/singularity_correction.hpp`.

### Main Source Files

The files most directly involved are:

- `src/ionic/singularity_correction.hpp`: implements `ionic::singularity_correction`, including the auxiliary function, the precomputed `f_k` values, `f_0`, and the per-k-point correction returned by `operator()(ik)`.
- `src/ionic/brillouin.hpp`: stores a lazy, cached `std::optional<ionic::singularity_correction>` on the Brillouin-zone object. `brillouin::sing_correction(cell, comm)` constructs the correction when first requested.
- `src/hamiltonian/exchange_operator.hpp`: owns the optional singularity correction in `hamiltonian::exchange_operator`. During exact exchange, it passes the correction as the Poisson `zeroterm` for each exchange pair.
- `src/solvers/poisson.hpp`: provides the Fourier-space Poisson kernels. In the 3D periodic kernel, the `G = 0` value is supplied externally through `zeroterm`; otherwise the kernel is `-1/G^2`, later multiplied by `-4 pi`.
- `src/hamiltonian/self_consistency.hpp`: decides the exact-exchange coefficient and builds the exchange operator in `update_exchange()`.
- `src/hamiltonian/ks_hamiltonian.hpp`: applies the adaptively compressed exchange (ACE) representation produced by the exact-exchange operator.
- `src/ground_state/calculator.hpp`: ground-state SCF driver; calls `self_consistency::update_exchange()` initially and periodically during SCF.
- `src/real_time/propagate.hpp` and `src/real_time/crank_nicolson.hpp`: real-time driver; exact exchange is supported through the Crank-Nicolson/parallel-transport path, which repeatedly calls `update_exchange()`. ETRS is rejected when exact exchange is present.
- `src/options/theory.hpp` and `src/hamiltonian/xc_functional.hpp`: define Hartree-Fock and hybrid-functional exact-exchange coefficients.

### Ground State vs Real Time

There is not a separate singularity-subtraction implementation for ground-state and real-time calculations. Both paths use the same exact-exchange stack:

```text
self_consistency::update_exchange()
  -> exchange_operator
  -> ionic::singularity_correction from brillouin::sing_correction()
  -> poisson::in_place(..., gshift, zeroterm)
  -> ACE orbitals stored on ks_hamiltonian
```

In ground-state calculations, `src/ground_state/calculator.hpp` calls `update_exchange()` before the SCF loop and then updates exact exchange every few SCF iterations. The convergence printout includes `dexx` only when `self_consistency::has_exact_exchange()` is true.

In real-time calculations, `src/real_time/propagate.hpp` calls `update_exchange()` during initialization. For exact-exchange functionals, ETRS is disallowed, and `src/real_time/crank_nicolson.hpp` updates exact exchange inside the Crank-Nicolson self-consistency loop and again at the end of the step.

### When It Is Used

Singularity subtraction is used only when the exact-exchange operator is enabled:

```cpp
bool enabled() const {
	return fabs(exchange_coefficient_) > 1.0e-14;
}
```

`exchange_operator` constructs the singularity correction only under that condition:

```cpp
if(enabled()) {
	sing_.emplace(electrons.brillouin_zone().sing_correction(electrons.states_basis().cell(), electrons.full_comm()));
}
```

The exact-exchange coefficient comes from two cases:

- Pure Hartree-Fock: `options::theory::hartree_fock(coeff)` sets `exchange_ = XC_HARTREE_FOCK` and uses `hf_coefficient_`.
- Hybrid libxc functionals: `xc_functional::exx_coefficient()` returns `xc_hyb_exx_coef(...)` for libxc hybrid functionals such as PBE0/PBEH or B3LYP.

For local/semi-local functionals, such as LDA, GGA, and non-hybrid MGGA, `exx_coefficient()` is zero and the exchange operator is disabled. In that case, no `ionic::singularity_correction` is constructed and no nonzero Poisson `zeroterm` is supplied by the exchange path.

The ordinary Hartree potential is different: in `self_consistency::update_hamiltonian()`, it is computed as:

```cpp
auto vhartree = solvers::poisson::solve(total_density);
```

That path uses the default Poisson zero term. It does not use the exact-exchange singularity correction.

### Implemented Scheme

The code comment in `src/ionic/singularity_correction.hpp` says:

```cpp
// Implements the method for the singularity correction of Carrier et al. PRB 75 205126 (2007)
// https://doi.org/10.1103/PhysRevB.75.205126
```

The correction follows the Carrier et al. auxiliary-function approach. In the code, reciprocal-space cell metrics are first formed as

```text
d1_j = 4 b_j . b_j
d2_j = 2 b_j . b_{j+1}
```

where `b_j` are reciprocal lattice vectors and `j + 1` is cyclic. The auxiliary function is labeled in the source as "the function defined in Eq. 16":

```text
F(q) = 4 pi^2 / [ 1/2 v1(q) + v2(q) ]
```

with

```text
v1(q) = d1_0 [1 - cos(p_0)] + d1_1 [1 - cos(p_1)] + d1_2 [1 - cos(p_2)]
v2(q) = d2_0 sin(p_0) sin(p_1)
      + d2_1 sin(p_1) sin(p_2)
      + d2_1 sin(p_2) sin(p_0)
```

Here `p_j = a_j . q_cart`, as computed by `cell_projection()`, with `a_j` the real-space lattice vectors. Note: the last term above reflects the source exactly; the code uses `dp2[1]` for both the second and third cross terms.

The code precomputes two quantities:

1. `f_k`, evaluated for the actual k-point grid:

```text
f_k(i) = sum_{j != i} [ 4 pi / Omega * w_j * F(k_i - k_j) ]
```

where `Omega` is the cell volume and `w_j` is the k-point weight.

1. `f_0`, evaluated by a numerical integration/refinement procedure over a dense mesh:

```text
f_0 ~= [8 pi / (2 pi)^3] Delta k^3 sum_q sum_{s=0}^{6} l_s^3 F(l_s q)
       + analytic tail
```

with

```text
l_s = 3^{-s}
Delta k^3 = [2 pi / (2 n_k + 1)]^3 / Omega
n_k = 60
```

The final per-k correction returned by `singularity_correction::operator()(ik)` is:

```text
C_i = -N_k Omega [ f_k(i) - f_0 ]
```

where `N_k` is the number of k points. This `C_i` is the value passed to the Poisson solver as the zero-frequency replacement for the exchange-pair solve.

### Where the Correction Enters the Fock Operator

In `exchange_operator::block_exchange()`, for each occupied/reference orbital `j`, INQ forms a pair density

```text
rho_ij(r) = phi_j(r)^* phi_i(r)
```

then solves Poisson with a k-point shift:

```cpp
solvers::poisson::in_place(rhoij, -phi.kpoint() + kpt[jj], (*sing_)(idx[jj]));
```

The third argument is the singularity-subtraction correction, installed by the Poisson solver as the special `G = 0` contribution of the shifted Coulomb solve.

The exchange action then accumulates

```text
(V_x phi_i)(r) += -1/2 * alpha_x * f_j * phi_j(r) * v_ij(r)
```

where `alpha_x` is the exact-exchange fraction and `f_j` is the occupation. The resulting Fock action is compressed into ACE orbitals by `exchange_operator::ace()` and later applied in `ks_hamiltonian::ace()`.

### References Found in Comments

- Carrier et al., Phys. Rev. B 75, 205126 (2007): explicitly cited in `src/ionic/singularity_correction.hpp` for the singularity correction.
- Fock (1930) and Hartree/Hartree (1935): cited in `src/hamiltonian/xc_functional.hpp` for Hartree-Fock references.
- Ihm, Zunger, and Cohen, J. Phys. C 12, 4409 (1979): cited in `src/hamiltonian/self_consistency.hpp` for a separate pseudopotential `G = 0` energy correction, not the exact-exchange singularity subtraction.
- Rozzi et al., Phys. Rev. B 73, 205119 (2006): cited in `src/solvers/poisson.hpp` for the finite/0D Coulomb kernel, not the exact-exchange singularity subtraction.

I did not find Baldereschi, Gygi, or Spencer-style singularity-subtraction references in the INQ source files involved in the exact-exchange path.
