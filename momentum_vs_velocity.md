# Why `use_momentum` divides by an energy difference and `use_velocity` does not

Units: Rydberg atomic units, so $\hbar=1$, $m=1/2$, $e^2=2$, $\hat{\mathbf p}=-i\nabla$.
Bloch states are $\psi_{n\mathbf k}=e^{i\mathbf k\cdot\mathbf r}u_{n\mathbf k}$, and the Bloch
Hamiltonian is $H_{\mathbf k}=e^{-i\mathbf k\cdot\mathbf r}He^{i\mathbf k\cdot\mathbf r}$, with
$H_{\mathbf k}u_{n\mathbf k}=\varepsilon_{n\mathbf k}u_{n\mathbf k}$.

## 1. What the code computes

The two operators are computed in `Common/mtxel_optical.f90`:

- **`mtxel_v`** (velocity):
  `s0 = sum_G conjg(cg_c(G)) * cg_vq(G) / qshift`. Here `cg_vq` comes from
  `WFNq_fi` and the two G-vectors are matched through `isorti`. So
  $$s^{\rm vel}_{cv\mathbf k}=\frac{\langle u_{c\mathbf k}|u_{v,\mathbf k+\mathbf q}\rangle}{|\mathbf q|}.$$
  No energies appear.
- **`mtxel_m`** (momentum):
  `s0 = 2*sum_G conjg(cg_c(G))*cg_v(G)*(e.G)`, then `s0 = s0/de`, where
  `de = eqp%eclda(ic,ik) - eqp%evlda(iv,ik)`. Both states come from `WFN_fi`.
  `eclda`/`evlda` are the mean-field eigenvalues read from `WFN_fi`
  (`BSE/input_fi.f90:398`: `eqp%eclda(...) = kp%elda(...)`). So
  $$s^{\rm mom}_{cv\mathbf k}=\frac{\langle u_{c\mathbf k}|\,\hat{\mathbf e}\cdot 2(\hat{\mathbf p}+\mathbf k)\,|u_{v\mathbf k}\rangle}{\varepsilon_{c\mathbf k}-\varepsilon_{v\mathbf k}}.$$
  The code omits $\mathbf k$, which is harmless: $\langle u_c|u_v\rangle=0$ for $c\neq v$.

## 2. Which quantity the absorption formula needs

The comment in `BSE/absp0.f90:38-44` derives the prefactor from RMP 74, 601 (2002), eq. (4.12):
$$\varepsilon(\mathbf q,\omega)=1-v(\mathbf q)\sum_S\Big|\sum_{cv\mathbf k}A^S_{cv\mathbf k}\langle c\mathbf k|e^{-i\mathbf q\cdot\mathbf r}|v\mathbf k+\mathbf q\rangle\Big|^2\frac{1}{\omega-\Omega_S+i\eta}+\dots,\qquad v(\mathbf q)=\frac{4\pi e^2}{q^2}.$$
Moving the $1/q^2$ of $v(\mathbf q)$ into the matrix elements gives the $q\to0$ limit of
$\langle c\mathbf k|e^{-i\mathbf q\cdot\mathbf r}|v\mathbf k+\mathbf q\rangle/q
=\langle u_{c\mathbf k}|u_{v\mathbf k+\mathbf q}\rangle/q$.

That is exactly what `mtxel_v` evaluates at finite $q$. The velocity option therefore computes the required quantity directly. The momentum option computes something else and has to convert it.

## 3. Why the conversion requires an energy denominator

Taylor-expand in $\mathbf q$ and use $\langle u_{c\mathbf k}|u_{v\mathbf k}\rangle=0$:
$$\frac{\langle u_{c\mathbf k}|u_{v\mathbf k+\mathbf q}\rangle}{|\mathbf q|}=\hat{\mathbf q}\cdot\langle u_{c\mathbf k}|\nabla_{\mathbf k}u_{v\mathbf k}\rangle+O(q). \tag{A}$$
This is a derivative of the wavefunctions with respect to $\mathbf k$. It contains no energies and makes no reference to $H$. `use_velocity` evaluates it by finite differences, and that is why it needs a second file, `WFNq_fi`.

`use_momentum` reads only `WFN_fi`, so it cannot differentiate $u$. Instead it uses an identity that holds for **any** Hermitian $H_{\mathbf k}$. Differentiate $H_{\mathbf k}u_{v\mathbf k}=\varepsilon_{v\mathbf k}u_{v\mathbf k}$ and project onto $u_{c\mathbf k}$ ($c\neq v$):
$$\langle u_c|\nabla_{\mathbf k}H_{\mathbf k}|u_v\rangle+\varepsilon_c\langle u_c|\nabla_{\mathbf k}u_v\rangle=\varepsilon_v\langle u_c|\nabla_{\mathbf k}u_v\rangle
\;\Longrightarrow\;
\langle u_c|\nabla_{\mathbf k}u_v\rangle=\frac{\langle u_c|\nabla_{\mathbf k}H_{\mathbf k}|u_v\rangle}{\varepsilon_v-\varepsilon_c}. \tag{B}$$

Equation (B) turns the wavefunction derivative into a matrix element of an operator, $\nabla_{\mathbf k}H_{\mathbf k}=i\,e^{-i\mathbf k\cdot\mathbf r}[H,\mathbf r]e^{i\mathbf k\cdot\mathbf r}$ (the velocity operator). The price is the energy denominator. That is the whole answer:

- **Velocity:** the code differentiates $u$ numerically, so no denominator is needed.
- **Momentum:** the code replaces the derivative with (B). This introduces $1/(\varepsilon_c-\varepsilon_v)$ and also requires an expression for $\nabla_{\mathbf k}H_{\mathbf k}$.

Two constraints follow from (B):

- The $\varepsilon$ in (B) must be eigenvalues of the **same** $H_{\mathbf k}$ whose eigenvectors are $u$, i.e. the mean-field eigenvalues, not the GW ones. The code uses `eclda`/`evlda`, which satisfies this.
- (B) is singular when $\varepsilon_c=\varepsilon_v$. The code sets $s^{\rm mom}=0$ when $|\varepsilon_c-\varepsilon_v|<$ `TOL_Degeneracy`.

`mtxel_m` approximates $\nabla_{\mathbf k}H_{\mathbf k}$ as $2(\hat{\mathbf p}+\mathbf k)$, which is the kinetic term only. Comparing with (A) and (B):
$$s^{\rm mom}=-\hat{\mathbf e}\cdot\langle u_c|\nabla_{\mathbf k}u_v\rangle\quad\text{iff}\quad \langle u_c|\nabla_{\mathbf k}V_{\mathbf k}|u_v\rangle=0,$$
where $V$ is the full non-kinetic part of $H$. This holds exactly when $V$ commutes with $\mathbf r$, i.e. when $V$ is a multiplicative potential. The overall sign is opposite to $s^{\rm vel}$, which has no effect on $|d_S|^2$.

## 4. Case 1: semilocal functional (PBE) with norm-conserving pseudopotentials

$H=\hat p^2+V_{\rm loc}^{\rm PP}(\mathbf r)+V_H(\mathbf r)+v_{xc}^{\rm PBE}(\mathbf r)+\hat V_{\rm NL}^{\rm PP}$.

- $V_{\rm loc}$, $V_H$ and $v_{xc}^{\rm PBE}=\delta E_{xc}/\delta n(\mathbf r)$ are multiplicative functions of $\mathbf r$, so they commute with $\mathbf r$.
- $\hat V_{\rm NL}^{\rm PP}$ (Kleinman–Bylander projectors) is not multiplicative, so $[\hat V_{\rm NL},\mathbf r]\neq0$.

**Momentum:** the denominator is the exact PBE eigenvalue difference, consistent with (B). The numerator misses $\langle u_c|\nabla_{\mathbf k}V^{\rm NL}_{\mathbf k}|u_v\rangle$.

- Code evidence: `mtxel_m` uses only plane-wave coefficients and $\mathbf G$, with no projector data.
- Documentation evidence: `documentation/input_files/absorption.inp:142-146` says: *"When you use the momentum operator, you are throwing away the contribution from the non-local part of the pseudopotential."*
- I have not quantified the size of this term for your system.

**Velocity:** (A) involves no $H$, so the nonlocal part is included automatically, provided `WFN_fi` and `WFNq_fi` are eigenstates of the same $H$ (see the next paragraph). `absorption.inp:136-139` says the same: velocity *"includes effects such as the non-local part of the pseudopotential."* The error is $O(q)$ from truncating (A).

**Condition for velocity to be correct:** if `WFNq_fi` is an eigenstate of $H+\delta H$ instead of $H$, first-order perturbation theory adds $\langle u_c|\delta H|u_v\rangle/[(\varepsilon_v-\varepsilon_c)|\mathbf q|]$ to $s^{\rm vel}$. Because it is divided by $|\mathbf q|$, it is not suppressed as $q\to0$. With PBE, generating both files by non-self-consistent runs from one SCF density fixes $V_{\rm KS}$, so $\delta H=0$.

## 5. Case 2: hybrid functional (PBE0)

$H_{\rm PBE0}=H_{\rm PBE}+\tfrac14\big(\hat V_x^{\rm HF}-v_x^{\rm PBE}\big)$, with
$(\hat V_x^{\rm HF}\psi)(\mathbf r)=-\sum_{j\in\rm occ}\phi_j(\mathbf r)\!\int\!\phi_j^*(\mathbf r')v(\mathbf r-\mathbf r')\psi(\mathbf r')\,d\mathbf r'$.

This is an integral operator, not a multiplicative one, so $[\hat V_x^{\rm HF},\mathbf r]\neq0$.

**Momentum:** the denominator is now the PBE0 eigenvalue difference. It is still consistent with (B), because $u$ are eigenvectors of $H_{\rm PBE0}$. The numerator still uses only $2(\hat{\mathbf p}+\mathbf k)$, so it now misses two terms:
$$\langle u_c|\nabla_{\mathbf k}V^{\rm NL}_{\mathbf k}|u_v\rangle+\tfrac14\langle u_c|\nabla_{\mathbf k}V^{x,\rm HF}_{\mathbf k}|u_v\rangle .$$
The second term is new relative to PBE.

**Velocity:** (A) still makes no reference to $H$, so both nonlocal terms are included, subject to the same condition as in section 4 ($\delta H=0$). Here that condition is harder to meet:

- $\hat V_x^{\rm HF}$ depends on the occupied orbitals $\phi_j$ on the k/q mesh used to build it.
- A separate SCF on a shifted mesh would in general build a different Fock operator, so $\delta H\neq0$. That error is not suppressed as $q\to0$ (section 4).
- I verified that your `WFN_fi` (`01-scf/scf.in`) is `calculation='scf'`, `input_dft='PBE0'`, 8×8×8 mesh with 1 1 1 shift.
- `WFNq_fi` links to `../02-scf_q1e-5/WFN_in.h5`, but that directory no longer exists. So I **could not verify** how it was generated, or whether its Fock operator matches the one in `WFN_fi`.

**What your data shows (PBE0):** the f-sum rules printed in the absorption outputs are:

| run | non-interacting sum rule | BSE sum rule |
|---|---|---|
| velocity | 1.035 | 0.886 |
| momentum | 0.768 (b1), 0.768 (b2), 0.768 (b3) | 0.660 (b1), 0.660 (b2), 0.660 (b3) |

These numbers alone cannot split the difference into its three possible sources:

- missing $V_{\rm NL}$,
- missing $\hat V_x^{\rm HF}$,
- a possible $\delta H$ between `WFN_fi` and `WFNq_fi`.

To separate them, repeat the same comparison with PBE, where the Fock term is absent and $\delta H=0$ can be guaranteed.

## Summary

| | computes | energy denominator | $[V_{\rm NL}^{\rm PP},\mathbf r]$ | $[\hat V_x^{\rm HF},\mathbf r]$ (PBE0 only) | requirement |
|---|---|---|---|---|---|
| velocity | $\langle u_c\|\nabla_{\mathbf k}u_v\rangle$ by finite difference, eq. (A) | none needed | included | included | `WFN_fi`, `WFNq_fi` eigenstates of the same $H$; $O(q)$ error |
| momentum | $\langle u_c\|\nabla_{\mathbf k}H_{\mathbf k}\|u_v\rangle/(\varepsilon_v-\varepsilon_c)$, eq. (B), kinetic part only | required by (B); mean-field eigenvalues | omitted | omitted | $\varepsilon_c\neq\varepsilon_v$ |
