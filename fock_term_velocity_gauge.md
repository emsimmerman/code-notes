# Fock Term Under Length → Velocity Gauge Transformation in TDKS

## Conventions (check these against the code first)

Everything below assumes:

- **Atomic units**, $\hbar = m_e = |e| = 1$, and $c$ absorbed into $\mathbf{A}$. If the code carries explicit $1/c$ factors, the phase is $e^{i\mathbf{A}\cdot\mathbf{r}/c}$ throughout; the cancellation structure is unchanged.
- **Field/potential relation** $\mathbf{E}(t) = -\partial_t\mathbf{A}(t)$, with no scalar potential ($\phi = 0$) in either gauge.
- **Charge sign** absorbed such that the length-gauge coupling is $+\mathbf{E}\cdot\mathbf{r}$ and the velocity-gauge kinetic term is $(p+\mathbf{A})^2/2$. These two must be consistent with each other; flipping one without the other breaks the cancellation in Step 3.
- **Gauge phase** $U(t) = e^{+i\mathbf{A}(t)\cdot\mathbf{r}}$, mapping velocity-gauge states to length-gauge states, $\psi = U\psi'$.

If the code instead defines $U = e^{-i\mathbf{A}\cdot\mathbf{r}}$ (or equivalently $\psi' = U\psi$), then every phase below flips sign: the kernel prefactor becomes $e^{-i\mathbf{A}\cdot(\mathbf{r}-\mathbf{r}')}$ and the covariance relation reads $\hat F[\{U\psi'\}] = U^\dagger \hat F[\{\psi'\}] U$. The pairwise cancellation is identical either way — the thing to verify is that a *single* consistent choice is used in the propagator, the kinetic term, and the exchange build.

## Notation

Fix the gauge transformation. With $U(t) = e^{i\mathbf{A}(t)\cdot\mathbf{r}}$, a multiplicative (local) unitary:

$$\psi_j(\mathbf{r},t) = U\,\psi_j'(\mathbf{r},t) = e^{i\mathbf{A}\cdot\mathbf{r}}\psi_j'(\mathbf{r},t)$$

unprimed = length gauge, primed = velocity gauge. The same $U$ acts on every orbital $j$; this is essential below.

**The functional.** Exchange is a map from a *set* of orbitals to an operator. For any set $\{\chi_k\}$:

$$\big(\hat F[\{\chi\}]\,\varphi\big)(\mathbf{r}) \;=\; -\sum_k \chi_k(\mathbf{r})\int d\mathbf{r}'\, v(|\mathbf{r}-\mathbf{r}'|)\,\chi_k^*(\mathbf{r}')\,\varphi(\mathbf{r}')$$

**The two kernels.** Same functional, different arguments:

$$K_L(\mathbf{r},\mathbf{r}') \equiv -\sum_k \psi_k(\mathbf{r})\,\psi_k^*(\mathbf{r}')\,v(|\mathbf{r}-\mathbf{r}'|) \qquad \text{(length-gauge orbitals in)}$$

$$K_V(\mathbf{r},\mathbf{r}') \equiv -\sum_k \psi_k'(\mathbf{r})\,\psi_k'^*(\mathbf{r}')\,v(|\mathbf{r}-\mathbf{r}'|) \qquad \text{(velocity-gauge orbitals in)}$$

So $K_L = $ kernel of $\hat F[\{\psi\}]$ and $K_V = $ kernel of $\hat F[\{\psi'\}]$. Note both are defined with *no* explicit $\mathbf{A}$ in them; the only $\mathbf{A}$-dependence enters through which orbitals you feed in. Substituting $\psi_k = U\psi_k'$ into $K_L$ gives immediately

$$K_L(\mathbf{r},\mathbf{r}') = e^{i\mathbf{A}\cdot(\mathbf{r}-\mathbf{r}')}\,K_V(\mathbf{r},\mathbf{r}') \quad\Longleftrightarrow\quad \hat F[\{U\psi'\}] = U\,\hat F[\{\psi'\}]\,U^\dagger$$

which is the covariance relation. $K_V \neq K_L$ as functions — that residual $e^{i\mathbf{A}\cdot(\mathbf{r}-\mathbf{r}')}$ is real and $\mathbf{A}$-dependent. Nothing yet about the equation of motion.

## The Fock term in TDKS, step by step

**Step 0 — length gauge, nothing done.**

$$i\partial_t\psi_j = \Big[\tfrac{p^2}{2} + V_{\rm ext} + \mathbf{E}\cdot\mathbf{r} + V_H\Big]\psi_j + \hat F[\{\psi\}]\psi_j$$

with the Fock term written out:

$$\big(\hat F[\{\psi\}]\psi_j\big)(\mathbf{r}) = -\sum_k \psi_k(\mathbf{r})\int d\mathbf{r}'\,v(|\mathbf{r}-\mathbf{r}'|)\,\psi_k^*(\mathbf{r}')\,\psi_j(\mathbf{r}')$$

**Step 1 — substitute $\psi = U\psi'$. Three places, no simplification.**

The orbital appears three times in this term: the outer factor at $\mathbf{r}$, the conjugate at $\mathbf{r}'$, and the acted-on orbital at $\mathbf{r}'$.

$$= -\sum_k \underbrace{e^{i\mathbf{A}\cdot\mathbf{r}}\psi_k'(\mathbf{r})}_{\text{outer}}\int d\mathbf{r}'\,v(|\mathbf{r}-\mathbf{r}'|)\,\underbrace{e^{-i\mathbf{A}\cdot\mathbf{r}'}\psi_k'^*(\mathbf{r}')}_{\text{kernel, }\mathbf{r}'}\;\underbrace{e^{i\mathbf{A}\cdot\mathbf{r}'}\psi_j'(\mathbf{r}')}_{\text{acted-on}}$$

Three phases. This is the "before any manipulation" object, and this is where the claim lives: the term carries $U$-dependence from the orbitals inside the functional *and* will shortly acquire one more from the conjugation.

**Step 2 — the pair at $\mathbf{r}'$ cancels in place.**

$e^{-i\mathbf{A}\cdot\mathbf{r}'}e^{+i\mathbf{A}\cdot\mathbf{r}'} = 1$. This works because both phases are evaluated at the *same* point, and because $U$ is state-independent so the $k$-sum doesn't care:

$$= -e^{i\mathbf{A}\cdot\mathbf{r}}\sum_k \psi_k'(\mathbf{r})\int d\mathbf{r}'\,v(|\mathbf{r}-\mathbf{r}'|)\,\psi_k'^*(\mathbf{r}')\,\psi_j'(\mathbf{r}') \;=\; e^{i\mathbf{A}\cdot\mathbf{r}}\big(\hat F[\{\psi'\}]\psi_j'\big)(\mathbf{r})$$

$$\text{i.e.}\qquad \hat F[\{U\psi'\}]\,U\psi_j' \;=\; U\,\hat F[\{\psi'\}]\,\psi_j'$$

One phase survives, sitting out front at $\mathbf{r}$.

**Step 3 — the equation of motion supplies the conjugation.**

$i\partial_t(U\psi_j') = U(i\partial_t\psi_j') + (i\partial_t U)\psi_j'$. Move the anomalous piece to the right-hand side and multiply the whole equation from the left by $U^\dagger$. Every term gets conjugated, including exchange:

$$U^\dagger\,\hat F[\{U\psi'\}]\,U\psi_j' \;=\; U^\dagger\,U\,\hat F[\{\psi'\}]\,\psi_j' \;=\; \hat F[\{\psi'\}]\,\psi_j'$$

**Result.**

$$i\partial_t\psi_j' = \Big[\tfrac{(p+\mathbf{A})^2}{2} + V_{\rm ext} + V_H\Big]\psi_j' + \hat F[\{\psi'\}]\psi_j'$$

The dipole term $\mathbf{E}\cdot\mathbf{r}$ was cancelled by $-iU^\dagger\dot U$, and the exchange term is the same functional with primed orbitals in it — kernel $K_V$, no explicit $\mathbf{A}$.

## Bookkeeping summary

Three phases entered, and they cancel pairwise across the internal/external divide, not internally:

| phase | origin | cancelled by |
|---|---|---|
| $e^{-i\mathbf{A}\cdot\mathbf{r}'}$ | $\psi_k^* \to \psi_k'^*$ inside kernel | the $U$ acting on $\psi_j'$ (Step 2) |
| $e^{+i\mathbf{A}\cdot\mathbf{r}}$ | $\psi_k \to \psi_k'$ outside integral | $U^\dagger$ from conjugation (Step 3) |
| $e^{+i\mathbf{A}\cdot\mathbf{r}'}$ | the acted-on orbital | pairs with row 1 |

So the claim is confirmed and now sharper: the double $U$-dependence is real, and *both* are needed. Keep only the conjugation (freeze the density matrix at its length-gauge or ground-state value) and you're left with $e^{i\mathbf{A}\cdot(\mathbf{r}-\mathbf{r}')}$ multiplying the kernel — gauge-dependent. Keep only the internal substitution (forget to conjugate) and you're left with a stray $e^{i\mathbf{A}\cdot\mathbf{r}}$ out front, which is not even a valid operator on the primed states.
