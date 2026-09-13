# k-space diagonalisation, Zak phase and cavity dressing of the SSH chain

Step-by-step computational recipe behind `k_space_diag.ipynb`. Every step states its **input**,
the **formulas**, the **output**, and a **check** that pins it down numerically. Written so the
whole thing can be re-implemented from the maths alone.

System: the `PBC=True` bosonic `PHOTONICChain` of
`packages/qcg-quantum-design/src/ssh_chain_sq.py` — an SSH chain of $N$ unit cells (two transmons
each) with capacitive and cavity-mediated Josephson coupling.

---

## Part 0 — Model, indices and conventions

### 0.1 Indices

* $N$ unit cells, $n = 2N$ transmon modes.
* Site index $l = 2j+\sigma$, with $j = 0\ldots N-1$ the cell and $\sigma\in\{0,1\}=\{A,B\}$ the
  sublattice.
* Nambu vector $\Phi = (b_0,\dots,b_{n-1},\,b_0^\dagger,\dots,b_{n-1}^\dagger)^T$, length $2n$.

### 0.2 Coupling constants

$$g_{c1}=\tfrac14\Omega_C(1+\chi_C),\qquad g_{c2}=\tfrac14\Omega_C(1-\chi_C),$$
$$g_{J1}=\Omega_J(1+\chi_J),\qquad g_{J2}=-\Omega_J(1-\chi_J),\qquad g_J \equiv g_{J1}+g_{J2}.$$

Note $g_{c1}+g_{c2}=\tfrac12\Omega_C$ always.

### 0.3 The two normalisation conventions (this trips everything up — fix it once)

**Honest / Nambu convention.** $H=\tfrac12\Phi^\dagger\mathcal H\Phi+\text{const}$. Here $A$ has
diagonal $\tfrac12$ and the bare transmon frequency reads $\tfrac12$.

**Physical convention.** Multiply every block by $2$ (this is `_nambu_prefactor = 2` in
`ssh_chain_sq`). The bare transmon then sits at frequency $1$. All *reported* frequencies, all
Hamiltonians fed to the cavity problem, and the analytic block formulas of §1.5 use this one.

Rule of thumb: **build** in the honest convention (so the $\tfrac12\Phi^\dagger\mathcal H\Phi$
algebra and $\tau_3$ diagonalisation are standard), **report** in the physical one.

---

## Part I — Momentum space and block diagonalisation by permutation

### Step 1 — Real-space Nambu blocks

**Input:** $N,\Omega_C,\chi_C,\Omega_J,\chi_J$.

Build the SSH hopping matrix $T$ ($n\times n$, real symmetric):

$$T_{l,l+1}=T_{l+1,l}=\begin{cases} g_{c1} & l \text{ even (intra-cell)}\\ g_{c2} & l \text{ odd (inter-cell)}\end{cases}
\quad l=0\ldots n-2, \qquad T_{0,n-1}=T_{n-1,0}=g_{c2}\ \ \text{(PBC wrap)}.$$

Then

$$A=\tfrac12\big(I+T\big),\qquad B=-\tfrac12 T .$$

Identically for Josephson, with $T_J$ built the same way from $g_{J1},g_{J2}$:

$$A_J=\tfrac12\big(g_J I-T_J\big),\qquad B_J=-\tfrac12 T_J-\tfrac12 g_J I .$$

Assemble the $2n\times2n$ Nambu matrices

$$\mathcal H_C=\begin{pmatrix}A&B\\B&A\end{pmatrix},\qquad
  \mathcal H_J=\begin{pmatrix}A_J&B_J\\B_J&A_J\end{pmatrix}.$$

**Output:** $\mathcal H_C,\mathcal H_J$, both $2n\times2n$ real symmetric.

**Check:** these are exactly what `PHOTONICChain.__init__` builds (`D_herm`, `H_J_herm`).

---

### Step 2 — The Fourier transform as a change of operator basis

**Input:** $N$.

Momenta on the FFT grid (this matters — it is the grid on which $-k$ lands on another grid point):

$$k_m = 2\pi\,\texttt{fftfreq}(N)_m \quad\longrightarrow\quad \texttt{fftshift}\ \text{onto } (-\pi,\pi].$$

The transform itself is

$$b_{j\sigma}=\frac{1}{\sqrt N}\sum_{m} e^{ik_m j}\,c_{m\sigma}.$$

As a matrix, define $W_{jm}=e^{ik_m j}/\sqrt N$ ($N\times N$) and

$$\mathcal F_{(2j+\sigma),\,(2m+\sigma')} = W_{jm}\,\delta_{\sigma\sigma'} \qquad (n\times n,\ \text{unitary}),$$

i.e. $\mathcal F$ is $W$ replicated on the two sublattice diagonal blocks. Then $b=\mathcal F c$
and $b^\dagger=\mathcal F^*c^\dagger$, so the Nambu lift is

$$\mathcal G=\begin{pmatrix}\mathcal F&0\\0&\mathcal F^*\end{pmatrix},\qquad
\boxed{\ \mathcal H^k=\mathcal G^\dagger\,\mathcal H\,\mathcal G\ }$$

**Output:** $\mathcal H_C^k,\mathcal H_J^k$, both $2n\times2n$.

**Check:** $\mathcal F^\dagger\mathcal F=I$; the spectrum of $\mathcal H^k$ equals that of
$\mathcal H$ to $\sim10^{-16}$.

---

### Step 3 — The permutation that block-diagonalises

**Input:** the momentum list $\{k_m\}$.

For each $m$, find $m'$ with $k_{m'}\equiv-k_m \pmod{2\pi}$. Then walk $m=0\ldots N-1$, skipping
already-used indices, and emit an index group:

$$\text{TRIM } (m'=m):\quad \big[\,2m,\ 2m{+}1,\ n{+}2m,\ n{+}2m{+}1\,\big] \quad (4\text{ indices})$$
$$\text{generic}:\quad \big[\,2m,\,2m{+}1,\,2m',\,2m'{+}1,\ n{+}2m,\,n{+}2m{+}1,\,n{+}2m',\,n{+}2m'{+}1\,\big]\quad (8)$$

Concatenate the groups into a permutation $P$, and apply it symmetrically:
$\mathcal H^{\rm bd}=\mathcal H^k[P][:,P]$.

**Why it works.** $\mathcal F^\dagger A\mathcal F$ is diagonal in $k$ (translation invariance);
$\mathcal F^\dagger B\mathcal F^*$ connects $k\leftrightarrow-k$ only (a zero-momentum pair). So
the sparsity graph of $\mathcal H^k$ splits into connected components labelled by the unordered
pair $\{k,-k\}$ — that is exactly what a permutation can expose, and nothing more.

**Output:** block-diagonal $\mathcal H^{\rm bd}$: $4\times4$ blocks at TRIMs ($k=0$, and $k=\pi$ if
$N$ even), $8\times8$ elsewhere. Block sizes sum to $2n$.

**Check:** Frobenius norm of everything outside the declared blocks $\sim10^{-15}$.

---

### Step 4 — The per-$k$ blocks in closed form

Define the Bloch factors

$$\delta_{C,k}=g_{c1}+g_{c2}e^{-ik},\qquad \delta_{J,k}=g_{J1}+g_{J2}e^{-ik},$$

and note $\delta_{C,-k}=\delta_{C,k}^*$ and $|\delta_{C,k}|^2=g_{c1}^2+g_{c2}^2+2g_{c1}g_{c2}\cos k$.

Carrying the sums $\sum_d \alpha_{\sigma\sigma'}(d)e^{-ikd}$ through gives, in the honest
convention,

$$\tilde A(k)=\tfrac12\begin{pmatrix}1&\delta_{C,k}\\ \delta_{C,k}^*&1\end{pmatrix},\qquad
  \tilde B(k)=-\tfrac12\begin{pmatrix}0&\delta_{C,k}\\ \delta_{C,k}^*&0\end{pmatrix},$$

$$\tilde A_J(k)=\tfrac12\begin{pmatrix}g_J&-\delta_{J,k}\\ -\delta_{J,k}^*&g_J\end{pmatrix},\qquad
  \tilde B_J(k)=-\tfrac12\begin{pmatrix}g_J&\delta_{J,k}\\ \delta_{J,k}^*&g_J\end{pmatrix}.$$

The generic-$k$ $8\times8$ sector, ordered $(c_{kA},c_{kB},c_{-kA},c_{-kB})$ then daggers, is

$$\mathcal H^{(k)}=\begin{pmatrix}\mathcal A&\mathcal B\\ \mathcal B^*&\mathcal A^*\end{pmatrix},\quad
\mathcal A=\begin{pmatrix}\tilde A(k)&0\\0&\tilde A(-k)\end{pmatrix},\quad
\mathcal B=\begin{pmatrix}0&\tilde B(k)\\ \tilde B(-k)&0\end{pmatrix}.$$

**All-annihilation / factor-2 form** (what `k_space_hamiltonians` returns; multiply by 2 and drop
the redundant creation half). With $\bar\delta\equiv\delta^*$:

$$A_C=\begin{pmatrix}1&\delta_C&&\\ \bar\delta_C&1&&\\ &&1&\bar\delta_C\\ &&\delta_C&1\end{pmatrix},\qquad
  B_C=\begin{pmatrix}&&&-\delta_C\\ &&-\bar\delta_C&\\ &-\bar\delta_C&&\\ -\delta_C&&&\end{pmatrix},$$

$$A_J=\begin{pmatrix}g_J&-\delta_J&&\\ -\bar\delta_J&g_J&&\\ &&g_J&-\bar\delta_J\\ &&-\delta_J&g_J\end{pmatrix},\qquad
  B_J=\begin{pmatrix}&&-g_J&-\delta_J\\ &&-\bar\delta_J&-g_J\\ -g_J&-\bar\delta_J&&\\ -\delta_J&-g_J&&\end{pmatrix}.$$

At a TRIM ($\sin k=0$, so $\delta$ real) the sector collapses to $2\times2$:

$$A_C=\begin{pmatrix}1&\delta_C\\ \bar\delta_C&1\end{pmatrix},\
  B_C=-\begin{pmatrix}0&\delta_C\\ \bar\delta_C&0\end{pmatrix},\
  A_J=\begin{pmatrix}g_J&-\delta_J\\ -\bar\delta_J&g_J\end{pmatrix},\
  B_J=-\begin{pmatrix}g_J&\delta_J\\ \bar\delta_J&g_J\end{pmatrix}.$$

**Three traps, all hit in the first hand-derivation:**

1. the TRIM $B_C$ carries $-\delta_C$, the same sign as the generic case;
2. the generic $B_C$ B-sublattice entries are $-\bar\delta_C$, **not** $-\delta_C$ — the anomalous
   coefficient of $c_{kB}^\dagger c_{-kA}^\dagger$ is $\tilde B_{BA}(k)+\tilde B_{AB}(-k)=-\delta_C^*$;
3. the $-g_J$ on-site Josephson pairing lives on the $k\leftrightarrow-k$ block-diagonal
   (entries $[0,2],[1,3]$ and transposes), **not** on the main diagonal — on-site pairing
   $b_l^\dagger b_l^\dagger$ Fourier-transforms to $c_{k\sigma}^\dagger c_{-k\sigma}^\dagger$. It
   only *looks* diagonal at a TRIM because there $k=-k$.

Also: the TRIM test must be $\sin k = 0$, not `k in {0, pi}` — with `fftshift` the second TRIM sits
at $k=-\pi$.

**Check:** these agree block-by-block with $\mathcal G^\dagger\mathcal H\mathcal G$ to $\lesssim10^{-15}$.

---

## Part II — Diagonalising $H_C$ sector by sector

Everything here is a **canonical (paraunitary) transform**: a matrix $\mathcal T$ with

$$\mathcal T^\dagger\tau_3\mathcal T=\tau_3,\qquad \tau_3=\operatorname{diag}(I_m,-I_m),$$

acting as $\mathcal H\to\mathcal T^\dagger\mathcal H\,\mathcal T$. ($m$ = 2 at a TRIM, 4 generic.)

### Step 5 — Phased sublattice rotation

**Input:** $k$, $g_{c1}$, $g_{c2}$. Set $\phi_k=\arg\delta_{C,k}$.

The naive symmetric/antisymmetric combination $(b_{kA}\pm b_{kB})/\sqrt2$ **fails** for generic $k$:
it leaves an $\alpha_k$–$\beta_k$ number-conserving coupling of magnitude $\tfrac12 g_{c2}|\sin k|$
(zero only where $\delta_{C,k}$ is real, i.e. at the TRIMs). The correct combination is the
eigenbasis of $\tilde A(k)$:

$$\boxed{\ \alpha_k=\tfrac1{\sqrt2}\big(e^{-i\phi_k/2}b_{kA}+e^{+i\phi_k/2}b_{kB}\big),\qquad
        \beta_k =\tfrac1{\sqrt2}\big(e^{-i\phi_k/2}b_{kA}-e^{+i\phi_k/2}b_{kB}\big)\ }$$

so that $\tilde A(k)\to\operatorname{diag}\!\big(\tfrac12(1+|\delta_{C,k}|),\ \tfrac12(1-|\delta_{C,k}|)\big)$.

As a matrix, rows $=(\alpha,\beta)$:

$$R(k)=\frac{1}{\sqrt2}\begin{pmatrix}e^{-i\phi_k/2}&e^{+i\phi_k/2}\\ e^{-i\phi_k/2}&-e^{+i\phi_k/2}\end{pmatrix},
\qquad R(-k)=R(k)\big|_{\phi_k\to-\phi_k},$$

and for a generic sector $M=\operatorname{blkdiag}\big(R(k),R(-k)\big)$ ($4\times4$).

**Nambu lift of a passive rotation.** If the new modes are $\mu=M b$ then $b=M^{-1}\mu$, so

$$\mathcal T_{\rm rot}=\begin{pmatrix}M^{-1}&0\\0&(M^{-1})^*\end{pmatrix}.$$

(Getting this dagger backwards silently produces a wrong, non-diagonalising rotation — worth an
explicit test that the number-off-diagonal really drops to $\sim10^{-16}$.)

---

### Step 6 — Two-mode squeeze

**Input:** the rotated sector.

After Step 5 each $\{k,-k\}$ pair is, per species $\gamma\in\{\alpha,\beta\}$,

$$\omega_\gamma\big(n_{k\gamma}+n_{-k\gamma}\big)+\lambda_\gamma\big(\gamma_k^\dagger\gamma_{-k}^\dagger+\text{h.c.}\big),$$
$$\omega_{\alpha}=\tfrac12(1+|\delta_{C,k}|),\ \ \lambda_{\alpha}=-\tfrac12|\delta_{C,k}|,\qquad
  \omega_{\beta}=\tfrac12(1-|\delta_{C,k}|),\ \ \lambda_{\beta}=+\tfrac12|\delta_{C,k}| .$$

A two-mode squeeze $S\gamma_iS^\dagger=\cosh z\,\gamma_i+\sinh z\,\gamma_j^\dagger$ has Nambu form

$$\mathcal T_{\rm sq}=\begin{pmatrix}U&V\\V^*&U^*\end{pmatrix},\quad
U_{ii}=U_{jj}=\cosh z,\ \ V_{ij}=V_{ji}=\sinh z,\ \ \text{rest } U=I,\ V=0 .$$

The generic $2\times2$ problem $\big(\begin{smallmatrix}\omega&\lambda\\ \lambda&\omega\end{smallmatrix}\big)$
is killed by $z=\tfrac14\ln\frac{\omega+\lambda}{\omega-\lambda}$, leaving $\sqrt{\omega^2-\lambda^2}$.
Substituting $\omega_\pm,\lambda_\pm$:

$$\boxed{\ z_{k,\alpha}=\tfrac14\ln\!\big(1+2|\delta_{C,k}|\big),\qquad
        z_{k,\beta}=\tfrac14\ln\!\big(1-2|\delta_{C,k}|\big)\ }$$

applied on the pairs $(\alpha_k,\alpha_{-k})$ and $(\beta_k,\beta_{-k})$. At a TRIM $k=-k$ and this
degenerates into a **single-mode** squeeze ($i=j$), same $z$.

The **factor 2 inside the log is essential**: $z_\beta$ requires $1-2|\delta_{C,k}|>0$, and since
$\max_k|\delta_{C,k}|=g_{c1}+g_{c2}$ that is exactly the stated existence bound

$$g_{c1}+g_{c2}<\tfrac12 .$$

(If one instead writes $\delta$ with the factor 2 absorbed, the same formulas read
$z=\tfrac14\ln(1\pm|\delta|)$ — the bound is unchanged, only the definition of $\delta$ moved.)

**Total transform:** $\mathcal T=\mathcal T_{\rm rot}\,\mathcal T_{\rm sq}$, and
$\mathcal T^\dagger\mathcal H_C^{(k)}\mathcal T$ is diagonal.

**Bands (physical convention, $\times2$):**

$$\boxed{\ \omega_{k,\alpha}=\sqrt{1+2|\delta_{C,k}|},\qquad \omega_{k,\beta}=\sqrt{1-2|\delta_{C,k}|}\ }$$

**Checks (all three must agree to $\sim10^{-15}$):**

1. $2\times$ the diagonal of $\mathcal T^\dagger\mathcal H_C^{(k)}\mathcal T$;
2. the positive eigenvalues of $\tau_3\,\mathcal H_C^{(k)}$ — the convention-independent definition
   of bosonic BdG frequencies, times the same factor 2;
3. `PHOTONICChain(PBC=True)._nambu_prefactor() * diag(chain.H_C)`.

Also check paraunitarity: $\|\mathcal T^\dagger\tau_3\mathcal T-\tau_3\|<10^{-9}$.

**The same $\mathcal T$ does not diagonalise $H_J$** — $\arg\delta_{J,k}\neq\phi_k$ — so

$$\tilde{\mathcal H}_J(k)=\mathcal T^\dagger\mathcal H_J^{(k)}\mathcal T
=\begin{pmatrix}\tilde{\mathcal A}_J&\tilde{\mathcal B}_J\\ \tilde{\mathcal B}_J^*&\tilde{\mathcal A}_J^*\end{pmatrix}$$

retains an $\alpha$–$\beta$ number mixing $\tilde{\mathcal A}_J$ *and* a pairing block
$\tilde{\mathcal B}_J$. These two objects drive Parts IV–V.

---

## Part III — Zak phase

### Step 7 — Discrete Berry phase (Wilson loop) with the orbital gauge

**Formula.** Use the single-log-of-the-product form — it is gauge invariant (each $|u_k\rangle$
appears once as bra and once as ket, so the arbitrary phases `eigh` assigns cancel) and folds the
answer into $(-\pi,\pi]$ automatically:

$$\boxed{\ \gamma=-\operatorname{Im}\ln\Big[\textstyle\prod_{j}\langle u(k_j)\,|\,u(k_{j+1})\rangle\Big]\ }$$

(Writing it as $-\prod_j\operatorname{Im}\ln\langle\cdot\rangle$ is wrong; the sum form
$-\sum_j\operatorname{Im}\ln$ is equivalent modulo $2\pi$ but exposed to branch cuts.)

**Grid.** $k_j = 2\pi j/N_k$, $j=0\ldots N_k-1$, endpoint excluded — $k=0$ and $k=2\pi$ are the
same point; the wrap is done by hand in the closure term.

**Inner product.** $\langle a|b\rangle=\sum_i \bar a_i b_i$ (conjugate the **bra**).

**Orbital gauge.** The loop is over the cell-periodic part $|u_k\rangle=e^{-ik\hat r}|\psi_k\rangle$
with $\hat r=\operatorname{diag}(s_A,s_B)$ the sublattice positions in the cell. Since
$u_{k+2\pi}=e^{-i2\pi\hat r}u_k$, the entire orbital winding can be lumped into the closing link:

$$\gamma=-\operatorname{Im}\ln\Big[\prod_{j=0}^{N_k-2}\langle\psi_j|\psi_{j+1}\rangle\ \cdot\
\langle\psi_{N_k-1}|\,\hat g\,|\psi_0\rangle\Big],\qquad
\hat g=\operatorname{diag}\big(e^{-2\pi i s_A},\,e^{-2\pi i s_B}\big).$$

Lumping is exact whenever the sublattice weights $|\psi_{k,\sigma}|^2$ are $k$-independent (true
for a chiral SSH band); otherwise distribute $e^{-i\Delta k\,\hat r}$ over every link instead
(difference is $O(1/N_k)$).

**States.** $|\psi_k\rangle$ = eigenvector of the SSH Bloch Hamiltonian

$$h(k)=\begin{pmatrix}0&\delta_{C,k}\\ \delta_{C,k}^*&0\end{pmatrix} \;=\; 2\tilde A(k)-I,$$

lower band = eigenvalue $-|\delta_{C,k}|$, upper band = $+|\delta_{C,k}|$. (Same eigenvectors as
$\tilde A(k)$; the squeeze never rotates sublattice space, so these are also the quasiparticle
sublattice spinors.)

**Result.** For the chiral SSH band $|\psi_{k,A}|^2=|\psi_{k,B}|^2=\tfrac12$, so the orbital term
contributes $-\pi(s_A+s_B)$:

| $s=(s_A,s_B)$ | $\gamma$ at $\chi_C>0$ | $\gamma$ at $\chi_C<0$ |
|---|---|---|
| $(0,0)$ periodic gauge | $0$ | $\pi$ |
| $(-\tfrac14,\tfrac14)$ symmetric cell | $0$ | $\pi$ |
| $(0,\tfrac12)$ | $\pi/2$ | $-\pi/2$ |

Exactly quantised and grid-independent ($N_k=10,13,50,400$ identical); transition at $\chi_C=0$,
where $g_{c1}=g_{c2}$ and the gap closes at $k=\pi$.

---

### Step 8 — Does the pairing change it? (BdG Wilson loop)

**Question:** the squeeze mixes $b$ and $b^\dagger$; does the Zak phase of the *quasiparticle* band
differ from that of the bare SSH band?

**Setup.** Single-$k$ BdG Bloch Hamiltonian on $\Psi_k=(b_{kA},b_{kB},b_{-kA}^\dagger,b_{-kB}^\dagger)$,
physical convention:

$$\mathcal H_{\rm BdG}(k)=\begin{pmatrix}h(k)&\Delta(k)\\ \Delta(k)&h(k)\end{pmatrix},\qquad
  h(k)=\begin{pmatrix}1&\delta_{C,k}\\ \delta_{C,k}^*&1\end{pmatrix},\qquad
  \Delta(k)=I-h(k).$$

Diagonalise the **non-Hermitian** $\tau_3\mathcal H_{\rm BdG}(k)$ ($\tau_3=\operatorname{diag}(1,1,-1,-1)$),
pick the eigenvalue nearest $+\sqrt{1-2|\delta_{C,k}|}$ (the $\beta$ band), and normalise with the
**para-norm**

$$\Phi^\dagger\tau_3\Phi=+1\qquad(\text{i.e. }|u|^2-|v|^2=1),$$

not the ordinary norm. Run the same Wilson loop on the 4-component $\Phi_\beta(k)$, with the hole
block taking the *opposite* orbital phase:

$$\hat g = \operatorname{diag}\big(e^{-2\pi i s_A},\,e^{-2\pi i s_B},\ e^{+2\pi i s_A},\,e^{+2\pi i s_B}\big),$$

because $v_\sigma$ transforms like $b_{-k,\sigma}^\dagger$. (Magnitudes drop out of
$\operatorname{Im}\ln$, so the non-unit ordinary norm is harmless.)

**Answer: unchanged, to $\sim10^{-12}$, in every orbital gauge and on both sides of the
transition.** Reason: $\Delta(k)=I-h(k)$ **commutes with** $h(k)$, so the Bogoliubov rotation is
(phased sublattice rotation) $\times$ (real squeeze on the $\{k,-k\}$ pair). The particle block
$(u_{kA},u_{kB})$ is exactly the SSH eigenvector, the hole block points the same way with a real,
non-winding amplitude $\sinh z_{k,\beta}$, and chiral symmetry survives
($|u_{kA}|^2=|u_{kB}|^2$ and $|v_{kA}|^2=|v_{kB}|^2$ at every $k$). Nothing in the squeeze can wind.

---

## Part IV — Schrieffer–Wolff correction

### Step 9 — What is actually eliminated

Add the **cavity resonator** (the mode mediating the inductive coupling — not a readout resonator):

$$H=\omega_r a^\dagger a+H_C-\gamma\,(a+a^\dagger)\,H_J .$$

`PHOTONICChain.schrieffer_wolff_correction` **does not eliminate the cavity.** It eliminates only
the *matter-non-number-conserving* part of the coupling — the pairing blocks
$\tilde{\mathcal B}_J$, which would take the single quasiparticle into the 3-QP sector — and folds
it into a 2nd-order term living in the **same composite space** as $H$,

$$(\text{cavity Fock }d_{\rm ph})\ \otimes\ (\text{single-QP matter }d_X),$$

added as `H = H + H_SW`. The number-conserving $\tilde{\mathcal A}_J$ stays in $H$ at first order.

**Ingredients, per $\{k,-k\}$ sector** (physical convention, $m$ = 4 generic / 2 TRIM QP modes):

* QP frequencies $w_i = 2\,[\mathcal T^\dagger\mathcal H_C^{(k)}\mathcal T]_{ii}$, $i=1\ldots m$,
  ordered $(\alpha_k,\beta_k,\alpha_{-k},\beta_{-k})$;
* pairing block $\mathcal B \equiv 2\,[\tilde{\mathcal H}_J(k)]_{[:m,\,m:]}$, giving the operator
  $\hat B=\sum_{ij}\mathcal B_{ij}\,\gamma_i^\dagger\gamma_j^\dagger$;
* single-QP basis $|X\rangle=\gamma_X^\dagger|0\rangle$ with $E_X=w_X$ ($d_X=m$ states);
* virtual basis: **3-QP** occupation states $|r\rangle$ with $E_r=\sum_i n_i^{(r)}w_i$
  (the $-2$ branch would need $-1$ excitations and is skipped);
* matrix elements $V_{rX}=\langle r|\hat B|X\rangle$, evaluated in the bosonic Fock basis;
* photon ladder $s = \gamma(a+a^\dagger)$, i.e.
  $s_{n,n'}=\gamma\big(\sqrt{n'}\,\delta_{n,n'-1}+\sqrt{n'+1}\,\delta_{n,n'+1}\big)$.

**The formula.** With composite index $(n,X)\to n\,d_X+X$ (photon outer, matter inner):

$$H_{\rm SW}\big[(n_1,X),(n_2,X')\big]=\sum_{n_{\rm mid}} s_{n_1 n_{\rm mid}}\,s_{n_{\rm mid} n_2}\;
\Big(\mathrm{T1}_{XX'}-\mathrm{T2}_{XX'}\Big),$$

$$\mathrm{T1}_{XX'}=\sum_r\frac{V_{rX}\,\overline{V_{rX'}}}{(n_1-n_{\rm mid})\,\omega_r+E_X-E_r},
\qquad
\mathrm{T2}_{XX'}=\sum_r\frac{V_{rX}\,\overline{V_{rX'}}}{(n_{\rm mid}-n_2)\,\omega_r+E_r-E_{X'}},$$

then symmetrise and halve:

$$H_{\rm SW}\ \leftarrow\ \tfrac12\big(H_{\rm SW}+H_{\rm SW}^\dagger\big),\qquad
  H_{\rm SW}\ \leftarrow\ \tfrac12\,H_{\rm SW}.$$

Guard denominators with $|{\cdot}|>\varepsilon$ (drop the term otherwise). Because $s_{nn'}$ is
purely off-diagonal in photon number, $s_{n_1n_{\rm mid}}s_{n_{\rm mid}n_2}\neq0$ only for
$|n_1-n_2|\in\{0,2\}$ — the correction has a dispersive ($\Delta n_{\rm ph}=0$) piece and a
cavity two-photon squeezing ($\Delta n_{\rm ph}=\pm2$, $\propto a^{\dagger2},a^2$) piece.

**Exact $\gamma^2$ scaling:** $\gamma$ enters only through $s$, quadratically. So
$H_{\rm SW}/\gamma^2$ is $\gamma$-independent by construction — build it once at $\gamma=1$ and
rescale.

### Step 9b — The split, and why the block structure is simpler than expected

> **Observation (read off directly from the composite matrix plots).**
> $H_{\rm SW}$ has **no $+k\to-k$ matrix element at all** — those entries are identically zero, to
> machine precision, at every $k$ and for every parameter set tried. Its off-diagonal part is
> *purely* $\alpha_k\!\leftrightarrow\!\beta_k$ **at the same momentum**.
>
> So "$k$-changing" is a misnomer for it: the correction never changes $k$, it only flips the band.
> A better label is **interband**.

The immediate consequence is that the block structure of $H_{\rm SW}$ is **finer than the
$\{k,-k\}$ sector** the rest of the calculation is organised around. Inside one $8\times8$ sector,
$H_{\rm SW}$ does not connect the $+k$ half to the $-k$ half at all, so it further decomposes into
two *independent* blocks, one per momentum:

$$H_{\rm SW}^{(\{k,-k\})}
=H_{\rm SW}^{(+k)}\ \oplus\ H_{\rm SW}^{(-k)},\qquad
H_{\rm SW}^{(\pm k)}\ \text{acting on}\ \{\alpha_{\pm k},\beta_{\pm k}\}\otimes(\text{cavity Fock}),$$

each of dimension $2d_{\rm ph}$. Nothing couples the two. (Sanity check to reproduce:
$\|H_{\rm SW}[\,\{\alpha_k,\beta_k\},\ \{\alpha_{-k},\beta_{-k}\}\,]\|=0$ exactly — not small, zero.)

Within one such block the split is then simply diagonal vs off-diagonal in the band index:

| | matrix elements | share of $\|H_{\rm SW}\|$ |
|---|---|---|
| **band-diagonal** ("$k$-preserving") | $\gamma_{k\gamma}^\dagger\gamma_{k\gamma}\otimes(\cdot)_{\rm cav}$ | $\approx99\%$ |
| **band-off-diagonal** ("$k$-changing", really *interband*) | $\gamma_{k\alpha}^\dagger\gamma_{k\beta}\otimes(\cdot)_{\rm cav}$ | $10^{-3}$–$10^{-2}$ |

The interband part is small because $\tilde{\mathcal B}_J$ predominantly pairs $\alpha_k$ with
$\beta_{-k}$; only its same-band entries feed a band flip.

**Where it comes from — the structure of $\tilde{\mathcal H}_J$.** The whole SW correction is
assembled from $\tilde{\mathcal H}_J(k)=\mathcal T^\dagger\mathcal H_J^{(k)}\mathcal T$, and that
matrix has a rigid selection rule, unchanged over the whole $(\Omega_C,\chi_C)$ plane:

* the number block $\tilde{\mathcal A}_J$ is **momentum block-diagonal** — only $\{+k\}\times\{+k\}$
  and $\{-k\}\times\{-k\}$ entries. It flips the band but never the momentum.
* the pairing block $\tilde{\mathcal B}_J$ is **momentum block-anti-diagonal** — only
  $\{+k\}\times\{-k\}$ entries. Every $\hat B$ creates exactly one $+k$ and one $-k$ quasiparticle,
  a **zero-momentum pair**.

So each insertion of the coupling either leaves $k$ alone or adds/removes a zero-momentum pair; a
single quasiparticle at $+k$ can never emerge at $-k$. Consistent with (and more concrete than)
crystal-momentum conservation — the cavity carries no momentum. Whether the cancellation is *exact*
for a deeper reason, whether it survives at higher order, and what it implies for the effective
theory is left open here; see the open points at the end. A dedicated notebook cell (just before the cavity-dressed-bands analysis) plots
$|\tilde{\mathcal H}_J(k)|$ for several $(\Omega_C,\chi_C)$ and confirms the checkerboard is
parameter-independent.

### Step 9c — Interband resonance condition

The interband gap and the resonant momentum, solved in closed form. From
$\sqrt{1+2x}-\sqrt{1-2x}=\omega_r$ (square twice):

$$x\equiv|\delta_{C,k^\star}|=\frac{\omega_r}{2}\sqrt{1-\frac{\omega_r^2}{4}},\qquad
\cos k^\star=\frac{x^2-g_{c1}^2-g_{c2}^2}{2g_{c1}g_{c2}} .$$

The gap sweeps $\big[\,\mathrm{gap}(\pi),\ \mathrm{gap}(0)\,\big]$ across the BZ with
$\mathrm{gap}(0)=\sqrt{1+2(g_{c1}+g_{c2})}-\sqrt{1-2(g_{c1}+g_{c2})}\approx 2(g_{c1}+g_{c2})=\Omega_C$.
Choosing $\omega_r=2(g_{c1}+g_{c2})$ puts $k^\star$ just inside the zone
(for the reference parameters, $k^\star=0.2315=0.0737\pi$).

### Step 9d — Sweep of the interband ("$k$-changing") part

For each $(\chi_C,\Omega_C)$ with $\omega_r=\Omega_C$, sweep $k$ over the BZ and record
$\max_k\big|H_{\rm SW}^{\,k\text{-chg}}(k)\big|/\gamma^2$ (log colour scale). Findings:

* it **vanishes identically on $\chi_C=0$** ($\lesssim10^{-14}$) — there $g_{c1}=g_{c2}$,
  $\tilde{\mathcal B}_J$ is exactly band-diagonal, and $\alpha_k\to\beta_k$ is forbidden;
* dimerisation switches it on; grows monotonically with $|\chi_C|$ and $\Omega_C$, symmetric under
  $\chi_C\to-\chi_C$;
* peak $\approx3.0$ (in $\gamma^{-2}$ units) at $(\chi_C,\Omega_C)=(-0.70,0.64)$.

---

## Part V — Cavity-dressed bands and the topological second band

### Step 10 — The effective Hamiltonian per $k$

Restrict to the $+k$ quasiparticles. By Step 9b this is an **exact block decomposition, not a
truncation**: neither the direct coupling nor $H_{\rm SW}$ has a single matrix element linking
$+k$ to $-k$, so the $-k$ half simply factors out. Work on

$$\{\alpha_k,\beta_k\}\otimes(\text{cavity Fock }d_{\rm ph}),\qquad
\text{index } = n\cdot 2+s,\ \ s\in\{0=\alpha,\,1=\beta\},$$

a $2d_{\rm ph}$-dimensional space. Then

$$\boxed{\
H_{\rm eff}(k)=\underbrace{\omega_r\,a^\dagger a+\operatorname{diag}(\omega_{k\alpha},\omega_{k\beta})
 -\gamma\,(a+a^\dagger)\otimes M}_{\texttt{build\_full\_hamiltonian},\ \text{single-exc subspace}}
 \;+\;H_{\rm SW}^{+k}(k)\ }$$

with

* $M=\tilde{\mathcal A}_J+\tilde{\mathcal A}_J^\dagger$ restricted to the $+k$ band block — this
  reproduces `build_HJ_fock`'s coefficient $A_J+D_J^{T}$, since the lower Nambu block satisfies
  $D_J=\mathcal A_J^*$;
* $(a+a^\dagger)_{nn'}=\sqrt{n'}\,\delta_{n,n'-1}+\sqrt{n'+1}\,\delta_{n,n'+1}$;
* $H_{\rm SW}^{+k}$ = Step 9's $H_{\rm SW}$ restricted to rows/cols whose QP mode index is
  $\alpha_k$ or $\beta_k$ — **all of it**, band-diagonal *and* interband, since by Step 9b nothing
  in $H_{\rm SW}$ leaves the $+k$ block.

$H_{\rm eff}$ is Hermitian; diagonalise it normally.

**What opens the gap:** the direct term $-\gamma(a+a^\dagger)\otimes M$, whose interband element
drives a Jaynes–Cummings anticrossing between $|\alpha_k,0\rangle$ and $|\beta_k,1\rangle$ exactly
at $k^\star$. Gap size $\approx 2\gamma\,|M_{\alpha\beta}(k^\star)|$.

### Step 11 — Composite inner product and Zak phase

An eigenstate is

$$|\psi_n(k)\rangle=|\phi^{(n)}_\alpha(k)\rangle_{\rm cav}\otimes|\alpha_k\rangle
                   +|\phi^{(n)}_\beta(k)\rangle_{\rm cav}\otimes|\beta_k\rangle .$$

The band kets are $k$-dependent, so overlaps at different $k$ do **not** reduce to the eigenvector
components alone:

$$\langle\psi_n(k)|\psi_n(k')\rangle=\sum_{s,s'}
\underbrace{\langle u_s(k)|u_{s'}(k')\rangle_{\rm sub}}_{2\times2}\;
\underbrace{\langle\phi^{(n)}_s(k)|\phi^{(n)}_{s'}(k')\rangle_{\rm cav}}_{\text{Fock, }k\text{-indep}} .$$

**Practical recipe — lift into a fixed basis.** Reshape the eigenvector to
$C^{(n)}(k)\in\mathbb C^{2\times d_{\rm ph}}$ with $C_{s,n_{\rm ph}}$, take
$U(k)=[\,u_\alpha(k)\,|\,u_\beta(k)\,]$ (columns = the SSH Bloch spinors of Step 7), and set

$$\boxed{\ \Psi_n(k)=U(k)\,C^{(n)}(k)\ \in\ \mathbb C^{2}_{\rm sub}\otimes\mathbb C^{d_{\rm ph}}_{\rm cav}\ }$$

flattened sublattice-outer. Then $\langle\psi_n(k)|\psi_n(k')\rangle=\Psi_n(k)^\dagger\Psi_n(k')$ —
an ordinary inner product. The BZ-boundary factor acts on the **sublattice half only** (the cavity
does not wind):

$$\hat g=\operatorname{diag}\big(e^{-2\pi i s_A},e^{-2\pi i s_B}\big)\otimes\mathbb 1_{d_{\rm ph}}
\quad\Longrightarrow\quad \hat g_{\rm flat}=\texttt{repeat}\big([e^{-2\pi i s_A},e^{-2\pi i s_B}],\,d_{\rm ph}\big),$$

$$\gamma^{(n)}_{\rm Zak}=-\operatorname{Im}\ln\Big[\prod_{j=0}^{N_k-2}\Psi_n(k_j)^\dagger\Psi_n(k_{j+1})\ \cdot\
\Psi_n(k_{N_k-1})^\dagger\big(\hat g_{\rm flat}\odot\Psi_n(k_0)\big)\Big].$$

**Multiband version** (for a near-degenerate set $\mathfrak B$, more robust than any single band):
build the link matrices $\mathcal M^{(j)}_{nn'}=\Psi_n(k_j)^\dagger\Psi_{n'}(k_{j+1})$, take
$W=\prod_j\mathcal M^{(j)}$ (with $\hat g$ on the last link) and $\gamma=-\operatorname{Im}\ln\det W$.

**Which band.** Sort eigenvalues ascending. Band 0 is the zero-photon lower state
$|\beta_k,0\rangle$; the **second band** is index 1, the lower polariton of the anticrossing.

### Step 12 — Cavity purity (diagnostic colouring)

Trace the band degree of freedom out of the eigenstate. With $C\in\mathbb C^{d_{\rm ph}\times2}$
indexed $[\,n_{\rm ph},\,s\,]$,

$$\rho_{\rm ph}=C\,C^\dagger \ (d_{\rm ph}\times d_{\rm ph}),\qquad
\text{purity}=\operatorname{Tr}\rho_{\rm ph}^2=\sum_{m,m'}\big|(\rho_{\rm ph})_{mm'}\big|^2 .$$

Range $[\tfrac12,1]$: $1$ = pure Fock (band and cavity unentangled, a clean photon replica),
$\tfrac12$ = maximally band–photon hybridised (right at an anticrossing). Plot
$\log_{10}(1-\operatorname{Tr}\rho_{\rm ph}^2)$ (clipped at $\sim10^{-6}$) to resolve the nearly-pure
regions.

---

## Results

### Reference parameters

$N=10$, $\Omega_C=0.2$, $\chi_C=0.5$, $\Omega_J=1.0$, $\chi_J=0$
$\Rightarrow g_{c1}=0.075,\ g_{c2}=0.025,\ g_{J1}=1,\ g_{J2}=-1,\ g_J=0$.
$\omega_r=2(g_{c1}+g_{c2})=0.2$, $\gamma=0.03$, $d_{\rm ph}=4$–$6$, $s=(-\tfrac14,\tfrac14)$.
Interband gap window $[\mathrm{gap}(\pi),\mathrm{gap}(0)]=[0.100,\,0.201]$; $k^\star=0.2315$.

### Zak phases

| object | $\gamma_{\rm Zak}$ |
|---|---|
| bare SSH lower band, $\chi_C>0$ / $\chi_C<0$ | $0$ / $\pi$ |
| BdG $\beta$ quasiparticle band (pairing included) | identical, to $10^{-12}$, in every gauge |
| dressed band 0 — $|\beta_k,0\rangle$ | $0$ |
| **dressed band 1 — second band (lower polariton)** | $\boldsymbol\pi$ |
| dressed band 2 — upper polariton | $\pi$ |
| dressed bands $\{1,2\}$ combined (det Wilson loop) | $0$ |

* $\gamma=0\Rightarrow$ band 1 is trivial; **any $\gamma>0\Rightarrow\pi$**. The cavity gaps the
  $|\alpha_k,0\rangle$–$|\beta_k,1\rangle$ crossing and the second band becomes topological — the
  photon-induced edge state.
* The invariant tracks the SSH phase: $\chi_C>0$ gives (band 0 trivial, band 1 topological),
  $\chi_C<0$ the reverse. The dressing hands the invariant from band 0 to band 1.

### Cavity frequency as a topological knob

$$\gamma_{\rm Zak}(\text{band 1})=\pi\ \ \text{for }\omega_r\lesssim\mathrm{gap}(0),
\qquad =0\ \ \text{for }\omega_r>\mathrm{gap}(0),$$

a clean step at $\omega_r\approx0.207$. Above it no $k$ satisfies the resonance condition, band 1
becomes pure $|\alpha_k,0\rangle$, and the invariant reverts to the trivial bare-SSH value. The
combined $\{1,2\}$ invariant stays $0$ throughout — the $\pi$ is genuinely split between the two
polaritons. Robust across $d_{\rm ph}=3$–$6$ and $N_k$.

### Photon replicas and the mechanism

The uncoupled ladder is the two SSH bands replicated at $\omega_{k\alpha/\beta}+m\omega_r$; the
dressed bands follow them except at replica crossings, where they anticross. Colouring by cavity
purity shows the mechanism directly:

* band 0 stays essentially pure across the whole BZ;
* the second band runs pure $|\beta_k,1\rangle$ near the zone centre $\to$ 50/50 hybridised at
  $k^\star\approx\pm2$ (exactly at the $\beta$/$\alpha$ replica crossing) $\to$ pure
  $|\alpha_k,0\rangle$ near $k=\pi$. **That character swap across the BZ is the $\pi$ Zak phase.**
* sharp purity spikes sit at $k=0,\pm\pi$ — at the TRIMs the interband coupling
  $M_{\alpha\beta}\to0$ and the bands cannot hybridise at all.

---

## Caveats

1. **The SW correction is not a light dressing.** With $\Omega_J=1$ the pairing block
   $\tilde{\mathcal B}_J$ is $O(1)$, so at $\gamma=0.03$ one gets $\|H_{\rm SW}\|\sim0.1$–$0.5$,
   comparable to the bandwidth, and it grows with $d_{\rm ph}$ because the $\Delta n_{\rm ph}=\pm2$
   terms carry $\sqrt{n}$. Conclusions were checked stable in $d_{\rm ph}$ and $N_k$, but the
   genuinely perturbative regime is $\gamma\sim0.01$, where the with/without-SW spectra are nearly
   identical.
2. **The band 1 / band 2 gap at $k^\star$ is small** ($10^{-3}$–$10^{-2}$, set by
   $\gamma|M_{\alpha\beta}(k^\star)|$, and $k^\star$ sits near a TRIM where $M_{\alpha\beta}$ is
   itself small). The single-band invariant is therefore more delicate than the combined
   $\{1,2\}$ one. Hardening it means widening the gap: larger $\gamma$, or an $\omega_r$ that puts
   $k^\star$ further from the TRIM.
3. **Not done yet:** an open-boundary real-space calculation to exhibit the predicted edge state and
   confirm bulk–boundary correspondence for the dressed second band.
4. **Open (being explored):** why $H_{\rm SW}$ has *exactly* zero $+k\leftrightarrow-k$ weight.
   Momentum conservation explains that a single quasiparticle cannot change $k$, but it is worth
   checking directly whether that is the full mechanism, whether the same holds for the
   number-conserving coupling and at higher order, and what it implies for the block structure of
   the effective theory.

---

## Map to the notebook

| object | function |
|---|---|
| Step 1 real-space Nambu blocks | `real_space_nambu` |
| Step 2 FFT momenta, $\pm k$ partner | `k_grid`, `neg_k_index` |
| Step 2 operator Fourier matrix $\mathcal F$, $\mathcal G^\dagger\mathcal H\mathcal G$ | `fourier_operator_matrix`, `to_k_space_nambu` |
| Step 3 permutation + verification | `sector_permutation`, `offblock_norm` |
| Step 4 closed-form per-$k$ blocks | `k_space_hamiltonians`, `reduced_k_blocks`, `sector_nambu_at_k` |
| Steps 5–6 rotation, squeeze, paraunitarity | `phased_sublattice_rotation`, `squeeze_nambu`, `passive_nambu`, `diagonalize_C_sector`, `para_defect`, `bogoliubov_freqs` |
| Step 7 Wilson loop | `compute_zak_phase`, `ssh_lower_vec` |
| Step 8 BdG spinor | `H_bdg_bloch`, `beta_bogoliubov_spinor` |
| Step 9 SW, composite | `H_SW_composite` |
| Step 9 SW, cavity-eliminated variant ($\delta A,\delta B$) | `sw_correction` |
| Step 9d parameter sweep | `peak_kchg` |
| Steps 10–11 dressed bands + composite Zak | `ssh_bloch_U`, `H_eff_dressed`, `H_eff_noSW`, `dressed_state`, `zak_dressed` |
| Step 12 purity colouring | `spec_purity` |
