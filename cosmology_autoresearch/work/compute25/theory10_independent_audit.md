# Independent audit of `work/theory10` early-time IVS IC report

## Disposition

**REVISE (leading algebra passes).** The background signs and powers, Eq. (19) reduction, regular \(D=0,\ g<0\) mode, and equal-ζ sign all check independently. Three presentation corrections are needed: the \(D>0\) series remainder is too strong if \(A\) denotes its full background value; the edge-branch result should be attributed to the displayed CDM perturbation system rather than to a fully specified vacuum closure; and the anisotropic-stress caveat can be stated more precisely. None changes the leading result.

## Sources, scope, and reproducibility

I read the primary PDF, `context/papers/interacting_de_desi_dr2_2026.pdf`, pp. 2–4: background Eqs. (4)–(5) are on p. 2, interaction Eq. (8) on p. 3, and IVS perturbation Eqs. (19)–(20) on p. 4. The paper defines δ_i=δρ_i/ρ_i, velocity divergence θ_i, synchronous-gauge trace \(h\), and conformal-time prime. It says Γ has Hubble-rate units and Γ>0 transfers energy from CDM to DE. The passage introducing zero anisotropic stress and zero momentum transfer is attached to its \(w_x\ne-1\) equations; the IVS equations are then given separately. For IVS the text says DE has no perturbative contribution, but supplies no independent δ_x evolution equation or explicit covariant transfer-perturbation completion.

Input hashes (SHA-256):

- Paper PDF: `0ca9dfa076a493eda6b4ab8fd7e5535ed410aa4d2e74f70442d835dd8d262a0c`
- Extracted source text: `a99a4afdd858cdce3cbc2c41f0a13098c30c361dc67891331ee7d9d5f824868f`
- Reviewed report: `302fa75cca4bf58f2dab8dbdb45f4aa671dd982ab1e0a463c626c2e219df7ef9`
- Reviewed script: `6e2726dd45e5e136226321ac5534cec870d71b2a8228c34ffd6b9bfd312dbe34`
- Reviewed JSON: `78ff37c9f7189fc23e95296083e8c85d588d5508fd69a88defddc0b53a3faa38`

This is an independent analytic derivation from the displayed equations, not an import or execution of the theory10 script. No Boltzmann solver was run. The independent derivation did not need a numerical integration; all conclusions below are leading asymptotic identities. The pinned SI radiation normalization in theory10 is a massless-neutrino *early-time* approximation, while the paper specifies Σmν=0.06 eV. The report labels this distinction and does not use the normalization to claim a full finite-redshift reproduction.

## Independent derivation

Write \(n=\ln a\), \(E=H/H_0\), \(g=\Gamma/H_0\), and Ω_r for the nonzero early relativistic density. In a flat radiation-era expansion, \(E\sim\sqrt{\Omega_{r0}}a^{-2}\). From the paper's proper-time background equations, with ρ_x>0 and \(X_* = \lim_{a\to0}\rho_x>0\),

\[
\frac{d\ln\rho_x}{dn}=\frac{\Gamma}{H}=\frac{g}{E}
\sim \frac{g}{\sqrt{\Omega_{r0}}}a^2,
\quad
\rho_x=X_*\left(1+\frac{g}{2\sqrt{\Omega_{r0}}}a^2+O(a^3)\right).
\]

The \(O(a^3)\) term can arise because the matter correction to \(E/(\sqrt{\Omega_{r0}}a^{-2})\) is \(O(a)\). For \(q=a^3\rho_c\), the CDM equation gives

\[
q_{,n}=-\frac{g}{E}a^3\rho_x
\sim-\frac{gX_*}{\sqrt{\Omega_{r0}}}a^5,
\quad
q=D-\frac{gX_*}{5\sqrt{\Omega_{r0}}}a^5+\cdots,
\quad D=\lim_{a\to0}q.
\]

Thus \(D<0\) is incompatible with positive early CDM. For \(D>0\), either sign of \(g\) has \(\rho_c\sim D a^{-3}\), and

\[
A\equiv\frac{\Gamma}{H}\frac{\rho_x}{\rho_c}
\sim K_5a^5,
\qquad K_5=\frac{gX_*}{\sqrt{\Omega_{r0}}D}.
\]

For \(D=0\), \(g>0\) gives negative leading ρ_c; \(g=0\) gives no CDM and hence no defined δ_c. For \(D=0,\ g<0\),

\[
\rho_c\sim-\frac{gX_*}{5\sqrt{\Omega_{r0}}}a^2>0,
\quad \beta_c\equiv\frac{d\ln\rho_c}{dn}\to2,
\quad A\to-5.
\]

For perturbations, \(\mathcal H=aH\). With the regular CDM-comoving choice \(\theta_c=0\), converting paper Eq. (19) from conformal time to \(n\) yields

\[
\delta_{c,n}=-\frac{h_{,n}}2+A\delta_c.
\]

The conversion coefficient is \(a\Gamma/\mathcal H=\Gamma/H=g/E\). Γ and H have proper-time inverse units; \(\mathcal H\), θ, \(h'\), and \(a\Gamma\) carry conformal-time inverse units; \(g/E\), A, h, and δ are dimensionless. Eq. (20) becomes \(\theta_{c,n}=-\theta_c\), so \(\theta_c\propto a^{-1}\); \(\theta_c=0\) is invariant and the nonzero solution is singular as \(a\to0\).

In the leading radiation metric \(h=C a^2\), \(S=h_{,n}/2=h\). On the ordinary branch the leading equation has \(A_5=K_5a^5\). Its regular solution begins δ_c=−h/2. Substituting δ_c=−h/2+cA_5h gives, at order \(A_5h\), \(7c=−1/2\), hence

\[
\delta_c=-\frac h2-\frac{A_5h}{14}+\cdots.
\]

This confirms the reported \(O(a^7)\) leading interaction correction. The displayed script tests the ideal truncated \(A=A_5\propto a^5,\ h\propto a^2\) equation, where the next residual is \(O(A_5^2h)\).

On the edge branch, \(A\to-5\). The homogeneous solution of Eq. (19) is \(a^{-5}\), while trying δ_c=C h gives \(2C=-1-5C\); the unique regular leading particular solution is δ_c=−h/7. For separately conserved radiation, β_γ=−4 and the standard growing synchronous radiation mode has δ_γ=−2h/3. The component curvature convention

\[
\zeta_i=-\psi-\frac{\delta_i}{\beta_i},\qquad
\beta_i=\frac{d\ln\rho_i}{d\ln a},
\]

then gives δ_γ/β_γ=h/6. Equal ζ_c=ζ_γ would require δ_c=β_c h/6=+h/3. The regular Eq. (19) solution instead gives

\[
\frac{\zeta_c-\zeta_\gamma}{h}
=-\left(-\frac1{14}-\frac16\right)=\frac5{21}.
\]

Equivalently, inserting δ_c=+h/3 into the leading edge equation leaves residual \((\delta_{c,n}+h-A\delta_c)/h=10/3\); inserting δ_c=−h/7 gives zero. These signs and coefficients match theory10.

The same time-shift argument for vacuum would give δ_x=β_x h/6, with β_x=g/E=O(a^2), hence δ_x=O(a^4) for nonzero \(g\). This is a conditional consistency observation only: the paper says IVS DE has no perturbative contribution but does not define an independent δ_x equation or a strict all-species ζ prescription. For \(g=0\), β_x=0 and ζ_x is undefined, not an additional constraint.

## Corrections to theory10

1. **Clarify the remainder in the ordinary-branch series.** The reported
   \(\delta_c=-h/2-Ah/14+O(A^2h)\) is valid for the *truncated leading system* \(A=A_5=K_5a^5,\ h\propto a^2\). If A means its full background value, its subleading expansion generally contains a \(K_6a^6\) term because matter is an \(O(a)\) correction to radiation in E. The first-order response to that term has order \(a^8\), larger than \(A_5^2h=O(a^{12})\); simply putting full A over 14 does not account for its changed logarithmic slope. Recommended exact wording: “\(\delta_c=-h/2-A_5h/14+O(a^8)\) for the leading radiation/matter asymptotic hierarchy; in the idealized pure-\(A_5\), fixed-\(h\propto a^2\) system the next interaction iteration is \(O(A_5^2h)\).” If theory10 intended the latter truncation throughout, explicitly define A in that formula as \(A_5\), not the full A defined earlier. This is a remainder-notation correction, not a change to the \(O(a^7)\) result.

2. **Narrow the edge-branch closure attribution.** The paper supplies Eq. (19) and (20) and says IVS DE has no perturbative contribution; it does not specify a complete perturbation of the transfer four-vector or an independent vacuum perturbation variable. Replace “the paper's stated vacuum-homogeneous/CDM-geodesic closure does not admit…” with: “Under the displayed IVS CDM equation, the regular (D=0,g<0) solution fails the usual CDM–radiation equal-ζ condition. This is an incompatibility of that displayed perturbation equation with that imposed IC criterion on the tuned edge, not a demonstrated inconsistency of a fully specified covariant IVS.” Keep the conditional δ_x caveat as written.

3. **Sharpen the anisotropic-stress note.** The paper's zero-anisotropic-stress phrase immediately introduces its \(w_x\ne-1\) system. It does not clearly impose zero neutrino shear on the separate IVS radiation hierarchy. Theory10 is right not to suppress free-streaming neutrino shear; phrase this as “the paper does not specify a total-radiation stress closure for IVS” rather than implying that the IVS clause itself is ambiguous about whether a stated assumption applies globally.

The report already correctly limits the standard mode to a leading radiation-era/fixed-metric statement, says the edge's global positivity is not checked, labels the vacuum relation conditional, distinguishes the paper's (0.06) eV sum from the asymptotically relativistic limit, and disclaims a full Einstein–Boltzmann or stability result. Those caveats should remain.

## Execution and limits

Read-only review commands included `pdftotext -f 2 -l 2 -layout context/papers/interacting_de_desi_dr2_2026.pdf -`, `pdftotext -f 3 -l 5 -layout context/papers/interacting_de_desi_dr2_2026.pdf -`, `sed -n` on the report/script/JSON, and `sha256sum` on the paper and reviewed artifacts. Python available in the environment reports 3.12.3. Final JSON parsing and output hashing were pinned to CPU 0 with `taskset -c 0`; the available affinity was CPUs 0–19. No theory10 code was imported or run, no Boltzmann solver was run, no numerical fit was attempted, and no file outside `work/compute25/` was edited. This analytic review is CPU-light and stayed within the 600 s / one-CPU task bound. Model: Codex GPT-5; runtime effort setting is not exposed to this worker.

**Bottom line:** the leading asymptotic equations and signs pass independent derivation. Revise the remainder notation and qualify the interpretation of the edge as a failure of the displayed Eq. (19) to realize equal CDM–radiation curvature there, not proof about an unspecified full vacuum-transfer completion.
