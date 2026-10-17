# V15 Q2: the Manko–Novikov equatorial B = A = 0 point is a curvature singularity, K ∝ d⁻⁶

Code: `code/v15_q2_curvature.py`. Numbers: `results/v15_q2_curvature.json` (Python 3.14.5, sympy 1.14.0,
mpmath 1.3.0, 100 digits, with a 150-digit repeat). Metric: GLM 2008 Eqs 2–7 / Brink 2008 App. C, re-typed
from the `V9_mn_obstruction_check/code/v9_nve.py` transcription. I did not read `v13_kretschmann.py` or
tabula's script. I opened the V13 K table only after my own numbers were on disk (the `--compare-v13` step).

Notation: a, b are the GLM functions (`A_`, `B_` in v9_nve.py), and 𝒜, ℬ, 𝒞 are the GLM A, B, C.
y = 0 is the equator. x_s is the root, δ = x − x_s, and d is the coordinate distance.
R = √(x²+y²−1), and R_s = √(x_s²−1) at the ring.
Each claim carries a label: **[derived]** (algebra done here), **[verified]** (checked numerically to the
precision stated), or **[assumed]** (used but not proven).

## Summary
* At B = A = 0 the Ernst potential ℰ has a simple pole. To leading order the pole is holomorphic in the Weyl
  variable W = (ρ−ρ₀) + i z, so B ∝ |W|² vanishes quadratically and isotropically in (ρ, z).
* The complex Weyl invariant I = (K − iK*)/16 behaves like 3/(ρ₀²Λ² W̄⁶). So |I| ∝ B⁻³, and
  K = 16 Re I ≈ 48 cos 6θ_W /(ρ₀²Λ²|W|⁶).
* **Leading exponent: K ∝ d⁻⁶ in coordinate distance**, which is ∝ s⁻³ in proper distance s ∝ d².
  On y = 0 the angular factor is exactly 1. Λ is given in closed form below.
* The exponent and the coefficient are **partially derived**. The structure up to the pole and the pp-wave
  leading order is derived. The final step, that the first correction to the pp-wave is the universal
  Kerr-like one, is **assumed**, and it is **verified** numerically at four rings: C0, F5-R, the
  (χ, q) = (0.6, 0.1024) root, and the Kerr ring itself.

## A. Curvature route (new code) [derived; verified]
The metric is a 2+2 orthogonally transitive split: ds² = G_AB dyᴬdyᴮ + h_ij dxⁱdxʲ, with yᴬ = (t, φ) and xⁱ = (x, y).

    G = [[−f, fω], [fω, ρ²/f − fω²]],   h = diag(H, H(x²−1)/(1−y²)),   H = e^{2γ}k²(x²−y²)/(f(x²−1)).

Since ∂_A = 0, the only Christoffels are ʰΓⁱ_jk, Γⁱ_AB = −½∂ⁱG_AB and Γᴬ_iB = ½(G⁻¹∂_iG)ᴬ_B. They give
the following closed-form blocks. Every component with an odd number of fibre indices vanishes.

    R_ijkl = K_h (h_ik h_jl − h_il h_jk)                       (2D Gauss curvature)
    R_ABCD = −¼ hⁱʲ(∂_iG_AC ∂_jG_BD − ∂_iG_AD ∂_jG_BC)         (Gauss equation of the orbits)
    R_ABij = −¼ (∂_iG G⁻¹∂_jG − ∂_jG G⁻¹∂_iG)_AB              (normal curvature)
    R_iAjB = −½ (D_iD_jG)_AB + ¼ (∂_jG G⁻¹∂_iG)_AB

The code builds exact sympy derivatives up to second order, lambdifies them to mpmath, assembles the full
R_abcd, and contracts it numerically. It computes K, K* = *R_abcd R^abcd, I, J = (CCC − i CC*C)/96, and the
speciality index S = 27J²/I³. It also computes the Ricci tensor and the first-Bianchi residual.

Controls:
* **Kerr, positive control.** Kerr (M = 1, a = 0.9) at 3 points matches the closed-form K to 98.6–100 digits.
  |K*| matches 96M²a r c(3r²−a²c²)(r²−3a²c²)/Σ⁶ to 99 digits and |Ψ₂| = |M/(r−iac)³| to about 100 digits.
  S − 1 is at most 3e−100 and R_abR^ab/K at most 1e−199. **[verified]**
  * My K* sign convention returns Ψ₂ complex-conjugated, so I carries the conjugate orientation. This does not
    affect K, |I| or any modulus.
* **Negative control.** C0 at (5, 0.3) gives K = 0.0835792611180587…. Under ±1e−6 shifts K moves by
  1.07e−6 (x) and 2.47e−6 (y), linearly; the second difference is 1.4e−12. Ricci²/K is about 3e−199. **[verified]**

## B. Exact Ernst structure [derived algebraically; verified to 1e−99]
Let 𝒲 = (1−ia)(1−ib), A_± = x(1+ab) + iy(b−a) ± 𝒲. Then:
* ℬ = |A₊|². Its real part is (x+1)+(x−1)ab and its imaginary part is −[(1+y)a+(1−y)b].
* 𝒜 = Re(A₋Ā₊) = |x(1+ab)+iy(b−a)|² − |𝒲|². This uses |𝒲|² = (1+ab)² + (a−b)².
* ℰ = e^{2ψ}A₋/A₊ has Re ℰ = e^{2ψ}𝒜/ℬ = f.
* ℰ satisfies the prolate Ernst equation Re ℰ [∂_x((x²−1)ℰ_x) + ∂_y((1−y²)ℰ_y)] = (x²−1)ℰ_x² + (1−y²)ℰ_y².
  Residuals are below 2e−99 at 6 points. This was checked numerically, not proven symbolically.

From the metric formulas, exactly:
* e^{2γ}/f = exp(Γ − 2ψ) ℬ / ((x²−y²)(1−α²)²), where Γ is the γ exponent. The 𝒜 cancels, so **the base metric
  is ∝ ℬ**.
* det G = −ρ², which stays finite.

## C. The equator [derived]
**Reflection symmetry.** u → −u and P_l → (−1)ˡP_l give a(x,−y) = −b(x,y). So on y = 0:
* a = −b, ab is even and a+b is odd in y;
* A₊ is real and 𝒲 = 1+a² is real;
* ∂_xA₊ is real and ∂_yA₊ is purely imaginary.

**Consequences on y = 0.**
* ℬ = A₊², and 𝒜 = A₊(A₊ − 2(1+a²)).
* **B = 0 ⇒ A = 0**, and near the root 𝒜 ≈ −2(1+a²)A₊ is linear and changes sign. So the ergosurface passes
  through the ring.
* f ≈ −2(1+a²)e^{2ψ}/A₊ ∝ 1/d. f is a pole, not a zero: f = −45.40 at d = 1e−3 and −452.6 at d = 1e−4 (C0).
* A₋ = −2(1+a²) ≠ 0 at the ring (−8.3019 at C0), so **ℰ has a simple pole**.

**Root condition.** Write E = exp 2β(1 − x/R + x/(2R³)), so ab = −α²E². Then A₊ = x+1 − (x−1)α²E² = 0 is the
Q1 equation α²E² = (x+1)/(x−1).

**Linear coefficient.** Using (ln E²)' = −6β/R⁵, at the root

    c₁ ≡ ∂_xA₊ = −2/(x_s−1) + 6β(x_s+1)/R_s⁵.

This is 14.310752812093860163 at C0 and 14.449887135042295765 at F5-R, and equals the numerical derivative to
25 digits. c₁ vanishes only if 3β = R_s³, which does not hold at C0 or F5-R. **[derived; verified]**

## D. The pole is holomorphic in Weyl coordinates [derived], so B ∝ d²
**Setup.** Write ℰ = c/ζ + O(1) with ζ = λ₁(ρ−ρ₀) + λ₂z. Near (ρ₀ > 0, 0), substitute into
(Re ℰ)(ℰ_ρρ + ℰ_ρ/ρ + ℰ_zz) = ℰ_ρ² + ℰ_z².
* The order-ζ⁻⁴ part, LHS − RHS, equals (λ₁²+λ₂²)|c|²/(ζ̄ζ³).
* No other term reaches this order (ℰ_ρ/ρ is O(ζ⁻²)).
* So λ₂ = ±iλ₁. The pole is in W_ε = (ρ−ρ₀) + iεz.

**Translation to (x, y).** At y = 0, ∂ρ/∂x = kx/R and ∂z/∂y = kx. So the derivation predicts |∂_yA₊/∂_xA₊| = R_s.
* Numerically: 1.6505263122370314366611 against R_s = 1.6505263122370314366611 (C0), and
  1.641275855529671525 against 1.641275855529671525 (F5-R). **[verified]**
* ε = +1 at the MN rings. At the Kerr ring the orientation is the conjugate one.

**Hence:**
* ℬ ≈ c₁²(δ² + R_s²y²) = |λ|²|W|², with |λ| = c₁R_s/(kx_s). The hint holds: B is a sum of two squares, each
  vanishing linearly, and it is positive definite. So B ∝ d² along every direction.
* The base metric is h ≈ Λ|W|²|dW|² = Λ|du|² with u = W²/2. So **proper distance s ≈ √Λ|W|²/2 ∝ d²**; the
  numerical log s / log d slope is 2.005 on y = 0. **[verified]**
* Λ in closed form [derived from the formulas in B and the linear A₊], at y = 0, with Γ₀ the γ exponent there:

      Λ = exp(Γ₀ − 2ψ₀) c₁² R_s² / (k² x_s⁴ (1−α²)²),
      Γ₀ = −3β²/(4R_s⁶) + 4βx_s/R_s − 2βx_s/R_s³ − 4β,   ψ₀ = −β/(2R_s³).

  This gives Λ = 6.9676524043966923 at C0 and 8.6906047338e−4 at F5-R. It matches e^{2γ}/(f|W|²) at
  d = 1e−9 along three directions to about 1e−8, which is the O(d) residual.

## E. Leading order is a vacuum pp-wave, so naive power counting fails [derived; verified]
**The leading-order model.** Take f + iχ = c/W with c real (rotate W if needed), so f = cξ/|W|² and
χ = −cη/|W|², where W = ξ + iη.
* The twist equations ∂_ξω = ρf⁻²∂_ηχ and ∂_ηω = −ρf⁻²∂_ξχ integrate at ρ ≈ ρ₀ to ω = ω₀ − (ρ₀/c)|W|²/ξ = ω₀ − ρ₀/f.
* With T = t − ω₀φ, the fibre metric collapses to −f dT² − 2ρ₀ dT dφ. So to leading order

      ds² ≈ −Re(c/W) dT² − 2ρ₀ dT dφ + Λ|du|²,   u = W²/2.

* ∂_φ is null and parallel, and Re(c/W) = Re(c(2u)^{−1/2}) is harmonic in u. So this is a **vacuum pp-wave**.

**Why it matters.** All its scalar invariants vanish. The code confirms this for the model: K ≈ 1e−174 to 1e−201,
while the mixed Riemann components reach 1.8e15. Naive counting (G ∝ 1/B, h ∝ B, G⁻¹ ∝ 1/B) would give
degree −8 in |W|. The true leading invariant comes instead from the first correction to the pp-wave. This is why
the uncancelled power has to be fixed beyond leading order.

## F. The exponent [derived up to one assumption; verified]
**Kerr (exact).** ℰ = (ξ−1)/(ξ+1) with ξ+1 = (r − ia cosθ)/M, so A₊ ∝ r − iac. Then
I = 3Ψ₂² = 3M²/(r−iac)⁶ and K = 48M² Re((r+iac)⁶)/Σ⁶, with Σ ∝ B. In local data, ρ₀ = a, Λ = a²/M⁴ and
W ∝ (M/a)(r ± iac). So **I = 3/(ρ₀²Λ² W̄⁶) holds exactly** for Kerr, with coefficient 3M⁸/a⁶.

**Assumption.** For MN, I take the first correction to the pp-wave to be fixed by the universal data. These are:
* ρ = ρ₀ + Re W, which is exact;
* the pole orientation ε;
* Λ.

Non-universal subleading Ernst terms are taken not to enter at that order. Under this assumption

    I ≈ 3/(ρ₀²Λ² W̄_ε⁶),   |I|·B³ → 3|λ|⁶/(ρ₀²Λ²) ≠ 0,   K ≈ 48 cos(6θ_W)/(ρ₀²Λ²|W|⁶),

with θ_W = arg W_ε.

**Exponent.** **K ∝ d⁻⁶ in every direction except the 12 nodal directions θ_W = (2n+1)π/12.** Along those,
the leading K cancels but |I| ∝ d⁻⁶ still holds. On y = 0, θ_W = 0 or π, so **K·δ⁶ → 48R_s⁴/(Λ²k⁸x_s⁶) > 0
on both sides**. Along x = x_s, θ_W = π/2, so **K·y⁶ → −48R_s⁴/(Λ²k⁸x_s⁶R_s⁶) < 0**.

**This is the K-H analogue.** The Weyl scalar is ∝ A₊⁻³, so ∝ B^{−3/2} in modulus, and |I| ∝ B⁻³. I did not read
the Kodama–Hikida paper; the "TS Weyl ∝ 1/B³" statement is taken from the brief and is not checked here.

**Numerical checks of the assumption's consequences** (C0 and F5-R at d = 1e−5; the third ring at d = 1e−6):

| check | C0 | F5-R | (0.6, 0.1024) ring |
|---|---|---|---|
| law 3/(ρ₀²Λ²) | 0.1193850088 | 969.45774 | 94.65843 |
| I·W̄⁶ on y = 0, mean of ± sides | 0.119385 | 969.457 | 94.658 |
| I·W̄⁶ along y | 0.1193850 − 1.9e−4 i | 969.4575 − 1.5 i | 94.65842 − 0.10 i |
| K·δ⁶ predicted / numerical (y = 0) | 109.001310 / 109.001272 | 1.7097608e−6 / 1.7097602e−6 | — |
| K·y⁶ predicted / numerical | −5.391344 / −5.391342 | −8.746733e−8 / −8.746730e−8 | — |
| S = 27J²/I³ | 1 − 9e−7 | 1 − 9e−7 | 1 − 6e−6 |

* The residuals scale as d (I·W̄⁶) or d² (the symmetric mean). With d = 1e−4 the K·δ⁶ residual is 100× larger.
* I·W⁶, the wrong orientation, is not direction-independent: on the diagonal it gives 0.1156 + 0.030 i.
* S → 1 means the ring is **asymptotically Petrov type D**, as in Kerr.
* The Kerr ring at x = −M/k gives I·W⁶ → 5.6450293 = 3M⁸/a⁶ and Λ → 0.81 = a²/M⁴.

## G. Caveats
* The pp-wave model and the Kerr exactness are derived. The universality of the first correction is
  **assumed**, and so far supported only by 4 rings. A full proof needs a perturbative solution of the Ernst
  equation plus the γ quadrature to first order beyond the pp-wave.
* The coordinate exponent −6 is the pre-registered quantity. The invariant statement is K ∝ s⁻³ in proper
  distance, or |I| ∝ B⁻³.
* The e^{2γ} normalisation does not show up in Ricci. Moving GLM's "−2" outside the l-sum multiplies e^{2γ} by
  e^{8β} (about 7e39 at C0). That metric is still Ricci-flat, but the axis is conical and K changes by e^{−16β}.
  The transcription used here has e^{2γ} = 1 to 20 digits on both axis half-lines and at x = 1e12. So Ricci
  flatness alone cannot settle the verifier's test (iii); axis regularity does.
* This code and V13 use the same V9 transcription. Their agreement checks the curvature code, not the
  transcription. The transcription is checked by Ricci = 0 and the Ernst equation.
