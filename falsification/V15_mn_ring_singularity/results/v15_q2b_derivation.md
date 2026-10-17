# V15 Q2b: the leading curvature invariant at a simple pole of the Ernst potential, derived

Pre-registered task: PREREGISTRATION.md, "Amendment Q2b". Code: `code/v15_q2b_check.py`. Numbers:
`results/v15_q2b_check.json`. Environment: Python 3.14.5, sympy 1.14.0, mpmath 1.3.0, at 90 digits; runtime 39 s on
one core at nice 10. Metric: f and e^{2γ} come from `V9_mn_obstruction_check/code/v9_nve.py::mn_functions`, imported
read-only. I re-typed ℰ = e^{2ψ}A₋/A₊ with the same a, b, ψ, and checked numerically that Re ℰ = f(v9). I did not
reuse any Q2 code. Q2's stored K values are used only as comparison data.

Labels: **[derived]** means algebra done here (by hand, or by exact sympy where stated). **[cited]** gives the
reference. **[verified]** means checked numerically, at the precision stated. **[hypothesis]** means a stated input,
and for MN it is checked.

## 0. Result
Let p = (ρ₀ > 0, z₀) be a point where the Ernst potential has a simple pole, and let W = (ρ−ρ₀) + i(z−z₀).
Assume (H1, §3) that ℰ = N/D near p, with N and D real-analytic, D(p) = 0, dD(p) ≠ 0 and N(p) ≠ 0. Then the
Ernst equation forces ℰ = c/W + O(1) (or the mirror image c/W̄), and

    I ≡ ½ 𝒬_ij 𝒬^ij = (C_abcd C^abcd − i C_abcd *C^abcd)/16 = 3 (c/c̄)² / (ρ₀² Λ² W⁶) · (1 + O(|W|)),
    Λ ≡ lim_{W→0} e^{2γ}/(f |W|²) = lim e^{2γ}/Re(c W̄).

In terms of K and |I|:

    K = 16 Re I = 48 cos(6θ_W − 4 arg c) / (ρ₀² Λ² |W|⁶),   |I| = 3/(ρ₀²Λ²|W|⁶) = 3/(8 ρ₀² Λ^{1/2} s³),

where θ_W = arg W and s ≈ √Λ |W|²/2 is the proper distance to the ring.
* **Free Laurent data drop out.** The coefficient depends on the residue c only through its phase. It does not
  depend on any subleading Laurent coefficient. Through the order needed, those are c₀, c₁ and one non-holomorphic
  mode b₃.
* Λ is the integration constant of γ, the regularised value e^{2γ̂(p)} with γ̂ = γ − ½ ln(f|W|²). The local field
  equations cannot fix it; axis regularity fixes it globally. That is the only non-local input.
* **Equatorial rings of reflection-symmetric solutions** (MN, Kerr) have c real (§6), so

      K ≈ 48 cos 6θ_W/(ρ₀²Λ²|W|⁶).

  On y = 0 this gives K·δ⁶ → 48 R_s⁴/(Λ² k⁸ x_s⁶). Along x = x_s it gives K·y⁶ → −48/(Λ² k⁸ x_s⁶ R_s²).
* This is exactly Q2's assumed law, I_Q2 = 3/(ρ₀²Λ²W̄⁶). Q2's sign convention for K* returns the complex
  conjugate.

**PASS (pre-registered):** the closed form reproduces the numerical limit of K·δ⁶ at all four rings. The worst
relative difference is **1.4e−27** (criterion ≤ 1e−6), and the Kerr ring is exact. Details are in §8.

## 1. Conventions and field equations [cited; verified]
Line element and field equations:
* ds² = −f(dt − ω dφ)² + f⁻¹[e^{2γ}(dρ² + dz²) + ρ²dφ²], with ℰ = f + iχ.
* Ernst equation (Ernst, Phys. Rev. 167, 1175 (1968)): f ∇²ℰ = ∇ℰ·∇ℰ, using flat cylindrical ∇ in (ρ, z, φ).
* γ equations (same paper): γ_ρ = (ρ/4f²)(|ℰ_ρ|² − |ℰ_z|²) and γ_z = (ρ/4f²)·2 Re(ℰ_ρ ℰ̄_z).

Complex form, with ∂ = ∂_W = ½(∂_ρ − i∂_z) and ℰ̄_W ≡ ∂_W(conj ℰ) [derived]:

    (E)  f [4ℰ_{WW̄} + (ℰ_W + ℰ_W̄)/ρ] = 4 ℰ_W ℰ_W̄,       (G)  γ_W = ρ ℰ_W ℰ̄_W / (2f²).

(E) and (G) hold for the transcribed MN metric to 1e−84 (at the ring, d = 1e−3) and to 1e−89 elsewhere, and for
Kerr to 1e−90. These are JSON `ernst_rel` and `gamma_eq_rel` **[verified]**. ℰ and ℰ̄ enter (E) and (G)
symmetrically, so the sign convention for χ is irrelevant to everything below.

## 2. An exact identity: I from ℰ alone, with γ only as an overall factor [derived; verified]
**2.1 Setup.** Let ξ = ∂_t, λ = −ξ·ξ = f and F_ab = ∇_aξ_b.
* Self-dual Killing form: 𝓕 = F + iF*.
* Killing identity: ∇_a∇_bξ_c = −R_bcad ξ^d [cited: Wald 1984, App. C, Eq. C.3.6]. In vacuum this gives
  ∇_a𝓕_bc = −𝒞_bcad ξ^d, with 𝒞 = C + iC*.
* Ernst one-form: σ_a ≡ 2ξ^b𝓕_ba = ∇_a(f − iχ′) [derived from ξ^bF_ba = ½∇_aλ and ξ^bF*_ba = −½ω_a]. Call this
  potential ℰ̃; it is ℰ or ℰ̄, and either choice works.

**2.2 The Hessian of ℰ̃.** Differentiate σ_c. With u = ξ/√λ and 𝒬_bc ≡ 𝒞_abdc u^a u^d (the electric plus i
magnetic Weyl part):

    ∇_bσ_c = 2F_b{}^a𝓕_ac + 2λ𝒬_bc.

Then use the self-dual identities [standard; re-derived here]:
* 𝓕_b{}^a𝓕_ac = −¼g_bc𝓕²;
* 𝓕̄_b{}^a𝓕_ac = F_b{}^aF_ac + F*_b{}^aF*_ac, which is the Maxwell stress of 𝓕 (the cross terms cancel by
  F_{ba}F*_c{}^a = ¼g_bc F·F*);
* 𝓕_abξ^b = −σ_a/2, so the electric vector of 𝓕 is v = −σ/(2√λ), with 𝓕² = −4v·v = −σ·σ/λ.

These give, on the quotient space with metric h = g + u⊗u (D is its Levi-Civita connection):

    𝒬_ij = (1/2f) [ D_iD_jℰ − (σ_iσ̄_j + σ̄_iσ_j)/(4f) + h_ij (σ·σ̄ − σ·σ)/(4f) ],   σ = dℰ.

Here the spatial projection of ∇∇ℰ is D D ℰ, because ℰ is ξ-invariant [cited: O'Neill, Michigan Math. J. 13
(1966) 459: the horizontal part of ∇_XY for basic fields is the base connection]. The relation between the 4D
and quotient d'Alembertians is □ℰ = Δ_hℰ + (σ·σ + σ̄·σ)/(4f). With it, **Tr 𝒬 = 0 is equivalent to the Ernst
equation** λ□ℰ = σ·σ [derived; consistency check]. Everything is rational in f, so the formulas continue
analytically to f < 0.

**2.3 Axisymmetric reduction.** The quotient metric is h = H(dρ² + dz²) + (ρ²/f)dφ² with **H ≡ e^{2γ}/f**.
* ℰ is φ-independent, so 𝒬_{Wφ} = 0.
* Conformal 2D Christoffels: Γ^W_WW = ∂_W ln H, and Γ^·_{WW̄} = 0.
* Tracelessness gives 𝒬^φ_φ = −(4/H)𝒬_{WW̄}.

Hence I ≡ ½ Tr 𝒬² is

    I = (4/H²) ( 𝒬_WW 𝒬_W̄W̄ + 3 𝒬_WW̄² ),                                                      (★)
    𝒬_WW  = (1/2f) [ ℰ_WW − (∂_W ln H) ℰ_W − ℰ_W ℰ̄_W/(2f) ]
          = (1/2f) [ ℰ_WW + ℰ_W²/(2f) − ρ ℰ_W² ℰ̄_W / f² ]                       (using (G))
    𝒬_W̄W̄ = the same with W → W̄ (ℰ̄_W̄ = conj ℰ_W)
    𝒬_WW̄ = (1/2f) [ ℰ_WW̄ − ℰ_Wℰ_W̄/(2f) ] = (1/2f) [ ℰ_Wℰ_W̄/(2f) − (ℰ_W+ℰ_W̄)/(4ρ) ]    (using (E))

**Two structural facts follow from (★).**
1. γ enters only as the overall factor H⁻² = f²e^{−4γ}. Its derivatives were removed by (G).
2. Nothing else in (★) depends on ω or on the γ constant.

**Verification [verified]:**

| Check | Points | Result |
|---|---|---|
| (★) against Kerr's closed-form K and 3Ψ₂² | 5 generic points, both prolate sheets | ≤ 1.2e−90 and ≤ 3.2e−90 |
| (★) against Q2's independent direct-Riemann K | C0 and F5-R, y = 0, d = 1e−2…1e−5, both sides (16 points) | ≤ 3.6e−40 (Q2 printed 40 digits) |
| (★) against Q2's negative control | C0, (5, 0.3) | 5e−26 (Q2 printed 25 digits) |
| Raw form (metric's own ∂ ln H and ℰ_WW̄) against substituted form | generic points; near the ring | ≤ 2e−90; 7e−82 |

So the normalisation K = 16 Re I = 8(E²−B²) is right.

## 3. Local hypothesis at the pole [hypothesis; derived and verified for MN]
**(H1)** ℰ = N/D near p, with N, D real-analytic, D(p) = 0, dD(p) ≠ 0 and N(p) ≠ 0.

**Holomorphy of the pole [derived].** Write D = αW + βW̄ + O(2) and ℰ ≈ c/ζ with ζ = αW + βW̄. At order |W|⁻⁴,
(E) gives 4cαβc̄/(ζ̄ζ³) = 0. So αβ = 0: the pole is holomorphic in W, or in W̄. Take W (the W̄ case is the
mirror z → −z, with the same c). This is Q2 §D, re-derived.

**Form of the expansion [derived].** Expanding N/D gives ℰ = Σ_{n≥−1} ℰ_n, with ℰ₋₁ = c/W and ℰ_n homogeneous of
degree n:

    ℰ_n = Σ_{k=0}^{2n+2} a_{n,k} W^{n−k} W̄^k.

The W̄²/W terms come from the quadratic part of D. As a robustness test, I widened the ansatz to k ≤ 2n+6. Every
extra coefficient is forced to zero.

**For MN [verified].** N = e^{2ψ}A₋ and D = A₊ are real-analytic at the rings, since R > 0 there.
* ∂_W̄A₊ = 0 to ≤ 4e−88 at all four rings, so the pole is holomorphic in W at the MN rings and in W̄ at Kerr's.
* A₋(p) = −8.30, −8.34, −56.7 and −2.79 at C0, F5-R, the (0.6, 0.1024) ring and Kerr's ring.
* |∂_WA₊| = 28.08, 0.316, 104.2 and 1.254.

## 4. Perturbative solution of the Ernst equation near the pole [derived; exact sympy]
Scale W → εW and solve (E)·ρ order by order, with ρ = ρ₀ + (W+W̄)/2 exact. W, W̄ are independent symbols, and
conj(a_{n,k}) is an independent symbol.
* **Order ε⁻⁴: identically zero.** Any holomorphic ℰ solves the planar part, since ℰ_W̄ = 0 makes both sides of
  4fℰ_{WW̄} = 4ℰ_Wℰ_W̄ vanish.
* **Order ε⁻³ fixes ℰ₀ completely, up to a constant c₀.** It is forced by the 1/ρ term:

      ℰ₀ = c₀ − (c/4ρ₀) W̄/W − (c²/8c̄ρ₀) W̄²/W².

* **Order ε⁻² fixes ℰ₁ up to two free modes.** These are c₁W (holomorphic) and one non-holomorphic homogeneous
  mode b₃·[(c̄W + cW̄)³/(c³W²) − (c̄/c)³W]:

      ℰ₁ = c₁W + b₃[(c̄W+cW̄)³/(c³W²) − (c̄/c)³W] − (Re c₀/2ρ₀)(W̄ + cW̄²/(c̄W))
           + (1/32ρ₀²)[(3c − 9c̄)W̄ + (c² − 6cc̄)W̄²/(c̄W) + c³W̄⁴/(2c̄²W³)].

  The wide ansatz finds the same free set {c₀, c₁, b₃}.
* **Consistency check: H from (G) [derived].** ∂_W ln H = 2γ_W − f_W/f = 1/W + O(1). So
  **H = Λ|W|²(1 + O(|W|))**, with Λ an integration constant.
  * The leading part can be checked by hand. γ_W = c̄/(2(c̄W + cW̄)) and f_W/f = c̄/(c̄W + cW̄) − 1/W.
  * So e^{2γ} ≈ Λ Re(cW̄) vanishes linearly at the ring. This is Q2's "e^{2γ}/f ∝ ℬ", now derived from the field
    equations rather than read off the MN formulas.

## 5. Expansion of (★): which orders vanish and why [derived; exact sympy, with the key step by hand]
Write D₀ = c̄W + cW̄ = 2|W|² f₋₁. The leading terms (exact rational functions; JSON `A_series_*`) are:

| component | naive order | first non-zero order | leading term | origin |
|---|---|---|---|---|
| 𝒬_WW  | ε⁻² | ε⁻² | 3cW̄/(W²D₀) = 3c/(2f₋₁W³) | ℰ₋₁ alone: the pp-wave (Q2 §E) |
| 𝒬_WW̄ | ε⁻² | ε⁻¹ | cW̄/(4c̄ρ₀W²) | 1/ρ term, first order |
| 𝒬_W̄W̄ | ε⁻² | **ε⁰** | 3cW̄D₀/(16c̄²ρ₀²W²) | 1/ρ term, second order |

**Why the naive leading order vanishes.** ℰ₋₁ = c/W is holomorphic. So at the naive order only 𝒬_WW survives.
Since h^{WW} = 0, that 𝒬 is null and Tr𝒬² = 0: type N, the vacuum pp-wave of Q2 §E. Naive counting would give
X ~ ε⁻⁴, so I ~ H⁻²ε⁻⁴ ~ |W|⁻⁸ (Q2's "degree −8").

**What breaks the degeneracy.** Only the cylindrical 1/ρ term of (E) puts W̄-dependence into ℰ:
* 𝒬_WW̄ first appears at relative order |W|/ρ₀.
* 𝒬_W̄W̄ has an order-ε⁻¹ term that cancels identically. By hand: with ℰ_W̄ = (κ₁W + 2κ₂W̄)/W², κ₁ = −c/4ρ₀ and
  κ₂ = −c²/8c̄ρ₀,

      2f𝒬_W̄W̄ ⊃ 2κ₂/W² + 4ρ₀c̄(κ₁W + 2κ₂W̄)²/(W²D₀²) = −c²/(4ρ₀c̄W²) + c²/(4ρ₀c̄W²) = 0.

  So 𝒬_W̄W̄ first appears at relative order (|W|/ρ₀)².

Hence, with X ≡ 𝒬_WW𝒬_W̄W̄ + 3𝒬_WW̄²:

    X = (9/16 + 3/16) c² W̄² / (c̄² ρ₀² W⁴) + O(ε⁻¹) = 3c²W̄²/(4c̄²ρ₀²W⁴) + O(ε⁻¹),
    I = 4X/H² = 3 (c/c̄)² / (ρ₀² Λ² W⁶) · (1 + O(|W|)).                                     ∎

* sympy confirms that X has no terms below ε⁻², and that X₋₂ equals the expression above exactly
  (`X_leading_equals...`: true).
* **The free symbols of X₋₂ are only {W, W̄, c, c̄, ρ₀}.** So c₀, c₁ and b₃ (and |c|) all drop out, in both the
  standard and the widened ansatz.
* c₀ and b₃ do appear in the O(ε⁻¹) parts of 𝒬_WW and 𝒬_WW̄. They enter I only at relative order |W|.
* **Which Laurent data enter the first non-vanishing term:** arg c, ρ₀ (through the 1/ρ forcing), and Λ (through
  H⁻²). This is precisely the "universal data" that Q2 assumed.

## 6. Reflection symmetry, real c, and the Kerr check
**Real c [derived].** Suppose the solution is symmetric under z → −z (ω even, χ odd, so ℰ(ρ,−z) = conj ℰ(ρ,z)),
and the ring sits at z = 0. Then c/W̄ = conj(c/W) = c̄/W̄, so c = c̄ and the phase factor is 1.
* Verified at all four rings: arg c = π, π, π and 0 (c = −0.02306, −2.050, −0.02230 and 2M²/a).
* **Phase test off the equator [verified; not part of PASS].** MN with χ 0.9, q = −1 has a B = A = 0 point at
  (x, y) = (1.72495, 0.80268), from Q1's q < 0 family. There arg c = −1.0233.
  * I·W⁶ equals the complex coefficient 3(c/c̄)²/(ρ₀²Λ²) = −0.60547 + 0.84986i to 1.7e−11 on 7 rays.
  * K|W|⁶ = 48 cos(6θ − 4 arg c)/(ρ₀²Λ²) holds there to 1.2e−11 (relative to the scale 48/(ρ₀²Λ²)).
  * Without the factor (c/c̄)² this test would fail by O(1).

**Kerr (β = 0), exact [derived; sympy].**
* Local data: near r = 0, cos θ = 0, W = −(M/a)(r + ia cos θ) + O(2). The Ernst potential
  ℰ = 1 − 2M/(r − ia cos θ) has its pole in W̄, with residue c = 2M²/a (real). Also ρ₀ = a and
  H = Σ/((r−M)² − k²cos²θ), so Λ = a²/M⁴.
* Hence I_lead = 3/(ρ₀²Λ²W̄⁶) = 3M⁸/(a⁶W̄⁶) = **3M²/(r − ia cos θ)⁶ = 3Ψ₂² exactly** (sympy ratio = 1).
* On the equator: K·δ⁶ = 48M²/k⁶. This is exact for all δ, since there K = 48M²/r⁶ with r = kδ.
* Numerically (M = 1, a = 0.9):

| Quantity | Closed form | Numerical | Agreement |
|---|---|---|---|
| Λ | 0.81 | 0.81 | — |
| Coefficient | 5.6450292695 = 3M⁸/a⁶ | — | — |
| K·δ⁶ | 6998.104679982505 = 48/k⁶ | from (★) | 3.6e−75 |

## 7. Λ in closed form (MN, equatorial ring) [derived; Q2 §D formula re-derived and verified]
**Generic statement.** Λ = lim e^{2γ}/(f|W|²) = e^{2γ̂(p)}, the regularised γ at the ring. It is fixed by the γ
integration constant, which is set by axis regularity, a global condition.

**For MN.** v9's formulas give exactly H = e^{2γ}/f = e^{Γ−2ψ}|A₊|²/((x²−y²)(1−α²)²): the 𝒜 cancels, and
ℬ = |A₊|². With A₊ = λW + O(2), at y = 0:

    Λ = e^{Γ₀−2ψ₀} |λ|² / (x_s² (1−α²)²),   |λ| = R_s |∂_xA₊| /(k x_s),   ∂_xA₊ = −2/(x_s−1) + 6β(x_s+1)/R_s⁵,
    Γ₀ = −3β²/(4R_s⁶) + 4βx_s/R_s − 2βx_s/R_s³ − 4β,   ψ₀ = −β/(2R_s³),   ρ₀ = kR_s.

**Through the Ernst residue.** With c = e^{2ψ₀}A₋(p)/λ and A₋(p) = −2(1+a_p²):

    Λ = e^{Γ₀+2ψ₀} A₋(p)² / (|c|² x_s² (1−α²)²).

So the coefficient of I is 3/(ρ₀²Λ²) = 3|c|⁴x_s⁴(1−α²)⁴/(k²R_s² e^{2Γ₀+4ψ₀}A₋(p)⁴).

**Checks.** Both forms equal the numerical limit H/|W|² at d = 1e−30 along 3 directions, to all 25–30 printed
digits, at every ring. Equatorial limits (W = (kx_s/R_s)δ on y = 0; W = ikx_s y on x = x_s):

    K·δ⁶ → 48 R_s⁶/(ρ₀²Λ²k⁶x_s⁶) = 48 R_s⁴/(Λ²k⁸x_s⁶),     K·y⁶ → −48/(ρ₀²Λ²k⁶x_s⁶).

## 8. Pre-registered check [verified]: PASS
**Numerical limit.** K is computed from the exact identity (★), using exact sympy derivatives of the exact ℰ and H
at 90 digits. (★) is independently tied to Q2's direct-Riemann code (§2). For each ring:
* average K(x_s ± δ)δ⁶ over the two sides at δ = 1e−5, 1e−6, 1e−7 and 1e−8;
* apply two Richardson steps in δ².

The closed form is §0 and §7, using the Laurent data from the exact functions.

| ring | ρ₀ | Λ | 3/(ρ₀²Λ²) | K·δ⁶ closed form | rel. error vs limit | raw δ = 1e−8 | K·y⁶ rel. error |
|---|---|---|---|---|---|---|---|
| C0 (χ 0.9, q 0.95) | 0.7194477399 | 6.9676524044 | 0.11938500880 | 109.00131048361711 | **9.7e−35** | 3.5e−13 | 1.1e−34 |
| F5-R (89, 80, q 19/20) | 64.009758366 | 8.6906047338e−4 | 969.45774291 | 1.7097608273991279e−6 | **1.0e−34** | 3.6e−13 | 1.1e−34 |
| (χ, q) = (0.6, 0.1024) | 0.3176862940 | 0.5603798487 | 94.658426034 | 14.601860639839696 | **1.4e−27** | 1.4e−10 | 1.1e−31 |
| Kerr ring (M 1, a 0.9) | 0.9 | 0.81 | 5.6450292695 | 6998.1046799825047 | **3.6e−75** | 3.5e−75 | 2.2e−59 |

**Comparison with Q2's stored finite-d numbers.**
* K·δ⁶ at d = 1e−5: 3.5e−7 (C0) and 3.6e−7 (F5-R).
* K·y⁶: 3.2e−7 and 3.2e−7.
* Third ring, at d = 1e−6: 1.4e−6.
* These residuals are the O(d²) truncation; the δ-sequence above shows them falling by 100 per decade.
* Q2's predicted values agree with the closed form here to 1.0e−15 and 1.2e−15, which is Q2's printed precision.

**Ingredient checks.** The individual leading terms of 𝒬_WW, 𝒬_WW̄ and 𝒬_W̄W̄ from §5 match the actual MN and Kerr
components at d = 1e−9, to O(d) residuals (1e−9 to 4e−6, JSON `ingredient_rel_err_d1e-9`).

## 9. Status of each step
| step | status |
|---|---|
| Ernst equation and γ equations (§1) | cited (Ernst 1968); verified for this metric to ≤ 1e−84 |
| Killing identity ∇∇ξ = −Rξ | cited (Wald 1984, Eq. C.3.6) |
| self-dual 2-form algebra, Maxwell-stress form of 𝓕̄𝓕 | derived (standard) |
| projected Hessian = quotient Hessian | cited (O'Neill 1966); verified through the identity checks |
| identity (★) and its substituted form | derived; verified (Kerr 1e−90; MN vs Q2 direct Riemann 4e−40) |
| H1 (ℰ = N/D, simple zero, N ≠ 0) | hypothesis; for MN verified (Q2 §B–D; here ∂_W̄A₊ = 0 to 4e−88) |
| holomorphy of the pole | derived (order |W|⁻⁴ of (E)) |
| perturbative ℰ through degree 1, free data {c₀, c₁, b₃} | derived (exact sympy; widened ansatz agrees) |
| H = Λ|W|²(1 + O(|W|)) | derived from (G); Λ = integration constant (non-local) |
| cancellations, and X₋₂ = 3c²W̄²/(4c̄²ρ₀²W⁴) | derived (exact sympy; the ε⁻¹ cancellation also by hand) |
| I = 3(c/c̄)²/(ρ₀²Λ²W⁶) and K = 48 cos(6θ − 4 arg c)/(ρ₀²Λ²|W|⁶) | derived |
| c real for reflection-symmetric equatorial rings | derived; verified (arg c ∈ {0, π}) |
| Λ closed form for MN | derived (from the transcribed metric); verified (25–30 digits) |
| Kerr reduction to 3Ψ₂² | derived (sympy, exact) |
| 4-ring coefficient match | verified, ≤ 1.4e−27 |

**Nothing is assumed beyond H1, which holds at the MN rings.** Q2's assumption, that the first correction to the
pp-wave is universal, is now a theorem under H1: the correction is driven by the 1/ρ term, and the free
subleading Laurent data cancel at the leading order.

## 10. Caveats
* **Scope.** The result is local and asymptotic: leading order, with relative corrections O(|W|). It needs ρ₀ > 0.
  An on-axis pole would need separate treatment.
* **Rigour of the series step.** The series is formal, carried to the order needed. The real-analyticity of N and
  D (H1) guarantees that the expansion exists. Its order-by-order coefficients are exact.
* **Λ.** Λ depends on the γ normalisation. The "−2 inside" reading of GLM Eq. 3h multiplies e^{2γ}, and so Λ, by
  e^{−8β}. Then (★) multiplies I and K by exactly e^{16β}, which is V's finding. The d⁻⁶ law and the cos 6θ_W
  structure are unchanged.
* **Proper distance.** The proper-distance form |I| = 3/(8ρ₀²Λ^{1/2}s³) uses s = √Λ|W|²/2 along rays. At C0 this
  gives s = 0.3429 d² on y = 0, matching V's 0.343 d².
