# V15 V: adversarial verifier report on the Manko–Novikov equatorial ring

Date: 2026-10-17. Verifier: Opus agent, working independently of the Q1–Q3 code (which I did not open).
Env: Python 3.14.5, sympy 1.14.0, mpmath 1.3.0 (conjecture_machine venv), single core per job, ≤ 2 jobs.

**Verdict: the ring stands.** I tried to break the claim on seven fronts and none of them gave way:
- My own transcription matches the fleet's v9 transcription to 77–79 digits.
- The metric is vacuum to working precision.
- The fleet's γ is the axis-regular reading.
- The multipoles are exactly GLM's.
- A curvature invariant that cannot vanish on any approach direction diverges as d⁻⁶ on every direction.
- The ring sits at finite proper distance outside the horizon rod, and timelike and null geodesics from far away reach it at finite proper time or affine parameter.

## Code and results

| Test | Code | Results |
|---|---|---|
| Transcription | `code/v15_V_mn.py` | (module) |
| Curvature engine | `code/v15_V_curv.py` | (module) |
| 2 | `code/v15_V_t2_compare.py` | `results/v15_V_t2_compare.{json,log}` |
| 3, 4 | `code/v15_V_t3_vacuum.py` | `results/v15_V_t3_vacuum.json`, `_part1.log`, `_part2.log` |
| 5 | `code/v15_V_t5_multipoles.py` | `results/v15_V_t5_multipoles.{json,log}` |
| 4 (curvature), 6, 7 | `code/v15_V_t6_singularity.py` | `results/v15_V_t6_singularity.{json,log}` |

## Sources read

I read all three papers from the arXiv LaTeX e-prints, which I downloaded and read before opening v9:
- **GLM 2008,** arXiv:0708.0628, e-print dated 2008-01-11:
  - Eq. 2 is the prolate line element.
  - Eqs. 3a–3l are the functions f, ω, e^{2γ}, A, B, C, ψ, γ′, a, b, R and P_n.
  - Eq. 4 gives the moments, Eq. 5 the parameters, and Eqs. 6–7 the Weyl coordinates.
- **Brink 2008,** arXiv:0807.1179, App. C ("Manko-Novikov Metric"). It gives the Ernst potential ℰ = e^{2ψ̃}A₋/A₊, the metric functions, and a, b, ψ̃, γ̃.
- **LGAC 2010,** arXiv:1003.3120, Eqs. (ffunc)–(fLeg) and (freepar).

I could not read the original MN 1992 paper (CQG 9, 2477; paywalled, not circumvented). I also did not check the journal versions, only the e-prints.

### Paper-vs-paper differences, and how the tests settled each one

| # | Item | GLM | LGAC | Brink App. C | Settled by |
|---|---|---|---|---|---|
| D1 | −2 in γ′ (GLM 3h) | Literal typesetting: −2 inside Σ_{ℓ=0}^{2}, so −6β in total | −2 outside. LGAC say they moved it "so as to avoid any confusion" | γ̃ ∋ −½(ln(a/−α)+ln(b/α)). This expands *exactly* to the −2-outside form (checked by hand) | Axis regularity (test 4): **outside is regular** |
| D2 | β² coefficient in γ′ | 9β²/(6R⁶) = 3/2 | 3β²/(2R⁶) = 3/2 | α₂²(3/2)² = **9/4** | Vacuum (test 3): **9/4 is not vacuum** (Ricci/√K up to 0.58). 3/2 is also the standard Weyl quadrupole value (n+1)/2 at n = 2 |
| D3 | P₂ term in b (GLM 3j) | P₂/R³ | P₂/R³ | **P₂/r²** | Vacuum: **Brink's r² is not vacuum** (Ricci/√K up to 5.7). It is a typo |
| D4 | Line element | k²f⁻¹e^{2γ}(x²−y²)(dx²/(x²−1)+dy²/(1−y²)) | Weyl form | k²e^{−2ψ}[e^{2γ}(dρ²+dz²)+R²dφ²] with R² = k²(x²−1)(1−y²), so k⁴ in g_φφ unless k = 1 (MN's convention, per Brink's comment). His "ρ, z" are x = cosh ρ, y = cos z, so his e^{2γ} carries an extra (x²−y²) | Notational; equivalent to GLM for k = 1. Not tested further |

The a/b exponents are the same in all three papers once D3 is fixed:
- a = −α exp[−2β(−1 + Σ(x−y)P_ℓ/R^{ℓ+1})]
- b = α exp[2β(1 + Σ(−1)^{3−ℓ}(x+y)P_ℓ/R^{ℓ+1})]

Two side observations, neither an error:
- GLM (Sec. II.B) and LGAC both describe a *different* singularity: a "circular line singularity" at x = 1, y = 0, the equator of the horizon rod, where R = 0. Neither paper, in the sections I read, mentions an equatorial singularity at x_s ≈ 1.93.
- GLM's Sec. V says the ergodic region "touches the CTC region at a single point". The ring is the outermost point of the CTC region (Q3), so this may be the same point. That is my reading, unverified.

## Results by test

### 1. Independent transcription: PASS (verified)
`v15_V_mn.build()` follows GLM 3a–3l line by line, with switches for D1, D2 and D3. P_n come from the Rodrigues formula, as GLM 3l states it. Parameters follow GLM Eq. 5.

### 2. Comparison with v9_nve.mn_functions: PASS (verified)
- **Parameters:** k, α and β (from q) are identical (difference 0, or 1.8e−279).
- **Numerics:** f, ω and e^{2γ} at 12 random points per set, at 80 digits. Worst relative difference: **5.0e−79 (C0) and 6.6e−77 (F5-R)**, with γ in the outside reading.
- **What v9 does with D1:** v9 types GLM's "−2 inside" literally and then adds +8β to the exponent. Algebraically that is the regular outside reading.
- **Symbolic:** the γ-exponent residual factors as (x√((x²+y²−1)/x²) − √(x²+y²−1))·(…). That is only v9's spelling of R, which is identical for x > 0.
- Against my GLM-literal (inside) reading, e^{2γ} differs by e^{−8β}, as expected.
- Minor: v9 applies `nsimplify` to its inputs. This is exact for the rational V13 POINTS, but would silently rationalise a float β.

### 3. Vacuum: PASS (verified)
**Full 4D Ricci.** Christoffel → Riemann in (t, x, y, φ), with exact sympy jets. Measure: max|R^a_b|/√|K|, at 17 points.
- Points: 10 generic C0 points, 4 C0 points at d = 1e−3 from the ring (+x, −x, +y, diagonal), and 3 F5-R points (2 of them at d = 1e−3).

| Points | At 60 digits | At 100 digits |
|---|---|---|
| Generic | ≤ 9.7e−58 | ≤ 2.3e−97 |
| Near the ring | ≤ 9.7e−46 | ≤ 8.6e−87 |

The residual falls with precision, so it is round-off and the metric is vacuum.

**Variants (C0, 5 points, 60 digits):**

| Variant | max Ricci/√K | Vacuum? |
|---|---|---|
| GLM-literal −2 inside | ≤ 1e−47 | yes, as predicted |
| Brink 9/4 | 0.58 | no |
| Brink P₂/r² | 5.7 | no |

So the test has teeth.

**Ernst route.** Brink's ℰ = e^{2ψ}A₋/A₊ at all 17 points:
- It satisfies the Ernst equation to ≤ 7e−54 relative.
- Re ℰ = GLM f to ≤ 2e−57.
- The twist relations ∂_ρΦ = ∓(f²/ρ)∂_zω and ∂_zΦ = ±(f²/ρ)∂_ρω hold with ratio exactly −1, a sign convention. On y = 0 both sides vanish.
- So GLM's ω and Brink's potential, two independent printed formulas, describe the same solution.

**Kerr control (β = 0):** K matches 48M²(r⁶−15a²r⁴c²+15a⁴r²c⁴−a⁶c⁶)/Σ⁶ to 1e−59, and the Pontryagin invariant matches its closed form.

**Ring location and V13 values:**
- My own root of x+1+(x−1)ab = 0 gives x_s = 1.929828258521150008… (C0) and 1.921922587916759772… (F5-R). These differ from the registered values by 1.2e−12 and −3.2e−12, i.e. at their 12 printed digits.
- My K at d = 1e−3 equals the Q2/V13 values to all 20 printed digits:

| Point | K at d = 1e−3 |
|---|---|
| C0, +x | 1.1911145153985316e20 |
| C0, −x | 9.8130896336490276e19 |
| F5-R, +x | 1.8701851395038790e12 |

### 4. γ ambiguity and axis regularity: PASS (verified). Outside is regular; the singularity survives in both readings
**Axis, exact (sympy):**

| Reading | e^{2γ} at y = ±1, for every x > 1, α, β |
|---|---|
| Outside | **1** |
| Inside | **e^{−8β}** |

**At infinity (C0, x = 1e6):**

| Reading | e^{2γ} |
|---|---|
| Outside | 1 − 5e−12 |
| Inside | 1.4e−40 (not asymptotically flat in the standard form) |

**Elementary flatness.** Circumference/(2π × proper radius) for small circles about the axis, at C0, x = 1.5, 3 and 10, both axes:

| Reading | ε = 1e−6 | ε = 1e−10 | Interpretation |
|---|---|---|---|
| Outside | 1 + 2.6e−5 | 1 + 2.6e−9 | → 1, regular |
| Inside | 8.4488e19 | 8.4488e19 | = e^{4β}: a huge conical excess |

**Consequences:**
- **The outside reading is the regular one.** It is LGAC's and Brink's form, and it is what v9 computes.
- The ring's existence (B = 0) does not involve γ at all.
- **K in the inside reading is exactly e^{16β} × K in the outside reading** (5.0955354e79 at C0 on every ray, to 12 digits). This is because the curvature scales as e^{−2γ}. So the divergence and its d⁻⁶ law survive in both readings; only the coefficient changes.

### 5. Asymptotics and multipoles: PASS (verified, exact)
Method: Brink's ℰ on the upper axis, ξ = (1−ℰ)/(1+ℰ) expanded in 1/z, read with Fodor–Hoenselaers–Perjés (P_n = m_n for n ≤ 3).
- m₀ = M₀, m₁ = i S₁, m₂ = M₂ and m₃ = i S₃ of GLM Eq. 4, all **exactly** (symbolic difference 0).
- In terms of (M, χ, q): M₀ = M, S₁ = χM², and **M₂ = −M³(χ² + q) = −M(a² + qM²)**. So the quadrupole formula in the brief is correct.
- Numeric far field, C0, at x = 10⁴, 10⁶, 10⁸:
  - M from f: 0.99977 → 0.9999999771
  - J from g_tφ: 0.89979 → 0.8999999794 (expected aM = 0.9)
  - F5-R behaves the same way: M → 89, J → 7120.

### 6. Is the ring a coordinate artefact? PASS (verified numerically; the argument is standard)
**K and the Pontryagin invariant P,** C0, d = 1e−2 … 1e−5, on four rays, 60 digits. The 90-digit repeat agrees to 15 digits or better.

| Ray | K at d = 1e−5 | Last local slope of K |
|---|---|---|
| +x (ergo side) | 1.09e32 | −5.996 |
| −x (CTC side) | 1.09e32 | −6.004 |
| +y | −5.39e30 | −6.00001 |
| diagonal (−,+) | 1.67e31 | −6.002 |

Ricci/√K stays ≤ 8e−38 at d = 1e−5.

**Every direction, including K's nodes.** K changes sign 12 times around a small circle, which gives Q2's 12 nodal directions. But |I| = √(K² + P²), which equals 48|Ψ₂|² for a type D field, does not vanish there. On 36 angles at each d:
- min |I|·d⁶ = 9.04, 5.44, 5.392, 5.392 and max = 182, 119, 110, 109, for d = 1e−2 … 1e−5.
- So |I| ≥ 5.39 d⁻⁶ uniformly in a punctured neighbourhood. That means **every curve that approaches the ring sees a scalar invariant diverge**, not just the four rays.

**Why no coordinate change can remove it.** |I| is a scalar built from the Riemann tensor, so its value at a point does not depend on coordinates. If some coordinate change or extension made the metric C² across the ring, |I| would be continuous there. It is unbounded instead. So the ring is a scalar-polynomial curvature singularity, not a coordinate degeneracy (such as g_xx → 0 or f → ∞).

**Reachability (outside reading):**
- Proper distance along y = 0, which is totally geodesic by the reflection symmetry (g(x,y) = g(x,−y); ∂_y g = 0 on y = 0):
  - from x_s + 0.5 to x_s: **0.15469**
  - from x_s − 0.5 to x_s: **0.004729**
  - near the ring √g_xx ≈ 0.68563 |x − x_s|, so s ≈ 0.343 d², matching Q2's values.
- Equatorial geodesics from x_s + 10³ (≈ 437 M in Weyl r) inward to the ring, with g_xx ẋ² > 0 along the whole path (sampled):

| Geodesic | Proper time / affine parameter to the ring | Last 0.1 in x |
|---|---|---|
| Timelike, E = 1.2, L = 0 | 641.23 | 3.1e−4 |
| Timelike, E = 1.05, L = 0.5 | 1258.15 | |
| Null, E = 1, L = 0 | affine parameter 435.26 | |

These are incomplete geodesics ending on the curvature singularity. Near the ring, ẋ² ∝ |f|(ωE − L)²/(ρ²g_xx) ∝ δx⁻³, so the approach takes finite proper time for generic E, L.

### 7. Other attacks: none succeeded (verified unless marked)
- **Outside the horizon:**
  - x_s − 1 = 0.92983, so the ring is at Weyl ρ_s = 0.719448, z = 0, well away from the rod ρ = 0, |z| ≤ k.
  - The coordinate map (x, y) → (ρ, z) is non-degenerate there.
- **Uniqueness on the equator:** x+1+(x−1)ab has exactly one sign change for x − 1 ∈ (1e−6, 1e6) (4001-point log grid).
  - a = −b holds identically on y = 0. This is symbolic: ln(−a/α) − ln(b/α) = 0, and numerically |a+b|/|a| = 0.
  - My script printed max |a+b| = 1.6e+221730776. That is a reporting artefact of the absolute value near x → 1, where a and b are individually astronomically large; the relative check above is 0.
  - Q1's closed form β(cosh s − 1)²(cosh s + 2) = s + artanh κ holds at my x_s to 5e−111.
- **Not hidden behind a CTC region:** g_φφ > 0 on y = 0 for x ∈ (x_s + 1e−8, x_s + 1e4) (minimum 10.42). The ergo-side approach is causal and well behaved, as the geodesics above show.
- **Not a branch artefact:** R = √(x²+y²−1) > 0 at the ring, and every function is real-analytic in its neighbourhood except through 1/B.

## What an outside expert should still check
1. The original MN 1992 paper (unread here), in particular:
   - their γ normalisation, which fixes D1 at the source rather than through LGAC's restatement plus axis regularity;
   - whether their general solution fixes the a and b signs/branches as GLM do.
2. Whether this ring at x_s is already known in the literature. That is Q4's question, not settled here. GLM, LGAC and Brink each mention only the x = 1, y = 0 line singularity.
3. Q2's assumption about the first-order correction in the local analysis. Here it is supported numerically by the uniform |I| ~ d⁻⁶ bound, but it is not derived.
4. Q1's off-equator "no zeros for q > 0" result, which is numeric and partial. This test only rechecked the equator.
5. Whether the 9/4 (D2) and P₂/r² (D3) typos in Brink App. C were used in any later work that cites Brink's appendix.
