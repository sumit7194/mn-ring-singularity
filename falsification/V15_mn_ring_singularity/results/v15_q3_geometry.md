# V15 Q3: how the ergosurface, the CTC boundary and B = 0 meet at the Manko–Novikov ring

Code: `code/v15_q3_geometry.py` (my own; it does not import the V9 or V13 modules). Data: `results/v15_q3_geometry.json`.
Figures: `results/v15_q3_c0.png` (C0, Weyl ρ ∈ (0, 1.6), |z| < 1.2, z pointing up), `results/v15_q3_c0_zoom.png`
(|ρ−ρ_s|, |z| < 0.02) and `results/v15_q3_c0.svg` (both images with axes and a legend). matplotlib is not installed in the venv,
so the PNGs come from my own zlib writer. Environment: Python 3.14.5, sympy 1.14.0, mpmath 1.3.0, numpy 2.4.6, scipy 1.18.0,
macOS-26.5-arm64. Precision: mpmath 50 digits near the ring, 30 digits for the global grid (signs rechecked at 70/80 digits).

Notation: small a, b are the MN functions; capital A, B, C are the GLM polynomials;
f = e^{2ψ}A/B, ω = 2k e^{−2ψ}C/A − w0, w0 = 4kα/(1−α²), ρ² = k²(x²−1)(1−y²).
g_φφ = −fω² + ρ²/f = N/B with N = e^{−2ψ}k²Q + 4k w0 C − w0² e^{2ψ}A, where AQ = (x²−1)(1−y²)B² − 4C². I re-derived Q by
polynomial division (remainder 0, 46 terms). Let λ = 1 − 3β/R_s³ with R_s = √(x_s²−1), and let d, θ be the Weyl distance
and angle from the ring (θ = 0 points outward along the equator).

## 1. The exact equatorial algebra (symbolic, holds for all parameters) [verified]
- On y = 0, a = −b holds identically (sympy). Substituting it into the polynomials gives exact factorisations
  with F = (x−1)a² − (x+1):
  A = F·[a²(x+1) − (x−1)], B = F², C = 2ax·F, Q = F·[a²(x−1)³ − (x+1)³].
  The ring is F = 0, i.e. a² = (x+1)/(x−1). The second square in B is zero on y = 0 for every x.
- So at **every** ring point, for every (χ, q): A = C = Q = 0, and therefore N = 0. The ergosurface always passes through
  the ring. On the equator, N and A vanish to first order in F and B to second order.
- At F = 0: A/F = 4x/(x−1), Q/F = −4x(x+1), and N/F = −(4x/(x−1)) e^{2ψ} ω_r², with ω_r = k e^{−2ψ}√(x²−1) − w0
  (sympy-verified). ω_r is the ring limit of ω. On the equator this gives g_φφ/f → −ω_r².
- Parity: under y → −y with a → −b and b → −a, the polynomials A, C, Q (and so N) are even (sympy). B = S1² + S2², where
  S1 = x+1+(x−1)ab is even and S2 = (1+y)a + (1−y)b is odd. So nothing at the ring has a term linear in y, and B has no
  mixed term (B_xy = 0).
- Transversality: on y = 0, dh/dx = −3/(2R⁵) and ∂_y(exponent of a) = 3/(2R³) (sympy). It follows that
  dF/dx = 2λ/(x−1), ∂S1/∂x = −2λ/(x−1) and ∂S2/∂y = 2aλ. Everything below is non-degenerate exactly when λ ≠ 0,
  i.e. R_s³ ≠ 3β. At R_s³ = 3β the root of F is double, which is where two ring roots merge.

## 2. Local forms at C0 and F5-R [verified: two derivative implementations, mpmath numeric and sympy symbolic, agree to 3e-49 (C0) and 3e-43 (F5-R); the closed forms of §1 agree to 6e-50]
| | C0 (M = 1) | F5-R (M = 89) |
|---|---|---|
| x_s (refined) | 1.929828258521150008 (given value +1.2e−12) | 1.921922587916759772 (−3.2e−12) |
| ρ_s | 0.7194477398696 | 64.00975836566 |
| λ = 1 − 3β/R_s³ | −6.65327 | −6.66084 |
| A ≈ A_x δx + ½A_yy y² | −118.806 δx − 1121.23 y² | −120.494 δx − 1130.87 y² |
| B ≈ ½B_xx δx² + ½B_yy y² (B_xy = 0) | 204.798 δx² + 557.917 y², positive definite | 208.799 δx² + 562.461 y², positive definite |
| C ≈ C_x δx | −98.046 δx | −98.882 δx |
| Q ≈ Q_x δx | 323.656 δx | 324.585 δx |
| N(x_s, 0) | 1e−47 (zero at working precision) | 8e−44 (zero) |
| N ≈ N_x δx + ½N_yy y² | 1126.16 δx + 1778.17 y² | 9.05376e6 δx + 1.42083e7 y² |
| ω_r | 11.0241 | 982.747 (= 11.042 M) |

In Weyl coordinates the leading forms are isotropic:
B ≈ K_B d², with K_B = 4λ²(x+1)/((x−1)k²x²); f ≈ K_f cosθ/d, with K_f = 2k x² e^{2ψ}/(R λ); and g_φφ ≈ K_g cosθ/d, with
K_g = −ω_r² K_f. At C0, K_B = 788.46, K_f = −0.023060 and K_g = 2.8025. At F5-R, K_f = −2.0504 and K_g = 1.9802e6.
At d = 1e−6 all three ratios to these forms are 1 within 1e−5.

## 3. Why g_φφ has a pole, not a zero, and its sign [verified at C0 and F5-R; the sign rule is derived in general]
- f = e^{2ψ}A/B has a **simple pole**: A vanishes linearly (A_x ≠ 0 because λ ≠ 0) and B quadratically (positive
  definite). Along every ray not tangent to the ring's vertical direction, C/A tends to C_x/A_x, so ω → ω_r stays finite
  and nonzero. Then g_φφ = −fω² + ρ²/f ≈ −ω_r² f, and the −fω² term dominates. Equivalently, g_φφ = N/B with N ~ d and
  B ~ d².
- What does vanish at the ring is ρ²/f. This is the norm of the Killing combination ∂_φ + ω∂_t
  (g_φφ + 2ωg_tφ + ω²g_tt = ρ²/f ~ d). The (t, φ) block keeps its determinant at −ρ², which stays finite. Its singular part
  has rank one: −f (dt − ω_r dφ)². So g_tt, g_tφ and g_φφ all diverge like 1/d, and nothing degenerates to zero.
- Numerical checks on rays at θ = 0°…180° and d = 1e−2…1e−6: the slope of log|g_φφ| against log d is −1.000 ± 0.002 on
  every non-tangent ray. g_φφ/(−ω_r² f) → 1. The regularised N/B and the raw −fω² + ρ²/f agree to ≤ 2e−36. Along the
  tangent direction θ = 90° (δx ≈ +0.71 y²), g_φφ stays finite (g_φφ/(−ω_r² f) → 0.225).
- **The sign depends on direction.** At C0 and F5-R (λ < 0), g_φφ → +∞ from the outer side (θ < 90°, where f → −∞, inside
  the ergoregion) and g_φφ → −∞ from the inner side (θ > 90°, where f → +∞, the CTC side). In general f·g_φφ ≈ −ω_r² f² < 0,
  so near the ring the CTC side and the f < 0 side are always opposite. sign(λ) decides which physical side is which.
- Kerr control (β = 0, ring at x = −M/k): the same code reproduces the Boyer–Lindquist g_φφ and f to 1e−45. It gives
  ω_r = a_spin = 0.9 and g_φφ/(−ω_r² f) → 1. So the MN ring has the same local structure as Kerr's ring singularity.

## 4. The triple meeting [the meeting at the ring is verified as generic algebra; the curve structure is verified at C0/F5-R and generic when λ ≠ 0 and ω_r ≠ 0]
- **Generic.** For every (χ, q) with a ring root, A = C = Q = N = 0 there (§1). So the ergosurface (A = 0), the CTC boundary
  (g_φφ = 0, i.e. N = 0) and B = 0 always meet at the ring.
  The algebraic reason: on the equator, a = −b makes all of A, B, C, Q contain the factor F = (x−1)a² − (x+1).
- **B = 0** is an isolated point in the (x, y) plane, i.e. a circle in space, because its quadratic form is positive definite.
- **The ergosurface and the CTC boundary both pass through the ring as smooth curves**, not just touching it at one point.
  Their gradients at the ring are ∝ (1, 0) with A_x ∝ λ and N_x ∝ λω_r², both nonzero, so the implicit function theorem
  gives one smooth branch each. Both curves are tangent to the direction normal to the equatorial plane (vertical in
  Weyl ρ–z) and so tangent to each other, but their curvatures differ. Write each curve as δx = c y² (equivalently δρ = κ δz²).
  At C0: A = 0 has c_A = −9.4375 (κ_A = −7.306) and N = 0 has c_N = −1.5790 (κ_N = −1.646); exact roots at y = 3e−4 match
  to 2e−5. At F5-R: c_A = −9.385 and c_N = −1.569. Both curves bend to the inner side.
- **Sectors at the ring (C0, F5-R).** Inner side, δx < c_A y²: CTC only (f > 0, g_φφ < 0). The thin lens between the curves:
  f < 0 **and** g_φφ < 0. Outer side, δx > c_N y²: ergoregion only. No region with f > 0 and g_φφ > 0 touches the ring;
  in Kerr such a region sits between the two curves. That c_A < c_N < 0 holds in general is **unverified**; I checked only
  these two points.
- Exceptions to the generic picture: λ = 0 (a double root, R_s³ = 3β), where B, A and N vanish to higher order; and
  ω_r = 0, where N_x = 0. Here ω_r is 11.0 M at both points.

## 5. The CTC region at C0 [partially verified: grid h = 0.005 in (ρ, z); extremal edges refined by bisection and golden-section search at 40 digits; signs rechecked with 0/300 random mismatches at 70 digits and 0/2410 near-rod mismatches at 80 digits]
- **It is bounded.** It lies in ρ < ρ_s = 0.71945 and |z| ≤ 0.91976 (M = 1). The supremum of ρ is the ring itself, approached
  but not attained. Outside it: g_φφ > 0 at all sampled angles on Weyl radii 2, 5, 20 and 100, and g_φφ/ρ² > 0 on the axis
  for |z| > k = 0.43589.
- Three components at grid resolution (the region is mirror-symmetric in z):
  1. **The equatorial lobe**, which has the ring at its outer tip. On the equator it covers x ∈ (1.62844, x_s), i.e.
     ρ ∈ (0.56022, 0.71945). Its arms wrap around a normal lens (f > 0, g_φφ > 0) and reach inward to ρ ≈ 0.05. It extends to
     |z| ≤ 0.24705 (maximum at ρ = 0.4933). Most of the lobe lies inside the ergoregion (f < 0 and g_φφ < 0, purple in the
     figure). Only the tip at the ring is CTC-only.
  2. **Two mirror off-equator lobes.** They reach ρ ≤ 0.60990 (maximum at |z| = 0.7092) and |z| ≤ 0.91976 (maximum at
     ρ = 0.4476). They reach the horizon rod (ρ → 0, |z| < k) for |z| ≲ 0.31–0.33. That crossing point drifts slowly as x → 1:
     y* = 0.752, 0.732 and 0.714 at x − 1 = 1e−8, 1e−12 and 1e−16, so it is not converged. On the rod g_φφ ∝ (x−1) → 0,
     with an enormous coefficient. Their closures contain the singular horizon-equator point x = 1, y = 0: a polar scan at
     r = 0.05 down to 0.005 shows g_φφ < 0 in the 36°–90° sector.
- Points with both f < 0 and g_φφ < 0 cover 5204 of the 153,920 grid pixels. The ergoregion's equatorial lobe reaches
  ρ = 1.24577 (x = 3.02788).
- Whether the three components really are disjoint, rather than joined by necks thinner than 0.005, and the fine structure
  very near (ρ, z) = (0, 0), are **unverified**. The polar scan suggests the equatorial lobe does not reach r = 0.05 there.
