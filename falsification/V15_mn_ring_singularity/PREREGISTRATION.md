# V15: the Manko–Novikov equatorial ring singularity — existence map, analytic argument, geometry, literature

Pre-registered 2026-10-17, before any V15 computation. The user approved this work as Opus groundwork for the planned
email to Manko. It replaces an earlier planned review that stopped before producing results.

## Background (claims under test, from V13 FINDINGS "Equatorial ring singularity" and tabula 1e13988)
MN in the GLM 2008 (arXiv:0708.0628, Eqs. 2–7) / Brink 2008 (arXiv:0807.1179, App. C) form, as transcribed in
`V9_mn_obstruction_check/code/v9_nve.py` (the fleet's single transcription chain; tabula's is a second script from the
same papers). Parameters (M, a, β) with α = (−1+√(1−χ²))/χ, k = √(M²−a²), β = qM³/k³.
- B = A = 0 at one equatorial point: x_s = 1.92982825852 (C0: χ 0.9, q 0.95) and 1.92192258792 (F5-R: M 89, a 80, q 19/20).
- Kretschmann ~ d⁻⁶ there (two fleet implementations, validated on Kerr).
- (χ, q) = (0.6, 0.1024): B = A = 0 at x ≈ 1.0759622 (inside R < 0.5).

## Questions and methods (one agent each; the agents must not read each other's outputs until merge)
- **Q1, existence map.** Reduce B = 0 on y = 0 in closed form. Expected (to be checked, not assumed): with y = 0,
  P₁ = 0 and P₂ = −1/2, a(x,0) = −b(x,0) identically, so B = 0 ⇔ α² e^{4β(1 − x/R + x/(2R³))} = (x+1)/(x−1),
  R = √(x²−1). Map the (χ, q) region with a real root x_s > 1 (χ ∈ (0,1), q of both signs), and whether each root
  is in R < 0.5 or outside. Count roots. Off the equator: B is a sum of two real squares, so B = 0 needs
  x+1+(x−1)ab = 0 AND (1+y)a+(1−y)b = 0. Search for off-equator solutions (analytic argument plus a numeric scan).
- **Q2, curvature singularity, analytically and numerically.** An analytic argument (Kodama–Hikida gr-qc/0304064 §3.1.3
  style) that a curvature invariant carries an uncancelled inverse power of B at B = A = 0, with the leading exponent
  in the distance d. A numeric d-slope from NEW code that does not reuse `v13_kretschmann.py` or tabula's script
  (a different route: e.g. the Ernst-potential / Weyl–Papapetrou curvature formulae, or a tetrad Ψ₂ computation).
- **Q3, geometry.** The ergosurface A = 0, the boundary of g_φφ < 0, and B = 0: how they meet at the ring, and why
  g_φφ has a pole (not a zero) there. Analytic where possible.
- **Q4, literature.** Full text where reachable (no paywall circumvention): Manko & Novikov 1992 (CQG 9, 2477), GLM 2008,
  Brink 2008 (what "No singularities outside" refers to), Kodama–Hikida, and the MN/bumpy-BH/shadow/chaos papers since.
  Verdict "identified by X" or "not identified in (list)", listing what could not be read. Quotes ≤ 15 words.
- **V, adversarial verifier (after Q1–Q4).** Try to show the ring is a transcription or coordinate artefact. Required
  tests: (i) vacuum: Ricci = 0 to ≥ 30 digits at ≥ 10 points including d = 1e−3 from the ring at C0 (a wrong γ or a/b
  exponent would break this); (ii) asymptotics: the mass, spin and quadrupole read off the far field match the
  published MN moments; (iii) the GLM Eq. 3h "−2 inside/outside the sum" γ ambiguity: which reading is Ricci-flat;
  (iv) a different coordinate (Weyl ρ, z or a curvature invariant) cannot be a removable degeneracy if a scalar diverges.

## Pass / fail criteria
- Q1 PASS if the closed-form equation reproduces x_s at C0 and F5-R to ≤ 1e−9 and the (0.6, 0.1024) root to ≤ 1e−6.
  Any other root found at C0 or F5-R is reported (it would change the V13 domain discussion).
- Q2 PASS if the analytic leading exponent and the new code's fitted slope (d = 1e−2 … 1e−5, both sides) agree within
  ±0.1, and the new code agrees with the V13 K values at d = 1e−3 (C0) to ≤ 1e−6 relative. A disagreement is a FAIL to
  be reported, not tuned away.
- V PASS (ring stands) only if (i)–(iii) all pass. If Ricci ≠ 0 anywhere, the ring claim is SUSPENDED and every V13
  statement that relies on the metric is flagged.
- Q3, Q4: descriptive; each claim marked verified / partially verified / unverified.

## Controls
- Positive (known answer): Kerr (β = 0) Kretschmann closed form, 48M²(r⁶ − 15a²r⁴cos²θ + 15a⁴r²cos⁴θ − a⁶cos⁶θ)/Σ⁶,
  at 3 points; Kerr's own ring singularity (B = 0 needs x < 1 for β = 0, since α² < 1) is where the closed form says.
- Negative: Kerr has NO B = 0 point with x > 1 (α² < 1 < (x+1)/(x−1)); a generic MN point far from the ring (C0, x = 5,
  y = 0.3) has finite, smooth K and Ricci = 0.

## Outputs
`code/`, `results/` (JSON with versions), `FINDINGS.md`. Precision: mpmath ≥ 50 digits for curvature near the ring.
Python env: a pinned venv (versions recorded in each result; see ENVIRONMENT/).

## Amendment Q2b (2026-10-17 evening, before any work on it): derive Q2's assumed step
- **Why:** Q2's leading exponent rests on an ASSUMPTION. The first correction to the local vacuum pp-wave was taken
  to behave as in Kerr, giving K ≈ 48 cos 6θ_W / (ρ₀²Λ²|W|⁶). It is supported numerically (K·δ⁶ matches to 3.5e−7
  at 4 rings) but not derived. An outside expert would ask for the derivation.
- **Task (attempted analytically first):** derive the first-order correction from
  the Ernst equation near a simple pole of ℰ (holomorphic in W). Show analytically that the Weyl scalar Ψ₂ (or I) has
  leading term ∝ W̄⁻³ with the stated coefficient, including Λ in closed form.
- **PASS:** the derived coefficient reproduces the numeric K·δ⁶ at the 4 rings of Q2 to ≤ 1e−6 relative, and
  reduces to the exact Kerr result at β = 0.
- **FAIL / INCONCLUSIVE** is reported as such. The singularity verdict does not depend on Q2b, because the invariant
  divergence is shown numerically by V. Q2b only decides whether the exponent is derived or empirical.

## Q4 (registered 2026-10-23, before any computation): Komar mass of the equatorial ring
Motivation: an expert reader pointed to the "massless ring singularities" in vacuum stationary solutions (Manko & Ruiz,
PTP 125 (2011) 1241, p. 1250; Manko, PTP 127 (2012) 1057 = arXiv:1110.6564). Is our MN ring of that kind?
- Method: Komar mass M = −(1/8π)∮ *dξ♭ for ξ = ∂_t, on tori {closed curve C in (x, y)} × φ, at χ = 0.9, q = 0.95
  (β = q/κ³), with M = 1. In Weyl–Papapetrou prolate form the integrand reduces exactly to
  (*F)_{xφ} = −2(g_φφ ∂_y g_tt − g_tφ ∂_y g_tφ)/(k(x²−1)) and (*F)_{yφ} = 2(g_φφ ∂_x g_tt − g_tφ ∂_x g_tφ)/(k(1−y²)),
  using √|g| = k³|B|e^{G−2ψ}/(1−α²)² and the regularised g_φφ of v13. Evaluated in mpmath at 30 digits.
- Controls (all must pass, else the result is void):
  - C1: the asymptotic contour (x = 200, y ∈ [−1, 1]) gives the total mass 1 to ≤ 1e-6. This also fixes the
    overall sign and orientation; it is done for MN AND for Kerr (β = 0).
  - C2: a small circle around a REGULAR point (x = 3, y = 0.3) gives 0 to ≤ 1e-12.
- Measurement: circles of radius r ∈ {1e-2, 3e-3, 1e-3} in (x, y) around the ring point (x_s, 0).
- Verdict:
  - "massless ring" if |M_ring| ≤ 1e-10 at every radius (conservation in vacuum means it is radius-independent);
  - "ring carries Komar mass m" if the value is radius-independent and nonzero;
  - otherwise "inconclusive" (e.g. the contour meets the ergoregion/CTC structure).
