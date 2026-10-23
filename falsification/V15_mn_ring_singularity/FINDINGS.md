# V15 FINDINGS: the Manko–Novikov equatorial ring singularity

Pre-registration: PREREGISTRATION.md (bf73ce2). Four independent Opus agents (Q1–Q4), plus the verifier (V) after them.
Status: COMPLETE. Q1–Q4, V and Q2b are in. **Verifier verdict: the ring stands. Q2b: the d⁻⁶ law is derived.**

## Q1: existence map — PASS (all registered checks)
- **Equatorial reduction (verified, sympy plus 50-digit numerics).**
  - On y = 0: P₁ = 0, P₂ = −1/2, and a(x,0) = −b(x,0) identically.
  - B(x,0) = (x+1−(x−1)b²)², and A(x,0) = −(x+1−(x−1)b²)((x+1)b²−(x−1)). So A = 0 wherever B = 0.
- **Closed form.** With s = arcoth x, κ = √(1−χ²) and β = q/κ³, B = 0 on the equator ⇔
  **β (cosh s − 1)² (cosh s + 2) = s + artanh κ.** The ratio LHS/RHS is strictly monotone in s, so:
  - **exactly one equatorial ring for every q > 0 and every χ ∈ (0, 1);**
  - none for q ≤ 0 (Kerr, q = 0, has none: negative control PASS).
  - Bridge spot-check of the closed form at the V13 roots: residual 2.8e−12 (C0) and −8.0e−12 (F5-R), consistent with
    the 12 printed digits.
- **Registered roots.** C0 1.929828258521150 (Δ 1.2e−12); F5-R 1.921922587916760 (Δ 3.2e−12); (0.6, 0.1024)
  1.075962201195 (Δ 1.2e−9). No other root at any of the three. PASS (≤ 1e−9, ≤ 1e−9, ≤ 1e−6).
- **R < 0.5 boundary.** The ring is at R < 0.5 iff q < q*(χ) = (1−χ²)^{3/2}[asinh 2 + artanh√(1−χ²)]/(2(1+√5)).
  At χ = 0.9, q* = 0.0245. Confirmed on a 103 × 401 grid with 0 mismatches.
- **Off the equator.**
  - q > 0: no zeros (numeric scan plus edge asymptotics; partially verified, not a proof).
  - q < 0: exactly one ±y pair of off-equator B = A = 0 rings for every χ. Verified on the unreduced B at 60 digits.
    Example: χ 0.9, q −1 gives x ≈ 1.69, y ≈ 0.80.
  - Direct minimisation of B finds exactly one zero at each of 8 test points.
  - The q < 0 rings do not affect V13 (all its points have q > 0), but the V13 domain discussion should mention them.

## Q2: curvature singularity — singularity CONFIRMED; the strict per-side slope criterion FAILS on one side (reported as registered)
- **Analytic, partially verified:**
  - B = |A₊|², A = Re(A₋Ā₊) and Ernst ℰ = e^{2ψ}A₋/A₊. ℰ satisfies the Ernst equation to ~1e−99.
  - On the equator A₋ = −2(1+a²) ≠ 0, so ℰ and f have simple poles at the ring.
  - The Ernst equation forces the pole to be holomorphic in W = (ρ−ρ₀)+iz. Predicted |∂_yA₊/∂ₓA₊| = R_s matches to
    25 digits.
  - The leading local metric is a vacuum pp-wave (all invariants vanish), so naive power counting is wrong.
  - ASSUMED, not derived: the first correction behaves as in Kerr. That gives
    K ≈ 48 cos 6θ_W / (ρ₀²Λ²|W|⁶), i.e. **K ∝ d⁻⁶ in coordinate distance (proper distance s ∝ d², so ∝ s⁻³)**, with
    12 nodal directions.
  - The coefficient law is checked numerically at 4 rings (C0, F5-R, (0.6, 0.1024), Kerr's ring): K·δ⁶ matches to
    3.5e−7. The ring is asymptotically type D.
- **Numerical route, new and independent of v13_kretschmann:**
  - 2+2 block-Riemann formulation with exact sympy derivatives, mpmath at 100 digits (a 150-digit repeat agrees to
    77 digits).
  - Kerr control: matches the closed form to ≥ 98 digits at 3 points. Negative control (C0, x 5, y 0.3):
    K = 0.08357926…, finite and smooth.
  - **Ricci:** R_abR^ab / K ≤ 7e−160 at the ring points (vacuum holds to working precision).
- **Registered criteria:**
  - Agreement with V13 at C0, d = 1e−3: **PASS** (Δ 1.2e−12 and 9.9e−14). All 16 V13 values agree to ≤ 1.4e−10.
  - Exponent vs fitted slope (d 1e−2…1e−5) within ±0.1:
    - x > x_s: −5.930 (C0), −5.929 (F5-R): pass.
    - **x < x_s: −6.128 (C0), −6.121 (F5-R): FAIL by ~0.03 beyond the band.**
    - Off-axis rays: x = x_s, ±y −6.049; diagonal (+,+) −5.967; **diagonal (−,+) −6.310 (outside the band).**
  - **Bridge reading:** the prereg says "both sides" without saying pooled or per side. We take the STRICTER, per-side
    reading, because choosing the lenient one after seeing the data would be tuning. **Criterion: FAIL as registered.**
  - **Cause (post hoc, labelled as such):** the d = 1e−2 point is pre-asymptotic. On the inner side K is still
    negative there. The local slopes over the last decade are −6.004, −5.996, −6.002, −5.998 and −6.00001 on every ray.
  - So the exponent −6 is supported by the local slopes and the analytic coefficient match. The registered fit window
    was simply too wide.
  - This does not change the verdict that the ring is a curvature singularity. K grows without bound on every ray, by
    ~6 orders of magnitude per decade in d.
- **Caveat carried to V:** Q2 and V13 share the V9 transcription. Their agreement checks the curvature codes, not the
  transcription.
  - Ricci-flatness cannot settle the GLM Eq. 3h "−2" ambiguity: moving the −2 outside the sum multiplies e^{2γ} by
    e^{8β}, which is still vacuum but has a conical axis.
  - So V must decide (iii) by **axis regularity** (e^{2γ} = 1 on the axis, elementary flatness) and an independent
    transcription source. Ricci alone is not enough. This sharpens the prereg's test (iii).

## Q3: geometry at the ring (verified at C0 and F5-R; generic structure verified symbolically)
- **Local forms at C0:**
  - A ≈ −118.8 δx − 1121 y² (linear in δx);
  - B ≈ 204.8 δx² + 557.9 y² (positive-definite; ≈ 788 d² in Weyl coordinates);
  - the numerator N of g_φφ: N(x_s, 0) = 0, N ≈ 1126 δx + 1778 y².
- **Pole, not zero:** f = e^{2ψ}A/B ∝ 1/d. ω tends to a finite nonzero limit (11.02 at C0; 11.04 M at F5-R), so
  g_φφ ≈ −ω² f ∝ 1/d.
  - g_φφ → +∞ from the outside (ergoregion side) and −∞ from the inside (CTC side), and stays finite along the
    vertical tangent.
  - What vanishes is ρ²/f, the norm of ∂_φ + ω∂_t. This is the same structure as Kerr's ring, and the code
    reproduces Kerr's g_φφ to 1e−45.
- **Generic triple meeting (symbolic):** on the equator, with F = (x−1)a² − (x+1):
  - A, C and Q are each divisible by F, and B = F², so A = C = Q = N = 0 at every ring, for every (χ, q).
  - The ergosurface and the g_φφ = 0 curve both pass vertically through the ring and are tangent to each other.
  - The construction fails only if λ = 1 − 3β/R_s³ = 0 or ω_ring = 0.
  - The sector ordering (CTC-only inside, ergo-only outside, a thin lens with both) is verified at C0 and F5-R only.
- **CTC region at C0** (partially verified): bounded, with ρ < ρ_s = 0.71945 and |z| ≤ 0.91976. The ring is its
  outermost point. Three lobes: one equatorial, and a mirror pair touching the horizon rod. Whether the lobes are
  disjoint, and the structure near x = 1, are unverified.

Files: code/v15_q{1,2,3}_*.py; results/v15_q1_existence.json, v15_q1_map.svg, v15_q2_curvature.json,
v15_q2_analytic.md, v15_q3_geometry.{json,md}, v15_q3_c0*.png/svg. Env: the conjecture_machine venv (Python 3.14.5,
sympy 1.14.0, mpmath 1.3.0, numpy 2.4.6), as recorded in each JSON.

## Q4: literature (results/v15_q4_literature.md, screening record results/v15_q4_screen.json)
- **Verdict: no paper we could legitimately read identifies a curvature singularity off the horizon in MN.**
  - Full text read closely [verified]: Gair–Li–Mandel 2008; Brink 2008; Lukes-Gerakopoulos–Apostolatos–Contopoulos
    2010; Contopoulos et al. 2011 and 2012; Lukes-Gerakopoulos & Contopoulos 2013; Bambi & Barausse 2011; Bambi &
    Lukes-Gerakopoulos 2013; Johannsen 2013; Destounis & Kokkotas 2021; Staelens et al. 2023; Wang–Chen–Jing 2018;
    Dey et al. 2022; the Yunes Living Reviews article; the Geyer 2013 thesis; four Manko-group papers.
  - Keyword screen of full text [partially verified]: 284 further papers. 40 non-arXiv works were checked on metadata
    or abstract only. Text inside figures cannot be screened.
- **Brink 2008 Fig. 10 ("No singularities outside")** [verified on the rendered page]: a label inside a cartoon,
  making a family-wide claim with no parameter point or computation behind it. Brink's Fig. 6 parameters work out to
  χ 0.9, q 0.95, the same as C0.
- **Not read:** MN 1992 (CQG 9, 2477; the IOP page is behind a captcha, not circumvented, and it is closed access
  elsewhere); Castejon-Amenedo & Manko 1990; Manko 1990 (CQG 7 L209); Manko et al. 1998.
  - **So "not identified in MN 1992" is UNVERIFIED.** Johannsen 2013 cites MN 1992 only for "naked singularity".
- **Closest precedents** [verified]:
  - Kodama & Hikida 2003 (Tomimatsu–Sato δ=2): Ψ₂ ∝ B⁻³ at the equatorial B = 0 point where the inner ergo-boundary
    meets the equator, with g_φφ having a pole there. This is the same structure as Q2/Q3.
  - Manko 2011 (arXiv:1110.6564): a vacuum ring singularity sits where the stationary-limit surface touches the
    CTC-region boundary. This is the Q3 geometry, but it has not been applied to MN by anyone.
- **Unrecognised traces at ρ ≈ 0.72** [our inference; unverified]: GLM's touching point of the CTC and bound regions;
  the E- and L-independent inner boundary of allowed motion in Contopoulos et al. 2011/2012; the deep
  effective-potential hole in Lukes-Gerakopoulos et al. 2010; Geyer's g_φφ blob.
- **For the Manko email:** the honest framing is "we could not find this in the literature we could read, and we
  could not read MN 1992, so we're asking whether it is already known". It is not "new".

## V: adversarial verifier — RING STANDS (all 7 tests PASS; results/v15_V_verifier.md)
1. **Independent transcription** (from the arXiv LaTeX of GLM 2008, Brink 2008 App. C and LGAC 2010, made before
   reading v9_nve.py): PASS. The papers disagree with each other in three places:
   - **GLM Eq. 3h "−2":** typeset inside the sum in GLM. LGAC move it outside "to avoid any confusion"; Brink's form
     is equivalent to outside.
   - **β² coefficient in γ:** GLM and LGAC have 3/2; **Brink has 9/4, a typo** (it fails vacuum: Ricci/√K up to 0.58).
   - **b exponent:** **Brink prints P₂/r² where GLM and LGAC have P₂/R³, a typo** (Ricci/√K up to 5.7).
2. **Comparison with v9_nve:** f, ω and e^{2γ} agree to 5e−79 (C0) and 7e−77 (F5-R). v9 types GLM's "−2 inside" and
   adds 8β, which equals the "outside" reading. PASS.
3. **Vacuum:** full 4-D Ricci at 17 points, including d = 1e−3 from the ring. Ricci/√|K| ≤ 9e−87 at 100 digits
   (round-off). Brink's Ernst potential satisfies the Ernst equation to 7e−54 and reproduces GLM's f and ω.
   Independent K equals Q2/V13 to 20 digits. PASS.
4. **Axis regularity:** the "outside" reading (the one v9 uses) has e^{2γ} = 1 on both axes and at infinity, so it
   is elementary-flat. The "inside" reading is conical: the circumference ratio is e^{4β} ≈ 8.4e19 at C0. **The
   singularity survives in BOTH** (K differs by exactly e^{16β}, with the same d⁻⁶ law). PASS.
5. **Multipoles:** from the axis Ernst potential, mass, spin, quadrupole −M(a² + qM²) and current octupole equal GLM
   Eq. 4 exactly. PASS.
6. **Not a coordinate artefact:**
   - Local slopes are −6.000 ± 0.004 on 4 rays.
   - K has 12 nodal directions, but the invariant √(K² + P²) (P = Pontryagin) stays ≥ 5.39·d⁻⁶ at all 36 angles
     sampled, so a scalar invariant diverges along every approach.
   - **The ring is reachable:** proper distance is finite (0.1547 from x_s + 0.5), and equatorial timelike and null
     geodesics reach it at finite proper time or affine parameter.
   - PASS.
7. **Other attacks:** the ring is outside the horizon (x_s − 1 = 0.930, ρ = 0.719); it is the only equatorial root;
   there are no CTCs on the approach from outside (g_φφ ≥ 10.4); it is not a branch-cut effect. GLM and LGAC mention
   only the x = 1, y = 0 singularity. PASS.

**Prereg V criterion ((i)–(iii) all pass): MET. The ring claim stands** for the metric as published in GLM and LGAC
(and in Brink, once its two typos are corrected).
- **Still for an outside expert:**
  - MN 1992's own γ normalisation and formulae (the original was not read);
  - whether the ring is already known;
  - Q2's assumption about the first-order correction (supported numerically here, not derived);
  - Q1's "no off-equator zeros for q > 0", which is numerical only.
- **Side finding (for any outreach):** Brink 2008 App. C has two typos (β² coefficient 9/4 instead of 3/2; P₂/r²
  instead of P₂/R³ in b). Both are caught by the vacuum test.

## Registered-criteria summary
- Q1: PASS.
- Q2: V13 agreement PASS; slope band FAIL as registered (strict per-side reading; the cause is a pre-asymptotic
  d = 1e−2 point; the verifier's local slopes are −6.000 ± 0.004).
- V: PASS.
- Q3, Q4: descriptive.

## Q2b: the assumed step is now DERIVED — PASS (results/v15_q2b_derivation.md, results/v15_q2b_check.json)
- **Closed form:** with ℰ = c/W + … at the pole and W = (ρ−ρ₀) + iz,
  I = 3(c/c̄)²/(ρ₀²Λ²W⁶)·(1 + O(|W|)) and K = 48 cos(6θ_W − 4 arg c)/(ρ₀²Λ²|W|⁶), with Λ = lim e^{2γ}/(f|W|²).
  - Only arg c, ρ₀ and Λ enter; the subleading Laurent data cancel.
  - c is real at equatorial rings (reflection symmetry), so this is exactly Q2's assumed law.
  - Λ is fixed globally by γ's axis normalisation.
- **Mechanism (derived):** an exact identity I = (4/H²)(𝒬_WW𝒬_W̄W̄ + 3𝒬_WW̄²) with H = e^{2γ}/f.
  - The holomorphic leading pole leaves only 𝒬_WW, a null pp-wave with zero invariants.
  - The 1/ρ term of the Ernst equation generates 𝒬_WW̄ at first order. So the invariant appears two orders below
    naive counting: |W|⁻⁶, not |W|⁻⁸.
- **PASS (criterion ≤ 1e−6 relative):** closed form vs Richardson-extrapolated numeric K·δ⁶:
  - C0 109.0013104836 (9.7e−35);
  - F5-R 1.70976083e−6 (1.0e−34);
  - (0.6, 0.1024) 14.60186064 (1.4e−27);
  - Kerr 6998.10468 = 48M²/k⁶ (3.6e−75; bridge arithmetic check 48/0.19³ = 6998.10). The symbolic Kerr reduction is
    exact.
- **Extra, unregistered phase test:** the (c/c̄)² factor holds at an off-equator q < 0 pole (χ 0.9, q −1;
  arg c = −1.02) to 1.7e−11.
- **Status of each step:**
  - derived: the identity, the perturbative solution, the cancellations, the coefficient, Λ, the Kerr reduction;
  - cited: the Ernst and γ equations (Ernst 1968; checked here to 1e−84), Wald 1984 Eq. C.3.6, O'Neill 1966;
  - hypothesis: ℰ = N/D with real-analytic N and D, a simple zero of D and N ≠ 0 (holds at all 4 rings).
- **Consequence:** the d⁻⁶ law is now DERIVED (under the stated hypothesis), not assumed. 
- **Note:** the derivation also gives the off-equator q < 0 poles the same |W|⁻⁶ law, which suggests they are
  curvature singularities too. They have not been checked separately; stated as a hypothesis for future work.

## Literature additions (2026-10-18; leads from an AI search answer supplied by the author, each verified here)
- **MN 1992 bibliographic record** verified on Crossref: "Generalizations of the Kerr and Kerr-Newman metrics
  possessing an arbitrary set of mass-multipole moments", CQG 9 (11), 2477–2487, Nov 1992, DOI
  10.1088/0264-9381/9/11/013. There is still no open copy.
- **GLM 2008 corrects a typo in MN 1992's Eq. (14)** (the multipole moments), "following [19]" = H. Fang, PhD thesis,
  Caltech 2007. Verified in the GLM arXiv PDF.
  - Our multipole check (V test 5) is against GLM's corrected Eq. 4.
  - Fang's thesis is a possible open-access route to MN 1992's formulae. Not yet retrieved.
- **GLM's only singularity statement** is the "circular line singularity at x = 1, y = 0 (ρ = z = 0)", credited to
  Brink by private communication. That is ON the horizon rod, not our ring at x_s ≈ 1.93. Verified, and consistent
  with Q4. That search answer misdescribed this as an "equatorial ring singularity discussion".
- **NEW and relevant:** Chruściel, Maliborski & Yunes, "The structure of the singular ring in Kerr-like metrics", PRD
  101, 104048 (2020), arXiv:1912.06020. Verified on the arXiv record (the search answer gave wrong authors).
  - For Kerr-like metrics, K vanishes at the singular set along particular directions, but a combination of
    invariants diverges from every direction. This is exactly the structure found here (V test 6: 12 nodal
    directions of K; √(K² + P²) diverges everywhere).
  - It treats Kerr-family metrics, not MN, so it doesn't bear on novelty. It is a good citation for the email and any
    write-up.
- Lukes-Gerakopoulos, Contopoulos & Apostolatos 2014 (arXiv:1408.4697, "Non-linear effects in non-Kerr spacetimes"):
  real, about bumpy-BH dynamics, and not about singularities.

## Fang 2007 thesis check (results/v15_fang2007_check.md)
- **Source:** H. Fang, *Topics in Gravitational Physics: Tidal Coupling in Gravitational Wave Searches and Mach's
  Principle*, Caltech PhD 2007, DOI 10.7907/D99H-J577. Open access at https://thesis.caltech.edu/1915/01/thesis.pdf.
  The MN material is in Appendix 3.A.1, pp. 74–78.
- **The MN 1992 Eq. 14 typo** (footnote 20, p. 77): in M₂ the last α₂ should be α². It is in the moment formula only,
  not the metric.
- **Agreement:** our ψ, a, b, f, α, k and the β² term of γ′ agree with Fang Eqs. (3.64)–(3.70), which Fang presents as
  evaluations of MN 1992 Eqs. (11)–(13). Combined checks agree to 40 digits at 3 points, and ω matches Fang's
  large-r expansion.
- **γ normalisation:** Fang's closed form equals β(Σ − 2), i.e. "−2 outside the sum". That is our axis-regular
  reading; GLM's literal print would give β(Σ − 6). **This is independent support for the normalisation V test 4
  chose.**
- **Not covered by Fang:** MN's closed-form A, B, C, the ω constant and the (1−α²)⁻² prefactor. These are still
  checked only against GLM and LGAC, and need MN 1992 itself.
- Honest wording: "agrees with Fang 2007 (his evaluation of MN 1992 Eqs. 11–13)". NOT "agrees with MN 1992".

## Prior art brought to our attention by an expert reader (2026-10-21), and Q4: the ring is MASSLESS (2026-10-23)
- **Prior art:** off-axis "massless ring singularities" in VACUUM stationary axisymmetric solutions are known:
  - Manko & Ruiz, Prog. Theor. Phys. 125 (2011) 1241 (arXiv:1101.4687), §IV.B, p. 1250: a massless equatorial ring
    in the Kinnersley–Chitre metric (two corotating extreme Kerr sources), called "a very unexpected result";
  - Manko, Prog. Theor. Phys. 127 (2012) 1057 (arXiv:1110.6564): δ = 2 Tomimatsu–Sato. Already cited in our query, for the
    geometry.
  - So the PHENOMENON is not new. What our V15 adds, as far as we know: in the Manko–Novikov metric as parametrised
    by GLM 2008, there is exactly one such equatorial ring for every q > 0, outside the horizon.
- **Q4, the Komar mass of the ring (pre-registered, 03b00a5; code/v15_q4_komar.py; results/q4_komar.txt):**
  - Setup: χ = 0.9, q = 0.95, M = 1; ring at x_s = 1.92982825852324290 (y = 0); mpmath at 30 digits.
  - C1 (asymptotic contour, x = 200): Komar mass 1.0 for both Kerr and MN. PASS, after a disclosed correction: as
    literally registered, C1 gave 2.0 for both. It caught a factor-2 error in my registered normalisation
    (−(1/8π)∮*dξ♭ instead of the standard Komar −(1/8π)∮*∇ξ = −(1/16π)∮*dξ♭). The corrected normalisation is the
    standard one and was fixed by the control, not by the ring value. The ring verdict below does not depend on the
    factor, because it is zero either way.
  - C2 (circle around a regular point): −1.6e−32. PASS.
  - Ring circles r = 1e−2, 3e−3, 1e−3: 7.0e−30, 2.3e−29, 7.0e−29 (numerical noise ∝ 1/r, as the integrand grows near
    the ring).
  - **VERDICT (registered rule |M_ring| ≤ 1e−10 at every radius): the MN equatorial ring is MASSLESS.** It is of
    the same kind as the Manko–Ruiz / Tomimatsu–Sato massless ring singularities.
- **Quadrupole sign (checking a suggested explanation, that the singularity might be due to a positive mass
  quadrupole moment):** our axis-Ernst multipoles (FINDINGS §5 above, matching GLM Eq. 4) give a mass quadrupole
  −M(a² + qM²). So q > 0 makes it MORE negative than Kerr's −Ma², and the equatorial ring exists exactly for q > 0.
  In our parametrisation the ring is therefore NOT tied to a positive quadrupole. (A positive quadrupole needs
  q < −χ², the regime of the off-equator ring pair, whose curvature is unexamined.) This is recorded as data.
