# V15: check of Fang 2007 (Caltech PhD thesis) against our Manko–Novikov transcription

Date: 2026-10-18. Literature check only. Sources were read, and nothing was contacted or sent.

## 1. The thesis: found, open access

- **Citation in GLM 2008** (arXiv:0708.0628, LaTeX source `exactdynamics.tex`): `\bibitem{HuaThesis}`, rendered as "Fang, H., Ph.D. thesis, California Institute of Technology, 2007." GLM cite it twice: once for the Eq. (14) typo, and once for the relations ρ = k(x²−1)^{1/2}(1−y²)^{1/2}, z = kxy.
- **Author:** Hua Fang. The advisor is Kip S. Thorne, per the Mathematics Genealogy Project id 141043.
- **Title:** *Topics in Gravitational Physics: Tidal Coupling in Gravitational Wave Searches and Mach's Principle* (Caltech, 2007).
- **DOI:** 10.7907/D99H-J577 (DataCite record). Resolver: https://resolver.caltech.edu/CaltechETD:etd-05212007-004257, which redirects to https://thesis.caltech.edu/1915/.
- **PDF:** https://thesis.caltech.edu/1915/01/thesis.pdf (1,937,126 bytes).
- **Access:** open. The record page says CaltechTHESIS "is an open-access repository". The PDF downloaded with no login. (`thesis.library.caltech.edu` refused connections from this box, but the `thesis.caltech.edu` mirror served the files.)
- **Relevant part:** Chapter 3, *Gravitational Waves from Intermediate-mass-ratio Inspirals for Ground-based Detectors*, Appendix 3.A, *SARSAF spacetimes that differ from Kerr*. The MN material is in §3.A.1, *Exact vacuum solutions: The Manko-Novikov Spacetime*, printed pp. 74–78 (PDF pp. 85–89). The abstract does not mention MN, so a title or abstract search will not find it.

## 2. What Fang reproduces, and what he does not

Page numbers are printed page numbers. Fang uses βₙ for MN's αₙ, which matches GLM's identification "α₂ = β".

| Fang eq. (page) | Content | MN 1992 eq. cited by Fang |
|---|---|---|
| (3.54) p.74 | Papapetrou line element in prolate spheroidal (x,y). Same form as GLM Eq. (2). | — |
| (3.55) p.74 | ρ = k(x²−1)^{1/2}(1−y²)^{1/2}, z = kxy | — |
| (3.56)–(3.58) p.75 | Ernst potential and Ernst equation | Eq. (4) of MN is cited for γ |
| (3.60)–(3.63) pp.75–76 | Schwarzschild and Kerr (Ernst form). sec λ = M/k, tan λ = a/k, k = (M²−a²)^{1/2}, r = kx + M, cos θ = y | — |
| (3.64a,b) p.76 | **MN Ernst potential:** E = e^{2ψ} A₋/A₊, A∓ = x(1+ab) + iy(b−a) ∓ (1−ia)(1−ib) | MN "solution generation scheme" (no eq. number given) |
| (3.65) p.76 | ψ = Σ_{n≥1} βₙ R^{−n−1} Pₙ(xy/R), R = (x²+y²−1)^{1/2} | — |
| (3.66) p.77 | Multipoles M₀…M₃, J₀…J₃ (Fang's own computation, using Quevedo's recurrence) | MN Eq. (14) (via footnote 20) |
| (3.67) p.77 | ψ₂ = β₂ R^{−3} P₂(xy/R) | — |
| (3.68a,b) p.77 | ψ₀ = ψ_S + ψ₂, and an **explicit closed form for γ₀** (the q-subclass case) | "auxiliary potentials defined in Eq. (11) and (12) of [MN]" |
| (3.69), (3.70) p.77 | a(x,y), b(x,y) | "according to Eq (13) of Ref. [MN]" |
| (3.71a–c) p.78 | Large-r expansions of f, ω, e^{2γ} in Boyer–Lindquist coordinates, to O(r⁻⁶), O(r⁻⁴), O(r⁻⁶) | from "expression (9) of [MN]" |
| (3.73)–(3.75) p.78 | ΔM₂, ΔM₂/M₂,K, ΔJ₃/J₃,K | — |

**Not reproduced:** the closed-form f, ω, e^{2γ} (MN Eq. 9), and the functions A, B, C. Fang only says that (3.68b), (3.69) and (3.70) "can then be inserted into expression (9)" of MN (p.77). He gives no quotation of MN's γ′ general formula, and no ω constant.

## 3. The MN Eq. (14) typo (Fang footnote 20, p.77)

Fang writes: "There seems to be a typo in the expression for M₂ in Ref. [32]" (Ref. [32] is MN 1992). He goes on: "the last α₂ should be α²". The typesetting was checked from PDF word bounding boxes. In the first symbol, the "2" sits below the baseline, so it is a subscript, α₂. In the second, the "2" sits above, so it is a superscript, α².

Corrected (Fang 3.66 = GLM Eq. 4): M₂ = −k³[β + 4α²(1+α²)(1−α²)^{−3}].

So, as far as Fang reports, MN 1992 printed an α₂ (the Weyl-quadrupole coefficient) where α² (the rotation parameter squared) belongs. That is the last α in this expression, i.e. the one inside (1−α²)^{−3}. **Inference, not verified:** MN 1992 itself was not seen. "The last α" is Fang's wording, and which factor it means is my reading of it.

The typo affects only the moment formula. It does not touch the metric functions.

## 4. Symbol-by-symbol comparison

Our transcription is `mn_functions` in `falsification/V9_mn_obstruction_check/code/v9_nve.py`, which I read only. GLM refers to the LaTeX source, Eqs. (2)–(5).

| Item | Fang 2007 | GLM 2008 | Ours (`mn_functions`) | Verdict |
|---|---|---|---|---|
| line element | (3.54) | Eq. (2), identical | used downstream (`axial_r`) | agree |
| R, Pₙ argument | (x²+y²−1)^{1/2}, Pₙ(xy/R) | same | `x*sqrt(1+(y²−1)/x²)`, u = xy/R | agree for x > 0 |
| ψ | β₂R^{−3}P₂ (3.67) | βR^{−3}P₂ | `beta*P[2]/R**3` | agree |
| a | −α exp(−2β₂(−1 + Σ_{ℓ=0}^{2}(x−y)P_ℓ/R^{ℓ+1})) (3.69) | identical | identical | agree |
| b | α exp(2β₂(1 + Σ(−1)^{**1−ℓ**}(x+y)P_ℓ/R^{ℓ+1})) (3.70) | (−1)^{**3−ℓ**} | (−1)^{3−ℓ} | agree: (−1)^{1−ℓ} = (−1)^{3−ℓ}, a notational difference only |
| f | Re E with (3.64) | e^{2ψ}A/B | e^{2ψ}A/B | agree: we checked that \|A₊\|² = B exactly, and that f matches Re E at 3 points to 40 digits |
| ω | not given in closed form. Large-r form (3.71b): ω = ω_K[1 − β₂k³(5cos²θ−1)/(Mr²) − β₂k³(47cos²θ−11)/(4r³) + O(r⁻⁴)] | 2k e^{−2ψ}C/A − 4kα/(1−α²) | identical to GLM | consistent: at M=1, a=0.9, β=0.3, the relative correction for our ω is 3.46942e-7 vs Fang 3.46885e-7 at r=200 (cosθ=0.3), and −1.248061e-7 vs −1.248032e-7 at r=400 (cosθ=−0.6). Residuals scale as O(r⁻⁴) |
| A, B, C | not given | Eq. (3) | identical to GLM | agree with GLM; nothing in Fang to compare against |
| γ′ / γ₀, log term | ½ log((x²−1)/(x²−y²)) | same | prefactor (x²−1)/(x²−y²) | agree |
| γ′, β² term | (3β₂²/8R¹⁰)[x²y²(5x²y²−3R²)²/R² − (3x²y²−R²)²] | (9β²/6R⁶)(P₃²−P₂²) | exponent 3β²(P₃²−P₂²)/R⁶ (= 2γ′) | agree: sympy shows Fang's form = (3β²/2R⁶)(P₃²−P₂²) exactly |
| **γ′, linear term (the −2 question)** | **β₂x(2x⁴−5x²+5x²y²+3−3y²)/R⁵ − 2β₂** | printed as βΣ_{ℓ=0}^{2}(…P_ℓ/R^{ℓ+1} − 2), with the −2 inside the sum (a total of −6β) | 2βΣ(… − 2) + 8β = 2β(Σ − 2) (the "outside" reading) | **Fang agrees with ours, not with GLM's literal print.** sympy: Fang − β(Σ−2) = 0 exactly, and Fang − β(Σ−6) = 4β |
| e^{2γ} prefactor | not given (sent to MN Eq. 9) | e^{2γ′}A(x²−1)^{−1}(1−α²)^{−2} | same | agree with GLM; Fang cannot test this |
| α, k | M₀ = k(1+α²)/(1−α²), J₁ = −2αk²(1+α²)/(1−α²)². Kerr: k = (M²−a²)^{1/2} | Eq. (5): α = (−1+√(1−χ²))/χ, k = M(1−α²)/(1+α²) | α the same, k = √(M²−a²) | agree, since (1−α²)/(1+α²) = √(1−χ²) |
| β | M₂ = −k³[β₂ + …] (3.66) | β = qM³/k³ | `beta` input directly | agree |

**Combined numeric check.** I evaluated e^{2γ} built from Fang's γ₀ (3.68b), Fang's a and b (3.69–3.70), and GLM's A and (1−α²)^{−2} prefactor. It equals our `e2g` to 40 significant digits at (x,y) = (3/2,1/3), (5/2,−7/10) and (4,1/10), for M=1, a=0.9, β=0.3. Our f equals Re of Fang's Ernst potential at the same points to 40 digits. On the axis (y=1), Fang's linear term is 2β₂ − 2β₂ = 0. That is the regularity condition our reading satisfies, and the literal "inside" GLM print does not satisfy it.

**Two small internal points in Fang, not affecting us:**
- (3.73) writes ΔM₂ = β₂(M²−a²)^{3/2} with no minus sign, while (3.66) gives ΔM₂ = −k³β₂. This is a sign slip, or an unstated sign convention. GLM's q = −ΔM₂/M³ is consistent with (3.66).
- In the text extraction, (3.66) appears to give M₁ = −k²β₁ and J₂ = 4αβ₁k³/(1−α²). These matter only for β₁ ≠ 0, and we have β₁ = 0.

## 5. What we can honestly claim

Fang paraphrases MN 1992 and re-derives from it. He does not quote MN's metric formulae verbatim. He prints MN's Ernst-potential form, and he applies MN Eqs. (11)–(13) to the α₂-only case. So the following statement is **not** supported: "our transcription agrees with MN 1992 as quoted by Fang 2007". Fang never prints MN's f, ω, γ, A, B, C.

Supportable statements:

1. Our a, b and ψ agree with Fang 2007, Eqs. (3.67), (3.69) and (3.70). Fang presents these as MN 1992 Eq. (13) specialised to α₂ = β.
2. Our γ′ agrees exactly, including the "−2 outside the sum" normalisation, with the explicit closed form in Fang 2007 Eq. (3.68b). Fang presents that form as his evaluation of MN 1992 Eq. (12). This independently supports the axis-regular reading over the literal GLM 2008 print.
3. Our f agrees with Fang's MN Ernst potential, Eq. (3.64). Our ω and f agree with Fang's large-r expansions, Eq. (3.71), at leading orders.
4. Our ω constant, A, B, C and the (1−α²)^{−2} prefactor agree with GLM 2008 Eq. (3) only. Neither GLM nor Fang lets us check these against MN 1992 itself. That still needs the MN 1992 paper, CQG 9, 2477.
5. The MN Eq. (14) typo, per Fang footnote 20: in M₂, an α₂ was printed where α² belongs. It concerns only the multipole formula.

Environment: Python 3.14.5 (`<venv>/`), sympy 1.14.0, mpmath 1.3.0. The check script and the PDF text extraction are in the session scratchpad and were not added to the repo. PDF text came from pdftotext (Homebrew poppler).
