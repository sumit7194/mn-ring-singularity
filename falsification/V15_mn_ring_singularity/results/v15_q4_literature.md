# V15 Q4: literature on an off-horizon singularity in Manko–Novikov

*Written 2026-10-17 by a Q4 subagent of the Bridge session. I did not read the Q1–Q3 outputs.*

*Access rules: full text only where it was legitimately reachable (arXiv, an open repository, INSPIRE metadata). I did
not circumvent any paywall or captcha.*

*Quotes are at most 15 words, one per paper. Section, figure and equation numbers are the arXiv versions'.*

**Status tags:**
- **[V]**: verified, read in full text, with the figure checked on the rendered page where relevant.
- **[P]**: partially verified (a plot reading, or a keyword screen only).
- **[U]**: unverified.

**The feature in question:** B = A = 0 on the equator, off the horizon rod. At C0 (χ 0.9, q 0.95):
- x_s = 1.929828;
- Weyl ρ = k√(x_s² − 1) = **0.7195 M**;
- Brink's prolate coordinate arccosh x_s = **1.2755**.

At F5-R: ρ = 0.7192 M (arccosh 1.2707). At (0.6, 0.1024): ρ = 0.318 M.

## Verdict

**No paper I could read identifies a curvature singularity off the horizon in MN.** This covers MN papers and
Manko-group follow-ups. **[V]** for the closely read papers below and **[P]** for the keyword screen.

- **The opposite is stated** (no singularity off the horizon):
  - **Brink 2008, Fig. 10**: a general schematic statement (resolved below). [V]
  - **Geyer 2013** (MSc thesis) §6.1. [V]
  - Implicitly, the many papers that place the only singularity at x = 1, y = 0. [V]/[P]
- **Compatible but not recognised as a singularity** (in every case the feature is plotted or described without being
  named):
  - GLM 2008: Fig. 1 and §II.B (the CTC region touches the bound region at a point at ρ/M ≈ 0.72), and a §III remark
    about singular metric behaviour near the CTC region;
  - Brink 2008 Fig. 5 (concentric contours on the equator at prolate ρ ≈ 1.2–1.3);
  - LGAC 2010 (a deep V_eff hole);
  - CLA 2011 and CHL 2012 (an inner boundary at ρ ≃ 0.72, independent of E and L_z);
  - Lukes-Gerakopoulos & Contopoulos 2013 (an extremely low V_eff minimum);
  - Geyer 2013 Fig. 6.2 (a pole-like blob in g_φφ);
  - Johannsen 2013 (a singular surface at Kerr's r₊; he states he did not compute Kretschmann).
- **Closest precedents (not MN):**
  - **Kodama & Hikida 2003 §III.A.3**: for the Tomimatsu–Sato δ = 2 solution, Ψ₂ ∝ B⁻³ and Ψ₀Ψ₄ ∝ B⁻⁶ diverge at
    the unique equatorial B = 0 point, where the inner ergo-boundary meets the equator. [V]
  - **Manko 2011 (PTP 127 1057)**: gives the general rule that in vacuum a ring singularity sits where the stationary
    limit surface touches the boundary of the CTC region. GLM's own description of MN at q > 0 matches this geometry,
    but nobody has applied the rule to MN. [V]
- **MN 1992 itself: NOT READ.** It is closed access (details below). So "not identified in MN 1992" is **[U]**.

## 1. Manko & Novikov 1992 and Manko follow-ups

**Manko & Novikov 1992, CQG 9, 2477** (DOI 10.1088/0264-9381/9/11/013). **Not read.** What I tried:
- IOP PDF: blocked by a Radware bot-manager captcha. Not circumvented.
- INSPIRE record 2932082 (created 2025): abstract only, no arXiv or fulltext link.
- OpenAlex: `is_oa: false`, `oa_status: closed`, no repository full text.
- Crossref: a similarity-checking link only.
- Copenhagen research portal page: metadata only, no file.
- ADS abstract page: WebFetch returned 405.
- Web search: no legitimate copy found.

The INSPIRE abstract describes the construction and the multipoles only and says nothing about singularities. [V]
Second-hand: Johannsen 2013 (arXiv:1304.7786 §V) cites MN 1992 (his ref. [34]) for the statement that the MN metric harbors a naked
singularity. What MN 1992 says about where that singularity is, and whether it gives a regularity or validity
statement, is **unknown**. [U]

**Other Manko papers on this family, not read (closed, no arXiv):**
- Castejon-Amenedo & Manko, PRD 41 2018 (1990);
- Manko, CQG 7 L209 (1990);
- Manko et al., PRD 57 3382 (1998).

**Manko-group papers read in full text:**

| Paper | What it says about singularities off the axis / horizon | Bearing on the MN ring |
|---|---|---|
| Manko, Mielke, Sanabria-Gómez 2000, gr-qc/0001081 (PRD 61 081501) | Locates the ring singularities of their solution on the equator between the axis and the ergosphere. Quote: "no singularity outside the infinite red shift surface arises". | Shows the group's standard practice: locate the off-axis rings relative to the stationary limit surface (SLS). Does not concern MN. [V] |
| Manko 2011, arXiv:1110.6564 (PTP 127 1057), on the TS2 interpretation | Defines the vacuum ring singularity as "a locus of points at which the SLS touches the boundary of a CTC region". Attributes it to negative-mass constituents. | **Closest Manko-side precedent.** GLM §II.A and Fig. 1 describe exactly this geometry for MN at q > 0: the outermost edge of the CTC region lies inside the ergoregion and touches it on the equator. MN is not mentioned. [V] |
| Manko & Ruiz 2011, arXiv:1110.6764 (PTEP 2013 103E01) | A massless ring singularity off the axis in Kerr–Newman with negative mass and in charged TS2. | Same mechanism. MN is not mentioned. [V, keyword plus targeted read] |
| Manko, Mejía, Ramos, Ruiz 2023, arXiv:2304.00067 (PRD 108 044042) | Off-axis ring singularities are roots of the Ernst denominator and lie "on the stationary limit surface (SLS) f = 0 in the pure vacuum case". Equatorial ones are found first. | Same structure as the MN ring, where B = 0 is the Ernst denominator and A = 0 is f = 0. MN is not mentioned. [V] |
| Manko et al. 2011, arXiv:1105.2646 (PRD 84 064005) | Massless ring singularities off the axis in Kinnersley–Chitre binaries. | Cites MN 1992 only for an ω formula. [V] |
| arXiv:1603.08054, 1909.09952, 1707.04953 | No MN singularity statement. MN is not cited in their text. | — [P] |

## 2. Gair, Li & Mandel 2008 (arXiv:0708.0628, PRD 77 024035): read in full [V]

**§II.A (Spacetime properties):**
- They say the spacetimes must either lack an event horizon or have CTCs outside one, and that both hold.
- A partial horizon at ρ = 0, |z| ≤ k is broken on the equator by a circular line singularity at x = 1, y = 0. This is
  credited to Brink by private communication (ref. [20]).
- Quote: "For χ = 0 the spacetime is otherwise regular". No regularity claim is made for χ ≠ 0.
- For χ ≠ 0 they describe an ergosphere (g_tt = 0) and a CTC region (g_φφ < 0).
- They call the structure very similar to TS δ = 2, citing **Kodama & Hikida** (ref. [21]). That is the paper which
  derives TS2's equatorial B = 0 ring singularity. GLM did not carry it over.
- For q > 0, the outermost edge of the CTC region lies inside the ergoregion.
- They note that CTC zones are sometimes regarded as unphysical, and suggest excising them with an inner boundary.

**Fig. 1 (χ = 0.9, q = −1, 0, 1)** is the "bad region" figure (our term): zeros of g_tt and g_φφ. Readings off the rendered
q = 1 panels [P, plot reading ±0.02]:
- the inner tip of the unshaded g_tt < 0 slit inside the equatorial ergo-lobe is at ρ/M ≈ 0.72–0.74;
- the outer tip of the equatorial CTC lobe is at the same ρ/M ≈ 0.72–0.74;
- the two boundaries meet at one equatorial point, which is where an A = 0 / B = 0 ring would sit;
- a second, inner unshaded slit inside the CTC lobe ends at ρ/M ≈ 0.57.

Nothing is marked singular in the figure.

**§II.B, at χ 0.9, q 0.95, E 0.95, L_z 3M:**
- The inner bound region is ρ/M ≈ 0.72–2.12.
- In the plane z = 0, the outermost edge of the CTC region touches the innermost edge of the bound region.

0.72 is the ring's ρ (0.7195) to the stated precision. **Compatible.**

**§III:**
- The ergodic region touches the CTC region at a single point.
- Singular metric behaviour as the CTC region is approached might explain the ergodicity.

This is speculation, not an identification of a curvature singularity.

**Parameter ranges:**
- χ = 0.9 throughout the structure plots, with q ∈ {−1, 0, 1} in Fig. 1;
- the orbit study is at q = 0.95;
- both regions grow as |q| increases.

**Overall:** GLM do not identify the ring, but they describe and plot its location (the CTC/bound-region touching point
at ρ ≈ 0.72). **Compatible.**

## 3. Brink 2008, "Spacetime Encodings II" (arXiv:0807.1179, PRD 78 102002): read in full, figures rendered [V]

**What "No singularities outside" refers to.**
- **Fig. 10** (p. 10) is a cartoon, captioned as the properties of the MN spacetime. The text annotations inside the
  image (paraphrased except the last) are:
  - vacuum, R_ij = 0;
  - any mass multipoles and spin can be specified;
  - the horizon area is fixed by the multipole moments;
  - an "Event Horizon" pointer;
  - an arrow to the horizon's equator labelled "**Curvature Singularity on Horizon Equator (No singularities
    outside)**".
- These labels are embedded in the image, so pdftotext misses them. I checked them on the page rendered at 110 dpi.
- It is a **general, qualitative claim about the MN family**. No parameter point, derivation, curvature computation or
  supporting sentence is given.
- The text (§VII, p. 6) says only that some properties are sketched in Fig. 10.
- App. C says the same metric generated Figs. 5 and 6 and has its properties sketched in Fig. 10. So the cartoon is at most tied to the App. C quadrupolar solution and, by implication, the Fig. 5–6 parameters.

**Parameters:**
- The Fig. 6 caption gives α = 0.626789 and α₂ = 11.4708. These are |α| and β for **χ = 0.9, q = 0.95** (GLM sign
  convention α = −0.626789). [V, arithmetic]
- The Fig. 5 caption gives no parameters. [V]

**Brink's own Fig. 5 shows the feature.** Brink's ρ is prolate (x = cosh ρ, y = cos z; App. C), so the equator is
z = π/2. Readings [P, plot reading ±0.05; the parameters are assumed to be Fig. 6's]:
- the g_tt panel has a set of concentric contours centred on the equator at prolate ρ ≈ 1.22–1.27, where several zero
  contours converge;
- the g_φφ panel shows the CTC lobes meeting the equator in a cusp at ρ ≈ 1.25–1.3;
- for comparison, arccosh x_s = 1.2755.

The figure labels only the ergoregions and the closed timelike curves.

**Other notes:**
- The App. C metric has the same A, B and C as GLM [V]. As already noted in our earlier K20 notes, his second ln(b/α)
  line prints P₂/r² for P₂/r³.
- §VII says the chaotic region A lies outside the Fig. 5 ergoregions.

**Resolution:** the "no singularities outside" label is a blanket schematic assertion for the MN family. It sits beside a
plot at χ 0.9, q 0.95 that shows the ring as unlabelled concentric contours. It **states the opposite** of the V13 /
V15 claim, with no computation behind it.

## 4. Kodama & Hikida 2003 (gr-qc/0304064, CQG 20 5121), §III.A.3 "Singularity": read in full [V]

**Setup:**
- For TS2, f = A/B and R² ∝ D/B (Eq. 33), with explicit polynomials A, B, C, D (Eq. 34) and g(x) = p²x⁴ + 2px³ − 2px − 1.
- §III.A.1: the ergoregion in x ≥ 1 is x₀² − (x₀² − 1)y² < x² < x₁² − (x₁² − 1)y², where x₀ is the positive root of
  g(x). So x₀ is the inner ergo-boundary on the equator.
- §III.A.2: the CTC region (D < 0) contains the rod and touches the inner ergo-boundary (Fig. 7).

**Method (§III.A.3):**
1. For the tetrad adapted to dt − ω dφ, Ψ₁ = Ψ₃ = 0, so every polynomial invariant is a polynomial in Ψ₂ and Ψ₀Ψ₄.
2. Symbolic computation gives Ψ₂ = P₁/B³ and Ψ₀Ψ₄ = P₂/B⁶ (Eq. 38), with P₁ and P₂ coprime to B.
3. The equatorial closed forms (Eq. 39) carry g(x)³ and g(x)⁶ in their denominators.
4. Quote: "the curvature invariants diverge at the zero of g(x), x = x0".
5. B vanishes for x ≥ 1 **only** at (x₀, 0), apart from the poles P±.
6. x₀ is a single root of D and a double root of B, so R² = D/B → ∞: **a ring singularity of infinite circumference**
   (KH call it the famous ring singularity of Tomimatsu–Sato 1972/73).

**Precedent value:**
- It is the same structure claimed for MN: B = 0 at one equatorial point off the horizon rod, A = 0 there (inner
  ergo-boundary), the CTC region touching, and g_φφ with a pole, not a zero.
- It gives the uncancelled-B argument that Q2 needs.
- GLM cite this paper for the TS2–MN similarity. No one has carried the result over to MN.

## 5. Later MN literature

**Method:**
- **Primary set:** the INSPIRE citation lists of MN 1992 (recid 2932082, 142 citing) and GLM 2008 (recid 757591,
  170 citing). Together that is **252 unique papers**.
  - **216** have arXiv full text: downloaded and **full-text keyword screened**.
  - **36** without arXiv were screened on **metadata/abstract only**: books, theses, proceedings, some journal-only
    papers.
- **Supplementary set:** INSPIRE `fulltext:"Manko-Novikov"`, abstract search, and the citers of Brink 2008 (recid
  790157, 43) and LGAC 2010 (recid 849051, 101). That adds **72 more**: 68 full-text screened, 4 metadata only.
- The INSPIRE MN-1992 record was created in 2025 and may under-count citations. OpenAlex gives 156 against INSPIRE's
  142.
- **Totals:** 324 papers; 284 full-text screened, 40 metadata only. Plus 7 Manko papers, the Geyer thesis and
  Kodama–Hikida. The per-paper record is in `v15_q4_screen.json`.
- **Screen:** pdftotext, then a regex for "singular" or "naked" within a sentence of Manko/MN/bumpy, plus patterns for
  "singular … outside/off … horizon", "ring singularity", B = 0, Kretschmann/curvature invariant and 0.72.
  - In the primary set, 42 papers had singularity sentences near MN. All flagged hits were read in context.
  - About 15 papers were read closely: the figures and sections cited below.
- **Limitation:** text embedded in figures is invisible to the screen. Brink's Fig. 10 shows this can matter.

**Closely read papers:**

| Paper (full text unless marked) | Statement about singularities outside the horizon (≤ 15-word quote) | Ring identified? |
|---|---|---|
| Lukes-Gerakopoulos, Apostolatos, Contopoulos 2010, arXiv:1003.3120 (PRD 81 124005) §2 and §3.4 | Horizon broken on the equator by a circular line singularity. Quote: "outside some central region where singular behavior of the metric functions shows up". In the inner region, a deep hole in V_eff (to −10⁵) scatters orbits stochastically. | No. The deep hole and the vague "central region" are **compatible**. [V] |
| Contopoulos, Lukes-Gerakopoulos, Apostolatos 2011, arXiv:1108.5057 (IJBC 21 2261) §2–3 and Figs. 15–17 | Line singularity only, at ρ = z = 0. Quote: "the inner boundary remains at an almost constant position". At ρ ≈ 0.72 for both E 0.95 → 1 (L_z = 3) and L_z sweeps, and still 0.72 at E ≈ 0.19 (Fig. 17). | No. A CZV boundary that does not move with E or L_z points to a metric feature, not dynamics: **compatible**. [V] |
| Contopoulos, Harsoula, Lukes-Gerakopoulos 2012, arXiv:1203.1010 §4 and Fig. 16 | Quote: "The inner boundary is at an almost constant ρ ≃ 0.72." CTCs exist outside the horizon. | No. **Compatible**: ρ_ring = 0.7195. [V] |
| Lukes-Gerakopoulos & Contopoulos 2013, arXiv:1304.7612, Fig. 2 (χ 0.9, q 0.95, E 0.95, L_z 3) | Quote: "The local minimum between C′ and A is extremely low". The inner root C′ is at ρ ≈ 0.7 on the equatorial V_eff plot. | No. **Compatible.** [V text; P plot] |
| Bambi & Barausse 2011, arXiv:1108.4740 (PRD 84 084034) §II | Quote: "a naked curvature singularity on the equatorial plane (i.e. at x = 1, y = 0". CTCs at small radii. | No. Only the horizon singularity. [V] |
| Bambi & Lukes-Gerakopoulos 2013, arXiv:1302.0565 (PRD 87 083006) §II | Quote: "MN spacetimes with q ≠ 0 have regions with CTCs just outside the event horizon". Also gives the x = 1, y = 0 singularity. | No. [V] |
| Johannsen 2013, arXiv:1304.7786 (PRD 87 124010) §IV–V, Figs. 10–11 | In his Boyer–Lindquist-like rescaled form the metric is singular at Kerr's r₊. Quote: "I did not calculate its Kretschmann scalar to identify the nature of this singularity". Fig. 11 gives the outermost equatorial CTC radius against q. Thanks V. Manko for conversations. | No. The outermost equatorial CTC point may coincide with the ring for q > 0 (GLM Fig. 1), but that link is **[U]** pending Q3. [V] |
| Destounis & Kokkotas 2021, arXiv:2108.02782 (PRD 104 064023) §II, at exactly χ 0.9, q 0.95, E 0.95, L_z 3 | Horizon broken on the equator by a circular line singularity. Quote: "Regardless of its pathological nature, the MN metric is a prime candidate". | No. [V] |
| Destounis, Suvorov, Kokkotas 2020/21, arXiv:2009.00028 and 2103.05643 | No MN singularity statement. | No. [P] |
| Lukes-Gerakopoulos & Witzany 2021, arXiv:2103.06724 (review) | Quote: "this manifests as a ring singularity on the horizon". | No. Horizon only. [V] |
| Dey, Barausse, Basak 2022, arXiv:2212.10725 (QNMs / ringdown) | Quote: "it has (in general) a naked singularity on the equatorial plane (at x = 1, y = 0)". | No. Horizon only. [V] |
| Staelens, Mayerson et al. 2023, arXiv:2303.02111 (photon rings and imaging) §2.3 | Quote: "the event horizon has a singularity (only) at the equator". | No. **Implicitly the opposite**: the only singularity is placed on the horizon. [V] |
| Wang, Chen, Jing 2018, arXiv:1801.02118 (MN shadow) §II | Quote: "Although there exist the naked singularity and the closed timelike curves". | No. The singularity is not located. [V] |
| Wu, Zhang, Shao et al. 2023, arXiv:2303.04907 (MN lensing) | Quote: g_φφ > 0 holds "in accordance with most relevant axisymmetric metrics". | No. [V] |
| Yunes, Siemens et al. Living Review 2013 / 2025, arXiv:1304.3473 and 2408.05240 | Quote: "the event horizon is broken along the equator by a ring singularity". Their mention of a naked equatorial ring singularity refers to the Johannsen–Psaltis metric, not MN. | No. [V] |
| Johannsen–Psaltis / Johannsen 2011–2016 (1105.3191, 1105.5645, 1304.8106, 1501.02815, 1512.03818, …) | Generic: MN-type metrics have singularities or CTC regions outside the event horizon, with no location. | No. Generic only. [P] |
| Bambi et al. continuum-fitting papers (1012.2007, 1103.5135, 1112.4663, 1509.03884, 1511.07587) | Naked singularities and CTCs at small radii, said to be covered by an exotic object. | No. [P] |
| Geyer 2013, MSc thesis, Stellenbosch (hdl 10019.1/80306; open PDF) §6.1 and Fig. 6.2 | Quote: "an event horizon which has singularities only on the equator and is well behaved elsewhere". Fig. 6.2 (q 0.95; prolate coordinates) shows concentric g_tt contours and a dense g_φφ blob on the equator at ρ ≈ 1.28–1.35. | No. The text **states the opposite**; the figure is **compatible** and unflagged. [V text; P plot] |

Other MN-heavy papers screened with no singularity statement beyond the horizon: Apostolatos et al. 2009
(0906.0093), Lukes-Gerakopoulos et al. 2014 (1408.4697), Seyrich & Lukes-Gerakopoulos 2012 (1207.3175; its "almost
singular" refers to the integrator and the MSM metric), and ray-tracing / QPO / disc papers. [P]

**Unrecognised traces that line up with the ring:** the ring sits at ρ = 0.7195 M at C0. Five sources show something
at that radius or prolate coordinate, and none flags a singularity:
1. GLM: the 0.72 touching point between the CTC region and the inner bound region.
2. CLA/CHL: the 0.72 inner CZV boundary, which does not move with E or L_z.
3. LGAC: the deep V_eff hole.
4. Brink Fig. 5: the concentric g_tt contours at prolate 1.22–1.27.
5. Geyer Fig. 6.2: the g_φφ blob.

These are consistent with a metric singularity at the ring, but **that connection is our inference [U]**.

## What could not be read
- **Manko & Novikov 1992**, CQG 9 2477: closed. The IOP route is captcha-blocked, and there is no repository copy
  (OpenAlex says closed).
- **Castejon-Amenedo & Manko 1990** (PRD 41 2018), **Manko 1990** (CQG 7 L209) and **Manko et al. 1998** (PRD 57
  3382): closed, no arXiv.
- **40 non-arXiv citing works:** metadata or abstract only. Examples: Bambi's 2017 book "Black Holes: A Laboratory for
  Testing Strong Gravity", the LRR 2025 journal record (its arXiv version 2408.05240 was screened), several PhD theses and
  proceedings, and a 1994 NUT–Curzon paper. The list is in `v15_q4_screen.json` (mode = metadata/abstract only).
- **Text inside figures** of the 284 screened papers: not machine-screened. Only the figures discussed above were
  inspected visually.
