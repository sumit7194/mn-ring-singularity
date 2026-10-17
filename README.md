# Manko–Novikov spacetime: an equatorial ring singularity outside the horizon (curated snapshot)

This is a curated, read-only snapshot of the code, pre-registration and results behind one question.

**The question:** in the Manko–Novikov (MN) vacuum metric (as written by Gair, Li & Mandel 2008 and
Lukes-Gerakopoulos, Apostolatos & Contopoulos 2010), the functions B and A vanish together at one equatorial point
outside the horizon for every q > 0. Our computations say this is a curvature ring singularity (Kretschmann ∝ d⁻⁶).
We could not find it in the literature we could read, and we could not read Manko & Novikov 1992 itself, so we are
asking whether it is already known. **We do not claim it is new.**

The one-page summary is [`SUMMARY.md`](SUMMARY.md). The full write-up, with every registered criterion and its
outcome (including one that failed as registered), is
[`falsification/V15_mn_ring_singularity/FINDINGS.md`](falsification/V15_mn_ring_singularity/FINDINGS.md), and the
pre-registration is `PREREGISTRATION.md` in the same folder.

**Status:** computer-assisted, done with AI assistance (several Claude-based agents writing and cross-checking code
under the author's direction). Every computation was pre-registered, and failures are kept on record. No professional
physicist has reviewed it yet. Corrections are very welcome.

## What is here

| path | what it is |
|---|---|
| `falsification/V15_mn_ring_singularity/` | The study: `PREREGISTRATION.md`, `FINDINGS.md`, `code/`, `results/` |
| `.../code/v15_q1_existence.py` | Q1: closed-form equatorial condition, (χ, q) existence map, off-equator search |
| `.../code/v15_q2_curvature.py` | Q2: Kretschmann and Weyl invariants near the ring (2+2 block-Riemann route, 100 digits); `--compare-v13` compares with an earlier, separate implementation's table |
| `.../code/v15_q2b_check.py` | Q2b: the analytic derivation of the d⁻⁶ law (`results/v15_q2b_derivation.md`) checked against the numeric K·δ⁶ at four rings |
| `.../code/v15_q3_geometry.py` | Q3: ergosurface, g_φφ and B near the ring; the CTC region |
| `.../code/v15_V_*.py` | The adversarial verifier: an independent transcription (`v15_V_mn.py`), its own curvature code (`v15_V_curv.py`), and tests 2–6 (`t2_compare`, `t3_vacuum`, `t5_multipoles`, `t6_singularity`) |
| `.../results/` | Every result file (JSON with versions, logs, figures), the analytic notes `v15_q2_analytic.md`, `v15_q3_geometry.md`, the literature review `v15_q4_literature.md` and its screening record `v15_q4_screen.json`, the check against Fang's 2007 thesis `v15_fang2007_check.md`, and the verifier's report `v15_V_verifier.md` |
| `falsification/V9_mn_obstruction_check/code/v9_nve.py` | The earlier MN transcription that `v15_V_t2_compare.py` imports and compares against |
| `falsification/V13_mn_chaotic_sea/results/v13_kretschmann.json` | The earlier, separate implementation's Kretschmann table, read only by `v15_q2_curvature.py --compare-v13` |
| `ENVIRONMENT/` | `requirements.lock` (pip freeze) and `env.json` (interpreter, platform, package versions, BLAS/LAPACK) |
| `MANIFEST.sha256` | sha256 of every file, its source path and the source commit |

## Environment
CPython 3.14.5 on macOS 26.5 (arm64); sympy 1.14.0, mpmath 1.3.0, numpy 2.4.6, scipy 1.18.0. Each result JSON
records the versions it ran with. To recreate:

```
python3.14 -m venv .venv && .venv/bin/pip install -r ENVIRONMENT/requirements.lock
```

Only sympy, mpmath, numpy and scipy are used by this code; the other pinned packages are present in the interpreter
but not imported. Older or newer sympy/mpmath versions will very likely work, but these are the versions used.

## Re-running
Run from `falsification/V15_mn_ring_singularity/` (scripts find their inputs and outputs relative to their own
location, and **overwrite** the matching files in `results/`; work on a copy if you want to diff against ours).
Pin thread pools to one core first: `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`.

| command | what it reproduces | cost (single core, Apple M-series) |
|---|---|---|
| `python code/v15_q1_existence.py` | Q1: the closed form, the registered roots, the (χ, q) map, the off-equator analysis | ~4 min (original run 220 s) |
| `python code/v15_q2_curvature.py` then `python code/v15_q2_curvature.py --compare-v13` | Q2: K near the ring at C0, F5-R and the third ring; Kerr controls; slope fits; the comparison table | ~30 s (original run 32 s) |
| `python code/v15_q2b_check.py` | Q2b: closed-form coefficient vs numeric K·δ⁶ at C0, F5-R, the third ring and Kerr | ~40 s (original run 39 s) |
| `OMP_NUM_THREADS=2 python code/v15_q3_geometry.py` | Q3: the local forms, the pole of g_φφ, the CTC region figure (uses 2 worker processes) | ~2 min (original run 126 s) |
| `python code/v15_V_t2_compare.py` | V test 2: independent transcription vs `v9_nve.py` (f, ω, e^{2γ}, both readings of the GLM Eq. 3h "−2") | see `TESTRUN` below |
| `python code/v15_V_t3_vacuum.py` | V test 3: full 4-D Ricci at 17 points (100 digits); the slowest script | not recorded; the slowest (it was run in two parts) |
| `python code/v15_V_t5_multipoles.py` | V test 5: axis multipoles vs GLM Eq. 4, far-field M and J | see `TESTRUN` below |
| `python code/v15_V_t6_singularity.py` | V test 6: local slopes, √(K² + P²) at 36 angles, proper distance, geodesic reachability | not recorded |

Parameter points: **C0** is χ = 0.9, q = 0.95 (M = 1); **F5-R** is M = 89, a = 80, q = 19/20; the third ring is
(χ, q) = (0.6, 0.1024).

**TESTRUN.** Q1, V test 2 and V test 5 were re-run from a copy of this snapshot on one core before it was committed;
see the bottom of this file for the outcome.

## Provenance
- Files are copied from the author's private research repository; the source commit for each file is given in
  `MANIFEST.sha256`. Paths in the manifest are relative to that repository's root.
- **Five code files were edited for this snapshot**, only to replace absolute paths on the author's machine by paths
  relative to the script (`v15_V_t2_compare.py`, `v15_V_t3_vacuum.py`, `v15_V_t5_multipoles.py`,
  `v15_V_t6_singularity.py`, `v15_q2_curvature.py`). The manifest gives both the original and the snapshot sha256
  for these.
- **Light redactions (2026-10-18), recorded per file in the manifest:** absolute paths on the author's machine in two
  log files and in `PREREGISTRATION.md` were replaced by `./`, `<venv>/` or `<home>/`; two sentences in
  `PREREGISTRATION.md` and `FINDINGS.md` about the author's internal tooling and review budget were shortened; and
  `ENVIRONMENT/env.json` no longer records the BLAS/LAPACK backend. No result, number, criterion or outcome was
  changed. Every other file is byte-identical to its source.
- `results/v15_q1_map_arrays.npz` is an output of `v15_q1_existence.py` that is not tracked in the source repository;
  it is included for convenience.
- The documents mention other projects of the author ("V9", "V13", other session names). These are internal labels:
  V9 and V13 are earlier studies whose only files needed here are included.
- Commit hashes cited inside the copied documents refer to the source repository.

Author: Sumit Pal (independent researcher). License: MIT (see `LICENSE`).

## TESTRUN (2026-10-17, before the snapshot was committed)
A copy of this snapshot, at a different absolute location, was run on one core (thread pools pinned to 1) with the
interpreter in `ENVIRONMENT/`:
- `v15_q1_existence.py`: 140 s. `Q1_PASS: true` (C0, F5-R, (0.6, 0.1024) and both Kerr negative controls). The output
  JSON equals the committed one in every field except `runtime_s`; `v15_q1_map_arrays.npz` and `v15_q1_map.svg` are
  identical.
- `v15_V_t2_compare.py`: 82 s. Every computed field equals the committed JSON. (The committed JSON also carries a
  `symbolic_note` key that the verifier added by hand after its run, explaining the symbolic residual; it is not
  written by the script.)
- `v15_V_t5_multipoles.py`: 4 s. Identical JSON.
