# Quantum rigidity of prime-order circulants via two-basepoint Terwilliger algebras

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22182428.svg)](https://doi.org/10.5281/zenodo.22182428)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Verification scripts and data for the paper

> M. F. Marashdeh, *Quantum rigidity of prime-order circulants via two-basepoint
> Terwilliger algebras*.

Everything here runs in exact arithmetic. **No floating-point arithmetic is used
anywhere.**

## What is verified

`verify.py` has four independent parts, each keyed to a numbered result in the paper.

| Part | Paper | What it checks |
|------|-------|----------------|
| **A** | Section 7 | The four capture certificates, for `p = 13, 17, 31, 41`. Recomputes each cyclotomic class, each block `B` and `D`, and the convolution profile, and confirms the isolated point mass. Also confirms that none of the four pairs has a singleton block, so that each certificate is of optimal depth. |
| **B** | Lemma 6.1, Theorem 6.2, Remarks 6.3–6.4, Proposition 7.3 | The confinement lemma (every vertex in a block of size ≥ 2 has the form `(1-e')/(e-e')` with `e, e' ≠ 1`) and the resulting bound `(k-1)(k-2)`, on all 173 pairs with `p < 200`; that the bound is attained, so it cannot be improved by this argument; the quadratic threshold `p > (k-1)(k-2) + 2`; the enumeration of the eight residual pairs in types ≤ 10; and their four singleton-block witnesses. |
| **C** | Theorem 8.1 | The saturation of the two-basepoint module for all **214** pairs `(p, E)` with `p ≤ 250` — over **Q** in exact rational arithmetic for `p ≤ 60`, and over **F_q** with `q = 1000003` throughout. Confirms `dim W(0,1) = p` in every case. |
| **D** | Definition 5.14, Section 8.2, Appendix A | The capture depths, **recomputed from scratch** in exact integer arithmetic for all 214 pairs and then compared with `depths.json`: the depth distribution; the implications `d ≥ 2 ⇒ r ≤ 7`, `d ≥ 3 ⇒ r ≤ 4`, `d ≥ 4 ⇒ r = 2`; the 48 pairs with `r ≥ 3` and `d ≥ 2`; the non-decreasing Paley depth sequence; the β-stability of every capture set (Lemma 5.15); and the harmonic observations of Remark 8.4. |

### Two points of principle

**Why the modular arithmetic in Part C is a proof.** The rank over **Q** of an integer
matrix is at least its rank over **F_q**, since a nonvanishing `p × p` minor modulo `q`
is a nonzero integer. The failure mode is therefore one-sided: the computation can
report failure for a family that does span, but never success for one that does not.

**Why Part D does not use it.** A capture is a *membership* assertion, not a spanning
one, and the argument above does not apply to membership: a vector can lie in the
reduction of a lattice modulo `q` without lying in its rational span. Part D is
therefore carried out in exact integer arithmetic, by echelon reduction over **Z** with
primitive rows.

### The capture depth

`V_0 = span{δ_0, δ_1}` and `V_d = Δ · B · V_(d-1)`, where `B` is the Bose–Mesner algebra
of the cyclotomic scheme and `Δ` is the algebra of diagonal matrices constant on the
blocks of the basepoint partition. The capture depth `d(p,E)` is the least `d` for which
`V_d` contains a point mass `δ_z` with `z ∉ {0,1}`. By Lemma 5.13 of the paper,
`V_1` is the span of the block indicators, so `d = 1` holds exactly when some cyclotomic
number equals 1 — that is, exactly when the depth-one criterion applies.

## Running it

```bash
pip install -r requirements.txt

python3 verify.py A      # certificates            (< 1 second)
python3 verify.py B      # threshold and 8 pairs   (about 10 seconds)
python3 verify.py C      # saturation, 214 pairs   (about 5 minutes)
python3 verify.py D      # capture depths          (about 3 minutes)
python3 verify.py all    # everything
```

Each part prints one line per check, ending with `OK` or `FAIL`, and the run ends with
`ALL VERIFICATIONS PASSED` or `*** FAILURE ***`. Parts C and D are the slow ones: in
Part C the modular sweep over the full range takes about seventy seconds and the
rational sweep over `p ≤ 60` about four minutes, on one core.

## Regenerating the data and the table

`depths.json` and Table 1 of Appendix A are both reproducible from scratch:

```bash
python3 compute_depths.py depths.json     # about three minutes
python3 make_table.py > appendix_table.tex
```

`make_table.py` emits the LaTeX `tabular` environment that appears verbatim in
Appendix A of the paper.

## Files

| File | Contents |
|------|----------|
| `verify.py` | The four verification parts described above. |
| `compute_depths.py` | Regenerates `depths.json` from scratch. |
| `make_table.py` | Regenerates the LaTeX source of Table 1 from `depths.json`. |
| `depths.json` | 214 rows `[p, k, d, z]`: the type, the capture depth, and the least captured vertex. |
| `requirements.txt` | `sympy`, `numpy`. |
| `CITATION.cff` | Machine-readable citation metadata. |
| `.zenodo.json` | Metadata for the Zenodo deposit. |

## Citing

Please cite the paper. If you use the code or data directly, cite this archive as well
via the DOI badge above; the DOI shown is the concept DOI and resolves to the latest
archived version.

## License

MIT — see [LICENSE](LICENSE).
