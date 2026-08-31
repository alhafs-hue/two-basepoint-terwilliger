# Two-basepoint Terwilliger algebras and the quantum rigidity of circulant graphs of prime order

<!-- After the first Zenodo release, replace the two placeholders below with the
     concept DOI (the "all versions" DOI) and delete this comment. -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Verification scripts and data for the paper

> M. F. Marashdeh, *Two-basepoint Terwilliger algebras and the quantum rigidity of
> circulant graphs of prime order*.

Everything here runs in exact arithmetic. **No floating-point arithmetic is used
anywhere.**

## What is verified

`verify.py` has four independent parts, each keyed to a numbered result in the paper.

| Part | Paper | What it checks |
|------|-------|----------------|
| **A** | Section 6.2 | The four capture certificates, for `p = 13, 17, 31, 41`. Recomputes each cyclotomic class, each block `B` and `D`, and the convolution profile, and confirms the isolated point mass. |
| **B** | Lemma 6.1, Theorem 6.2, Proposition 6.4 | The confinement lemma (every vertex in a block of size ≥ 2 has the form `(1-e')/(e-e')`) and the resulting bound `(k-1)^2`, on all 173 pairs with `p < 200`; the quadratic threshold `p > (k-1)^2 + 2`; the enumeration of the eight residual pairs in types ≤ 10; and their four singleton-block witnesses. |
| **C** | Theorem 7.1 | The saturation of the two-basepoint module for all **214** pairs `(p, E)` with `p ≤ 250` — over **Q** in exact rational arithmetic for `p ≤ 60`, and over **F_q** with `q = 1000003` throughout. Confirms `dim W(0,1) = p` in every case. |
| **D** | Appendix A | The capture-depth table: 214 rows, the quoted depth distribution, and the observation that for every Paley graph `P_p` with `p ≤ 241` the first captured point mass is `δ_2` or `δ_{1/2}`. |

The modular arithmetic in part C is sound because the rank over **Q** of an integer
matrix is at least its rank over **F_q**: a nonvanishing `p × p` minor modulo `q` is a
nonzero integer.

## Running it

```bash
pip install -r requirements.txt

python3 verify.py A      # certificates            (< 1 second)
python3 verify.py B      # threshold and 8 pairs   (about 5 seconds)
python3 verify.py C      # saturation, 214 pairs   (about 5 minutes)
python3 verify.py D      # appendix table          (< 1 second)
python3 verify.py all    # everything
```

Each part prints one line per check, ending with `OK` or `FAIL`, and the run ends with
`ALL VERIFICATIONS PASSED` or `*** FAILURE ***`. Part C is the slow one: the modular
sweep over the full range takes about seventy seconds and the rational sweep over
`p ≤ 60` about four minutes, on one core.

## Regenerating the data and the table

`depths.json` and Table 1 of Appendix A are both reproducible from scratch:

```bash
python3 compute_depths.py depths.json     # about 80 seconds
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
| `depths.json` | 214 rows `[p, k, ℓ, z]`: the type, the capture depth, and a captured vertex. |
| `requirements.txt` | `sympy`, `numpy`. |
| `CITATION.cff` | Machine-readable citation metadata. |
| `.zenodo.json` | Metadata for the Zenodo deposit. |

## Citing

Please cite the paper. If you use the code or data directly, cite this archive as well
via the DOI badge above; the DOI shown resolves to the latest version.

## License

MIT — see [LICENSE](LICENSE).
