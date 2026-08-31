#!/usr/bin/env python3
"""
Regenerate the LaTeX source of Table 1 (Appendix A) from depths.json.

Usage:  python3 make_table.py > appendix_table.tex
"""
import json
import sys

NCOL = 4


def main(path="depths.json"):
    rows = json.load(open(path))
    per = -(-len(rows) // NCOL)
    cols = [rows[i * per:(i + 1) * per] for i in range(NCOL)]
    out = [r"\begin{tabular}{" + "rrrr" * NCOL + "}", r"\toprule",
           " & ".join([r"$p$ & $k$ & $\ell$ & $z$"] * NCOL) + r"\\", r"\midrule"]
    for i in range(per):
        cells = []
        for c in cols:
            if i < len(c):
                p, k, d, z = c[i]
                cells.append(f"${p}$ & ${k}$ & ${d}$ & ${z}$")
            else:
                cells.append(" & & & ")
        out.append(" & ".join(cells) + r"\\")
    out += [r"\bottomrule", r"\end{tabular}"]
    print("\n".join(out))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "depths.json")
