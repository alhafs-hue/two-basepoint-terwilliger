#!/usr/bin/env python3
"""
Regenerate the LaTeX source of Table 1 (Appendix A) from depths.json.

    python3 make_table.py > appendix_table.tex

Emits the table environment exactly as it appears in the paper.
"""
import json, sys

NCOL = 4

def main(path="depths.json"):
    rows = sorted(tuple(t) for t in json.load(open(path)))
    n = len(rows); nrow = -(-n // NCOL)
    cols = [rows[i*nrow:(i+1)*nrow] for i in range(NCOL)]
    out = []
    out.append(r"\begin{table}[htbp]")
    out.append(r"\centering\footnotesize")
    out.append(r"\begin{tabular}{" + "rrrr"*NCOL + "}")
    out.append(r"\toprule")
    out.append(" & ".join([r"$p$ & $k$ & $d$ & $z$"]*NCOL) + r"\\")
    out.append(r"\midrule")
    for i in range(nrow):
        cells = []
        for c in cols:
            if i < len(c):
                p, k, d, z = c[i]
                cells += [f"${p}$", f"${k}$", f"${d}$", f"${z}$"]
            else:
                cells += [" "]*4
        out.append(" & ".join(cells) + r"\\")
    out.append(r"\bottomrule")
    out.append(r"\end{tabular}")
    out.append(r"\caption{Convolution depth $d=d(p,E)$ and the least vertex $z$ with "
               r"$\delta_{z}\in V_{d}$,")
    out.append(r"for all $214$ pairs $(p,E)$ with $p\le250$; here $k=|E|$. Rows with $d=1$ "
               r"are exactly those")
    out.append(r"settled by the depth-one criterion of Corollary~\ref{cor:chassaniol}.}")
    out.append(r"\label{tab:depths}")
    out.append(r"\end{table}")
    print("\n".join(out))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "depths.json")
