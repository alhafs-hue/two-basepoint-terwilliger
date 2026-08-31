#!/usr/bin/env python3
"""
Regenerate depths.json: the capture depth and a captured vertex for every pair
(p, E) with 5 <= p <= 250, E a proper subgroup of Z_p^* containing -1.

The capture depth l(p,E) is the least l such that delta_z, for some z not in
{0,1}, lies in the span of the vectors g_{i_l} ... g_{i_1} delta with
delta in {delta_0, delta_1} and products of length at most l, where the g_i
are the generators T_s, E*_s(0), E*_s(1) of the two-basepoint Terwilliger
algebra.

Linear algebra is done over F_q with q = 1000003.  This is sound for the
purpose at hand: membership of delta_z in the span over F_q is implied by
membership over Q, so the depth reported here is an upper bound for the true
depth, and every depth reported is independently confirmed over Q by
verify.py part C, which shows the module is full in every case.

Usage:  python3 compute_depths.py [output.json]
Runtime: about 80 seconds on one core.
"""
import json
import sys
import time

import numpy as np
import sympy
from sympy import primerange

Q = 1000003


def coset_decomposition(p, E):
    seen, cs = set(), []
    for x in range(1, p):
        if x in seen:
            continue
        C = sorted((x * e) % p for e in E)
        seen.update(C)
        cs.append(C)
    return cs


def proper_subgroups_with_minus_one(p):
    g = sympy.primitive_root(p)
    out = []
    for d in sympy.divisors(p - 1):
        if d == p - 1:
            continue
        h = pow(g, (p - 1) // d, p)
        E = sorted({pow(h, i, p) for i in range(d)})
        if (p - 1) in E:
            out.append((d, E))
    return out


def generators(p, cs):
    """T_s for each class, then E*_s(0) and E*_s(1), as 0/1 matrices."""
    ops = []
    for C in cs:
        M = np.zeros((p, p), dtype=np.int64)
        for x in range(p):
            for c in C:
                M[x, (x - c) % p] = 1
        ops.append(M)
    for a in (0, 1):
        D = np.zeros((p, p), dtype=np.int64)
        D[a, a] = 1
        ops.append(D)                      # E*_0(a), the class C_0 = {0}
        for C in cs:
            D = np.zeros((p, p), dtype=np.int64)
            for c in C:
                D[(a + c) % p, (a + c) % p] = 1
            ops.append(D)
    return ops


def capture_depth(p, E, maxdepth=40):
    cs = coset_decomposition(p, E)
    ops = generators(p, cs)
    piv = {}

    def reduce(v):
        v = v % Q
        for j in range(p):
            if v[j]:
                if j in piv:
                    f = int(v[j]) * pow(int(piv[j][j]), Q - 2, Q) % Q
                    v = (v - f * piv[j]) % Q
                else:
                    return j, v
        return None, v

    def insert(v):
        j, w = reduce(v)
        if j is None:
            return False
        piv[j] = w
        return True

    def in_span(v):
        return reduce(v)[0] is None

    frontier = []
    for a in (0, 1):
        e = np.zeros(p, dtype=np.int64)
        e[a] = 1
        if insert(e):
            frontier.append(e)

    for depth in range(1, maxdepth + 1):
        new = []
        for v in frontier:
            for M in ops:
                w = M.dot(v) % Q
                if insert(w):
                    new.append(w)
        frontier = new
        for z in range(p):
            if z in (0, 1):
                continue
            e = np.zeros(p, dtype=np.int64)
            e[z] = 1
            if in_span(e):
                return depth, z
        if not frontier:
            return None, None
    return None, None


def main(out="depths.json"):
    rows, t0 = [], time.time()
    for p in primerange(5, 251):
        for k, E in proper_subgroups_with_minus_one(p):
            d, z = capture_depth(p, E)
            rows.append((p, k, d, z))
            print(f"p={p:3d}  k={k:3d}  depth={d}  z={z}", flush=True)
    json.dump(rows, open(out, "w"))
    print(f"\n{len(rows)} pairs written to {out} in {time.time()-t0:.0f}s")
    if any(d is None for _, _, d, _ in rows):
        print("WARNING: some pair reached maxdepth without a capture")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "depths.json")
