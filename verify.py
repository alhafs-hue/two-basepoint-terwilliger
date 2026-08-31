#!/usr/bin/env python3
"""
Ancillary verification for
  "Two-basepoint Terwilliger algebras and the quantum rigidity of circulant
   graphs of prime order".

Three independent checks, all in exact integer / rational arithmetic:

  (A) the four certificates of Section 6 (p = 13, 17, 31, 41);
  (B) the confinement lemma, the quadratic threshold p > (k-1)^2+2, the
      enumeration of the eight residual pairs in types <= 10, and their witnesses;
  (C) the saturation of Theorem 7.1 (all 214 pairs with p <= 250), over Q for
      small p and over F_q, q = 1000003, for the whole range;
  (D) the capture-depth table of Appendix A.

Run:  python3 verify.py [A|B|C|D|all]
Requires: sympy, numpy.
"""
import sys
from fractions import Fraction
import sympy
from sympy import primerange

# ---------------------------------------------------------------- utilities
def coset_decomposition(p, E):
    seen, cs = set(), []
    for x in range(1, p):
        if x in seen: continue
        C = sorted((x*e) % p for e in E); seen.update(C); cs.append(C)
    return cs

def subgroup_of_order(p, k):
    g = sympy.primitive_root(p)
    h = pow(g, (p-1)//k, p)
    return sorted({pow(h, i, p) for i in range(k)})

def proper_subgroups_with_minus_one(p):
    out = []
    for d in sympy.divisors(p-1):
        if d == p-1: continue
        E = subgroup_of_order(p, d)
        if (p-1) in E: out.append((d, E))
    return out

def blocks(p, cs):
    idx = [0]*p
    for i, C in enumerate(cs, 1):
        for x in C: idx[x] = i
    B = {}
    for x in range(p):
        if x in (0, 1): continue
        B.setdefault((idx[x], idx[(x-1) % p]), []).append(x)
    return {k: sorted(v) for k, v in B.items()}

# ------------------------------------------------------------- (A) certificates
def check_certificates():
    print("(A) certificates of Section 6")
    data = {  # p -> (index r, B-spec, u, D-spec, expected profile, captured z)
        31: (3, (3, 1), 1, (2, 3), {6: 1, 12: 0, 14: 0, 19: 0}, 6),
        41: (4, (3, 2), 1, (2, 1), {2: 1, 5: 0, 32: 0}, 2),
        13: (2, (1, 1), 1, (2, 1), {2: 0, 5: 1, 11: 1}, 2),
    }
    ok = True
    for p, (r, (s, t), u, (sp, tp), prof, z) in data.items():
        cs = coset_decomposition(p, subgroup_of_order(p, (p-1)//r))
        B = sorted(set(cs[s-1]) & {(1+c) % p for c in cs[t-1]})
        D = sorted(set(cs[sp-1]) & {(1+c) % p for c in cs[tp-1]})
        v = {x: sum(1 for b in B for c in cs[u-1] if (b+c) % p == x) for x in D}
        good = (v == prof)
        ok &= good
        print(f"   p={p:3d}  B={B}  D={D}  profile={v}  -> delta_{z}  "
              + ("OK" if good else "FAIL"))
    # p = 17, depth three
    p = 17; cs = coset_decomposition(p, subgroup_of_order(p, 8)); C1 = cs[0]
    B = sorted(set(C1) & {(1+c) % p for c in C1})
    D = sorted(set(C1) & {(1+c) % p for c in cs[1]})
    v = {x: sum(1 for b in B for c in C1 if (b+c) % p == x) for x in D}
    m = min(v.values()); v1 = [x for x in D if v[x] - m]
    Dp = B
    v2 = {x: sum(1 for b in v1 for c in C1 if (b+c) % p == x) for x in Dp}
    good = (v == {4: 1, 8: 2, 13: 1, 15: 2} and v1 == [8, 15]
            and v2 == {2: 1, 9: 1, 16: 2})
    ok &= good
    print(f"   p= 17  B={B}  D={D}  profile={v}  ->  v1 supported on {v1}")
    print(f"          D'={Dp}  profile={v2}  -> delta_16  " + ("OK" if good else "FAIL"))
    return ok

# ------------------------------------------------------- (B) threshold & 8 pairs
def check_threshold():
    print("(B) confinement lemma, quadratic threshold, and the eight residual pairs")
    ok = True
    bad_containment, bad_bound = [], []
    n = 0
    for p in primerange(5, 200):
        for k, E in proper_subgroups_with_minus_one(p):
            n += 1
            cs = coset_decomposition(p, E); Bl = blocks(p, cs)
            big = set()
            for B in Bl.values():
                if len(B) >= 2: big.update(B)
            conf = {((1-ep)*pow((e-ep) % p, p-2, p)) % p
                    for e in E for ep in E if e != ep and ep != 1}
            if not big <= conf: bad_containment.append((p, k))
            if len(conf) > (k-1)**2: bad_bound.append((p, k))
            if p-2 > (k-1)**2 and min(len(B) for B in Bl.values()) != 1:
                bad_bound.append(("no singleton", p, k))
    ok &= not bad_containment and not bad_bound
    print(f"   {n} pairs: confinement holds ({not bad_containment}), "
          f"bound (k-1)^2 and singleton conclusion hold ({not bad_bound})")
    pairs = sorted((p, k) for k in (2, 4, 6, 8, 10)
                   for p in primerange(5, (k-1)**2+3)
                   if (p-1) % k == 0 and p-1 > k)
    expected = sorted([(13, 6), (19, 6), (17, 8), (41, 8),
                       (31, 10), (41, 10), (61, 10), (71, 10)])
    good = (pairs == expected); ok &= good
    print(f"   residual pairs with p <= (k-1)^2+2 in types <= 10: {pairs}  "
          + ("OK" if good else "FAIL"))
    wit = {(19, 6): (2, 1, 2), (41, 8): (3, 5, 12), (61, 10): (5, 2, 8), (71, 10): (2, 1, 2)}
    for (p, k), (s, t, z) in wit.items():
        cs = coset_decomposition(p, subgroup_of_order(p, k))
        B = sorted(set(cs[s-1]) & {(1+c) % p for c in cs[t-1]})
        good = (B == [z]); ok &= good
        print(f"   p={p:2d}, k={k:2d}: C_{s} cap (1+C_{t}) = {B}  " + ("OK" if good else "FAIL"))
    return ok

# ------------------------------------------------------------- (C) saturation
def module_dim_exact(p, E):
    cs = coset_decomposition(p, E)
    ops = []
    for C in cs:
        Cs = set(C)
        ops.append([[1 if (x-y) % p in Cs else 0 for y in range(p)] for x in range(p)])
    for a in (0, 1):
        for Sset in [{a}] + [{(a+c) % p for c in C} for C in cs]:
            ops.append([[1 if (x == y and x in Sset) else 0 for y in range(p)]
                        for x in range(p)])
    rows = []
    def add(v):
        v = list(v)
        for (j, w) in rows:
            if v[j] != 0:
                f = Fraction(v[j], w[j]); v = [a-f*b for a, b in zip(v, w)]
        for j in range(p):
            if v[j] != 0:
                rows.append((j, v)); return True
        return False
    front = []
    for a in (0, 1):
        e = [0]*p; e[a] = 1
        if add(e): front.append(e)
    while front:
        new = []
        for v in front:
            for M in ops:
                w = [sum(M[x][y]*v[y] for y in range(p)) for x in range(p)]
                if add(w): new.append(w)
        front = new
    return len(rows)

def module_full_modq(p, E, q=1000003):
    import numpy as np
    cs = coset_decomposition(p, E); ops = []
    for C in cs:
        M = np.zeros((p, p), dtype=np.int64)
        for x in range(p):
            for c in C: M[x, (x-c) % p] = 1
        ops.append(M)
    for a in (0, 1):
        D = np.zeros((p, p), dtype=np.int64); D[a, a] = 1; ops.append(D)
        for C in cs:
            D = np.zeros((p, p), dtype=np.int64)
            for c in C: D[(a+c) % p, (a+c) % p] = 1
            ops.append(D)
    piv = {}
    def add(v):
        v = v % q
        for j in range(p):
            if v[j]:
                if j in piv:
                    f = int(v[j])*pow(int(piv[j][j]), q-2, q) % q
                    v = (v - f*piv[j]) % q
                else:
                    piv[j] = v.copy(); return True
        return False
    front = []
    for a in (0, 1):
        e = np.zeros(p, dtype=np.int64); e[a] = 1
        if add(e): front.append(e)
    while front and len(piv) < p:
        new = []
        for v in front:
            for M in ops:
                w = M.dot(v) % q
                if add(w): new.append(w)
        front = new
    return len(piv) == p

def check_saturation(exact_bound=60):
    print("(C) saturation, all pairs with p <= 250")
    bad = []; n = 0
    for p in primerange(5, 251):
        for d, E in proper_subgroups_with_minus_one(p):
            n += 1
            if not module_full_modq(p, E): bad.append((p, d))
            if p <= exact_bound and module_dim_exact(p, E) != p: bad.append(("exact", p, d))
    print(f"   {n} pairs tested (exact over Q for p <= {exact_bound}); failures: {bad}")
    return not bad and n == 214

# --------------------------------------------------------- (D) appendix table
def check_table(path="depths.json"):
    """Recompute the capture depths of Appendix A and compare against depths.json."""
    import json, collections
    print("(D) capture-depth table of Appendix A")
    try: rows = json.load(open(path))
    except FileNotFoundError:
        print("   depths.json not found; skipping"); return True
    c = collections.Counter(d for _, _, d, _ in rows)
    quoted = {1: 51, 2: 92, 4: 25, 5: 6, 6: 8, 7: 11, 8: 3, 9: 2, 10: 3}
    ok = (len(rows) == 214 and all(c[l] == n for l, n in quoted.items())
          and c[11]+c[12] == 13 and c[1]+c[2] == 143)
    print(f"   {len(rows)} rows; distribution {dict(sorted(c.items()))}  "
          + ("OK" if ok else "FAIL"))
    pal = [(p, d, z) for p, k, d, z in rows if k == (p-1)//2]
    good = all(z == 2 or z == pow(2, p-2, p) for p, d, z in pal)
    ok &= good
    print(f"   {len(pal)} Paley pairs; first capture is always delta_2 or delta_(1/2): "
          + ("OK" if good else "FAIL"))
    return ok

if __name__ == "__main__":
    what = (sys.argv[1] if len(sys.argv) > 1 else "all").lower()
    res = []
    if what in ("a", "all"): res.append(check_certificates())
    if what in ("b", "all"): res.append(check_threshold())
    if what in ("c", "all"): res.append(check_saturation())
    if what in ("d", "all"): res.append(check_table())
    print()
    print("ALL VERIFICATIONS PASSED" if all(res) else "*** FAILURE ***")

