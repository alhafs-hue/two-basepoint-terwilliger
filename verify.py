#!/usr/bin/env python3
"""
Ancillary verification for
  "Quantum rigidity of prime-order circulants via two-basepoint Terwilliger
   algebras".

Four independent checks, all in exact integer / rational arithmetic:

  (A) the four certificates of Section 7 (p = 13, 17, 31, 41), and their
      optimality;
  (B) the confinement lemma, the quadratic threshold p > (k-1)(k-2)+2, the
      sharpness of the count (k-1)(k-2), the enumeration of the eight residual
      pairs in types <= 10, and their four singleton-block witnesses;
  (C) the saturation of Theorem 8.1 (all 214 pairs with p <= 250), over Q for
      small p and over F_q, q = 1000003, for the whole range;
  (D) the capture depths of Definition 5.14 and Appendix A, recomputed from
      scratch in exact integer arithmetic, together with the structural and
      harmonic observations of Section 8.2.

Run:  python3 verify.py [A|B|C|D|all]
Requires: sympy, numpy.
"""
import sys
from fractions import Fraction
from math import gcd
import sympy
from sympy import primerange

# ---------------------------------------------------------------- utilities
def coset_decomposition(p, E):
    """Cosets of E in Z_p^*, ordered by least element; the first is E itself."""
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
    """The blocks C_s cap (1+C_t) partitioning Z_p \\ {0,1}, keyed by (s,t)."""
    idx = [0]*p
    for i, C in enumerate(cs, 1):
        for x in C: idx[x] = i
    B = {}
    for x in range(p):
        if x in (0, 1): continue
        B.setdefault((idx[x], idx[(x-1) % p]), []).append(x)
    return {k: sorted(v) for k, v in B.items()}

def all_blocks(p, cs):
    """All blocks of the basepoint partition P, including {0} and {1}."""
    return [[0], [1]] + [B for B in blocks(p, cs).values()]

# ------------------------------------------------------------- (A) certificates
def check_certificates():
    print("(A) certificates of Section 7")
    data = {  # p -> (r, B-spec, u, D-spec, expected profile, captured z)
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
    # optimality: none of the four pairs has a singleton block, so d >= 2
    good = True
    for p, k in ((13, 6), (17, 8), (31, 10), (41, 10)):
        cs = coset_decomposition(p, subgroup_of_order(p, k))
        good &= min(len(B) for B in blocks(p, cs).values()) >= 2
    ok &= good
    print("   no singleton block for (13,6), (17,8), (31,10), (41,10), so d >= 2  "
          + ("OK" if good else "FAIL"))
    return bool(ok)

# ------------------------------------------------------- (B) threshold & 8 pairs
def check_threshold():
    print("(B) confinement lemma, quadratic threshold, and the eight residual pairs")
    ok = True
    bad_containment, bad_bound, tight = [], [], 0
    n = 0
    for p in primerange(5, 200):
        for k, E in proper_subgroups_with_minus_one(p):
            n += 1
            cs = coset_decomposition(p, E); Bl = blocks(p, cs)
            big = set()
            for B in Bl.values():
                if len(B) >= 2: big.update(B)
            # Lemma 6.1: e, e' both != 1 and e != e'
            conf = {((1-ep)*pow((e-ep) % p, p-2, p)) % p
                    for e in E if e != 1 for ep in E if ep != 1 and e != ep}
            if not big <= conf: bad_containment.append((p, k))
            if len(big) > (k-1)*(k-2): bad_bound.append((p, k))
            if len(big) == (k-1)*(k-2) and k > 2: tight += 1
            if p-2 > (k-1)*(k-2) and min(len(B) for B in Bl.values()) != 1:
                bad_bound.append(("no singleton", p, k))
    ok &= not bad_containment and not bad_bound
    print(f"   {n} pairs with p < 200: confinement holds ({not bad_containment}), "
          f"bound (k-1)(k-2) and singleton conclusion hold ({not bad_bound})")
    print(f"   the bound (k-1)(k-2) is attained for {tight} of them (Remark 6.3)")
    pairs = sorted((p, k) for k in (2, 4, 6, 8, 10)
                   for p in primerange(5, (k-1)*(k-2)+3)
                   if (p-1) % k == 0 and p-1 > k)
    expected = sorted([(13, 6), (19, 6), (17, 8), (41, 8),
                       (31, 10), (41, 10), (61, 10), (71, 10)])
    good = (pairs == expected); ok &= good
    print(f"   residual pairs with p <= (k-1)(k-2)+2 in types <= 10: {pairs}  "
          + ("OK" if good else "FAIL"))
    # the four singleton-block witnesses, in the multiplier form of Proposition 7.3
    wit = {(19, 6): (2, 1, 2), (41, 8): (4, 8, 12), (61, 10): (8, 2, 8), (71, 10): (2, 1, 2)}
    for (p, k), (a, b, z) in wit.items():
        E = subgroup_of_order(p, k)
        S = sorted({(a*e) % p for e in E} & {(1 + b*e) % p for e in E})
        good = (S == [z]); ok &= good
        aa = "E" if a == 1 else f"{a}E"; bb = "E" if b == 1 else f"{b}E"
        print(f"   p={p:2d}, k={k:2d}: ({aa}) cap (1+{bb}) = {S}  "
              + ("OK" if good else "FAIL"))
    return bool(ok)

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

# --------------------------------------------------------- (D) capture depths
class _Basis:
    """Echelon basis over Z with primitive rows; exact, no floating point."""
    __slots__ = ("p", "piv")
    def __init__(self, p): self.p = p; self.piv = {}
    @staticmethod
    def _prim(v):
        g = 0
        for x in v:
            if x: g = gcd(g, abs(x))
        if g > 1: v = [x//g for x in v]
        for x in v:
            if x:
                if x < 0: v = [-y for y in v]
                break
        return v
    def _reduce(self, v):
        for i in range(self.p):
            if v[i]:
                row = self.piv.get(i)
                if row is None: return self._prim(v), i
                g = gcd(row[i], v[i]); m1, m2 = row[i]//g, v[i]//g
                v = self._prim([m1*v[j] - m2*row[j] for j in range(self.p)])
        return v, -1
    def insert(self, v):
        v, i = self._reduce(list(v))
        if i < 0: return False
        self.piv[i] = v; return True
    def contains(self, v):
        return self._reduce(list(v))[1] < 0

def capture_depth(p, k, dmax=40):
    """Least d with a new point mass in V_d = Delta . B . V_{d-1}, and the
    full set Z of vertices captured at that depth.  Exact integer arithmetic."""
    E = subgroup_of_order(p, k); cs = coset_decomposition(p, E)
    Bl = all_blocks(p, cs)
    # V_1 = span of the block indicators (Lemma 5.13)
    sing = sorted(B[0] for B in Bl if len(B) == 1 and B[0] not in (0, 1))
    if sing: return 1, sing
    cur = []
    e0 = [0]*p; e0[0] = 1; e1 = [0]*p; e1[1] = 1
    base = _Basis(p)
    for v in (e0, e1):
        if base.insert(v): cur.append(v)
    for d in range(1, dmax+1):
        S = _Basis(p); Sv = []
        for v in cur:
            cand = [list(v)]
            for C in cs:
                w = [0]*p
                for c in C:
                    for x in range(p): w[x] += v[(x-c) % p]
                cand.append(w)
            for w in cand:
                if any(w) and S.insert(w): Sv.append(w)
        V = _Basis(p); Vv = []
        for w in Sv:
            for B in Bl:
                u = [0]*p
                for x in B: u[x] = w[x]
                if any(u) and V.insert(u): Vv.append(u)
        Z = []
        for z in range(p):
            if z in (0, 1): continue
            e = [0]*p; e[z] = 1
            if V.contains(e): Z.append(z)
        if Z: return d, Z
        cur = Vv
    return None, []

def check_table(path="depths.json"):
    """Recompute every capture depth from scratch and compare with depths.json."""
    import json, collections
    print("(D) capture depths of Definition 5.14 and Appendix A")
    rows = []
    for p in primerange(5, 251):
        for k, E in proper_subgroups_with_minus_one(p):
            d, Z = capture_depth(p, k)
            rows.append((p, k, d, min(Z), Z))
    ok = all(d is not None for _, _, d, _, _ in rows) and len(rows) == 214
    c = collections.Counter(d for _, _, d, _, _ in rows)
    quoted = {1: 143, 2: 31, 3: 19, 4: 5, 5: 11, 6: 5}
    good = (dict(c) == quoted); ok &= good
    print(f"   {len(rows)} pairs recomputed; distribution {dict(sorted(c.items()))}  "
          + ("OK" if good else "FAIL"))
    # structural implications of Section 8.2
    r_of = lambda p, k: (p-1)//k
    good = (all(r_of(p, k) <= 7 for p, k, d, _, _ in rows if d >= 2)
            and all(r_of(p, k) <= 4 for p, k, d, _, _ in rows if d >= 3)
            and all(r_of(p, k) == 2 for p, k, d, _, _ in rows if d >= 4))
    ok &= good
    print("   d>=2 => r<=7,  d>=3 => r<=4,  d>=4 => r=2  " + ("OK" if good else "FAIL"))
    n48 = sum(1 for p, k, d, _, _ in rows if r_of(p, k) >= 3 and d >= 2)
    print(f"   pairs with r>=3 and d>=2: {n48}  " + ("OK" if n48 == 48 else "FAIL"))
    ok &= (n48 == 48)
    # Paley family
    pal = [(p, d, Z) for p, k, d, _, Z in rows if k == (p-1)//2]
    seq = [d for _, d, _ in pal]
    good = (len(pal) == 24 and seq == sorted(seq) and max(seq) == 6)
    ok &= good
    print(f"   {len(pal)} Paley pairs, depths {seq} non-decreasing  "
          + ("OK" if good else "FAIL"))
    inv2 = lambda p: pow(2, p-2, p)
    good = all(inv2(p) in Z for p, d, Z in pal); ok &= good
    print("   2^{-1} captured at minimal depth in every Paley pair  "
          + ("OK" if good else "FAIL"))
    trip = [p for p, d, Z in pal if sorted(Z) == sorted({2, p-1, inv2(p)})]
    good = (trip == [5, 13, 17, 53, 61, 157, 173, 181]); ok &= good
    print(f"   Z = harmonic triple exactly for p in {trip}  " + ("OK" if good else "FAIL"))
    nall = sum(1 for p, k, d, _, Z in rows if inv2(p) in Z)
    print(f"   2^{{-1}} captured at minimal depth in {nall} of 214 pairs  "
          + ("OK" if nall == 197 else "FAIL"))
    ok &= (nall == 197)
    # beta-stability (Lemma 5.15)
    good = all(set(Z) == {(1-z) % p for z in Z} for p, k, d, _, Z in rows); ok &= good
    print("   every capture set is stable under z -> 1-z  " + ("OK" if good else "FAIL"))
    # agreement with the stored table
    try:
        stored = [tuple(t) for t in json.load(open(path))]
    except FileNotFoundError:
        print("   depths.json not found; skipping comparison"); return bool(ok)
    mine = [(p, k, d, z) for p, k, d, z, _ in rows]
    good = (sorted(stored) == sorted(mine)); ok &= good
    print(f"   agrees with {path} row by row  " + ("OK" if good else "FAIL"))
    return bool(ok)

if __name__ == "__main__":
    what = (sys.argv[1] if len(sys.argv) > 1 else "all").lower()
    res = []
    if what in ("a", "all"): res.append(check_certificates())
    if what in ("b", "all"): res.append(check_threshold())
    if what in ("c", "all"): res.append(check_saturation())
    if what in ("d", "all"): res.append(check_table())
    print()
    print("ALL VERIFICATIONS PASSED" if all(res) else "*** FAILURE ***")
