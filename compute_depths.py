#!/usr/bin/env python3
"""
Regenerate depths.json: for every pair (p, E) with 5 <= p <= 250 and
E < Z_p^* proper with -1 in E, the capture depth d(p,E) of Definition 5.14
and the least vertex z not in {0,1} with delta_z in V_d.

    python3 compute_depths.py [depths.json]

Exact integer arithmetic throughout; no floating point.  About three minutes
on one core.
"""
import json, sys, time
from sympy import primerange
from verify import proper_subgroups_with_minus_one, capture_depth

def main(path="depths.json"):
    rows, t0 = [], time.time()
    for p in primerange(5, 251):
        for k, E in proper_subgroups_with_minus_one(p):
            d, Z = capture_depth(p, k)
            assert d is not None, (p, k)
            rows.append([p, k, d, min(Z)])
            print(f"p={p:4d} k={k:4d} d={d} z={min(Z)}  |Z|={len(Z)}"
                  f"  [{time.time()-t0:6.1f}s]", flush=True)
    with open(path, "w") as f:
        json.dump(rows, f, indent=0)
    print(f"wrote {len(rows)} rows to {path} in {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "depths.json")
