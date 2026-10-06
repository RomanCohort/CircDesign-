#!/usr/bin/env python3
"""Pool the independent advantage-space runs into one estimate.

Each run is a separate draw from the same pre-registered procedure with a different
seed offset, so pooling is legitimate and is what makes the confidence bound usable:
a single n=2000 run puts the correctness-passing domination rate at 1.30-2.52% across
runs, and only the pooled n gives a bound tight enough to interpret.

Reads advantage_space_seed*.json, writes advantage_space_pooled.json.

Run:  python pool_advantage.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "advantage_space_pooled.json"


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h) * 100, min(1.0, c + h) * 100


def main() -> None:
    runs = []
    for path in sorted(ROOT.glob("advantage_space_seed*.json")):
        d = json.loads(path.read_text(encoding="utf-8"))
        d["_source"] = path.name
        runs.append(d)

    if not runs:
        raise SystemExit("no advantage_space_seed*.json found -- run advantage_space.py first")

    print(f"pooling {len(runs)} independent run(s):")
    for d in runs:
        print(f"  {d['_source']:<34} n={d['n_generated']:<5} "
              f"passing={d['n_correctness_passing']:<5} "
              f"dominators(all)={d['dominators_all']:<4} "
              f"dominators(passing)={d['dominators_correctness_passing']}")

    n_all = sum(d["n_generated"] for d in runs)
    k_all = sum(d["dominators_all"] for d in runs)
    n_f = sum(d["n_correctness_passing"] for d in runs)
    k_f = sum(d["dominators_correctness_passing"] for d in runs)
    n_comp = sum(d["n_mfe_competitive"] for d in runs)

    lo_a, hi_a = wilson(k_all, n_all)
    lo_f, hi_f = wilson(k_f, n_f)

    print("\n" + "=" * 74)
    print("POOLED")
    print("=" * 74)
    print(f"  generated                : {n_all}")
    print(f"  correctness-passing      : {n_f} ({100*n_f/n_all:.2f}%)")
    print(f"  MFE-competitive          : {n_comp} ({100*n_comp/n_all:.2f}%)")
    print(f"  dominators (all)         : {k_all}/{n_all} = {100*k_all/n_all:.2f}%  "
          f"95%CI [{lo_a:.2f}, {hi_a:.2f}]")
    print(f"  dominators (passing)     : {k_f}/{n_f} = {100*k_f/n_f:.2f}%  "
          f"95%CI [{lo_f:.2f}, {hi_f:.2f}]")

    rule = json.loads((ROOT / "advantage_space_seed4242.json").read_text(
        encoding="utf-8"))["pre_registration"]["decision_rule"]
    supported = (k_f <= 2) and (n_f >= 500) and (hi_f < 2.0)
    print(f"\n  rule: {rule}")
    print(f"  -> {'SUPPORTED' if supported else 'NOT SUPPORTED'}")
    if not supported and 100 * k_f / n_f < 2.0:
        print("  (point estimate meets the 2% threshold; the 95% upper bound does not. "
              "That is a power limit at this n, not an absence of effect.)")

    OUT.write_text(json.dumps({
        "what": "pooled estimate across independent advantage-space runs",
        "runs": [{k: v for k, v in d.items() if k != "results"} for d in runs],
        "n_generated": n_all,
        "n_correctness_passing": n_f,
        "n_mfe_competitive": n_comp,
        "dominators_all": k_all,
        "dominators_correctness_passing": k_f,
        "rate_all_percent": 100 * k_all / n_all,
        "rate_passing_percent": 100 * k_f / n_f,
        "ci_all_percent": [lo_a, hi_a],
        "ci_passing_percent": [lo_f, hi_f],
        "decision_rule": rule,
        "supported": bool(supported),
        "point_estimate_meets_threshold": bool(100 * k_f / n_f < 2.0),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
