#!/usr/bin/env python3
"""Compute accessibility for EVERY buildable design, then re-run the multi-objective
test with the accessibility panel added.

WHY THE FULL SAMPLE
-------------------
accessibility.py measured the construct against a 250-design background and found that
the accessibility panel removes dominators preferentially (8.6% of dominators pass all
four declared metrics vs 20.0% of the background, Fisher p = 0.055). That estimate of
the resulting domination rate was an extrapolation from 250 sequences. This script
removes the extrapolation: it computes accessibility for all buildable designs and
reports the real number.

THE HONEST CAVEATS, STATED UP FRONT
-----------------------------------
1. Adding objectives mechanically reduces the number of dominators. A drop from 1.93%
   to some lower figure is therefore NOT evidence by itself; the evidence is whether
   the drop is PREFERENTIAL -- i.e. whether dominators pass the new objectives at a
   lower rate than equally-buildable non-dominators do (Fisher exact between the two).
   That comparison is the headline here, not the raw rate.
2. The accessibility panel (TIR, SPACER, CDS, longest helix) carries directions that
   were declared before any result was seen, but it is still a FOURTH objective group
   added after the three-objective front was already known. It therefore needs
   out-of-sample replication before it is treated as established.
3. The construct is WORSE than the background median on four of the five declared
   metrics, INCLUDING the one that does the discriminating. So the finding is not
   "the construct is accessible"; it is "the designs that beat it on MFE and IRES are
   even less accessible than it is". Both readings are reported.

Outputs accessibility_full.csv and accessibility_full.json.

Run:  python accessibility_full.py [--workers 28]
"""

from __future__ import annotations

import argparse
import json
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import accessibility as A

ROOT = Path(__file__).resolve().parent
OUT_CSV = ROOT / "accessibility_full.csv"
OUT_JSON = ROOT / "accessibility_full.json"

ACCESS_COLS = ["acc_TIR", "acc_SPACER", "acc_CDS", "longest_helix"]
# direction: True = higher is better, False = lower is better (longest helix)
DIR_HIGHER = {"acc_TIR": True, "acc_SPACER": True, "acc_CDS": True, "longest_helix": False}


def load_pool():
    frames = []
    for s in (771, 4242):
        p = ROOT / f"advantage_space_seed{s}.csv"
        if p.exists():
            d = pd.read_csv(p, encoding="utf-8-sig")
            d["seed"] = s
            d["pos"] = range(len(d))
            frames.append(d)
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    A._MATURE_RNA = "".join(
        l.strip() for l in A.MATURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("T", "U")
    recs = A.read_fasta(A.CDS_FASTA)
    pool_df = load_pool()

    def seq_for(seed_offset, n=2000):
        return A.makeup_shuffles(recs["original"], n, 20260804 + int(seed_offset))

    cache = {int(s): seq_for(s) for s in sorted(pool_df["seed"].unique())}

    todo = []
    for _, r in pool_df.iterrows():
        todo.append((f"s{int(r.seed)}p{int(r.pos)}", cache[int(r.seed)][int(r.pos)]))
    print(f"folding ALL {len(todo)} candidates (both seeds) on {args.workers} workers...")

    with Pool(args.workers) as pool:
        rows = pool.map(A.worker, todo, chunksize=4)
    acc = pd.DataFrame(rows)
    acc.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print(f"wrote {OUT_CSV.name}")

    # ---- verify the regenerated sequences are the right ones --------------------
    acc["seed"] = acc.candidate_id.str.extract(r"s(\d+)p").astype(int)
    acc["pos"] = acc.candidate_id.str.extract(r"p(\d+)$").astype(int)
    chk = pool_df.merge(acc[["seed", "pos", "mfe"]], on=["seed", "pos"], how="inner")
    dev = (chk.mfe_x - chk.mfe_y).abs().max()
    print(f"guard: recomputed MFE vs recorded MFE, max |deviation| = {dev:.6f} kcal/mol")
    if dev > 1e-6:
        raise SystemExit("guard failed: regenerated sequences do not match the recorded "
                         "MFE -- the sample-order replay is wrong, refusing to continue")
    print("guard passed: every regenerated sequence reproduces its recorded MFE")

    # ---- join and test ----------------------------------------------------------
    df = pool_df.merge(acc[["seed", "pos"] + ACCESS_COLS], on=["seed", "pos"], how="left")
    o = A.worker(("original", recs["original"]))
    build = df[df.correct].copy()
    dom3 = build[(build.mfe <= A.CON_MFE) & (build.recall >= A.CON_REC)]
    non3 = build.drop(dom3.index)

    def passes(d):
        m = pd.Series(True, index=d.index)
        for c in ACCESS_COLS:
            m &= (d[c] >= o[c]) if DIR_HIGHER[c] else (d[c] <= o[c])
        return m

    pk_dom = passes(dom3)
    pk_non = passes(non3)
    kd, nd = int(pk_dom.sum()), len(dom3)
    kn, nn = int(pk_non.sum()), len(non3)
    _, pv = stats.fisher_exact([[kd, nd - kd], [kn, nn - kn]])

    print("\n" + "=" * 88)
    print("可及性面板（四项，方向预先声明，全部使用）")
    print("=" * 88)
    print(f"  构造可及性: " + "  ".join(f"{c}={o[c]:.4f}" for c in ACCESS_COLS))
    print(f"\n  三目标支配者 (n={nd})  通过可及性: {kd}  ({100*kd/nd:.1f}%)")
    print(f"  三目标非支配者(n={nn}) 通过可及性: {kn}  ({100*kn/nn:.1f}%)")
    print(f"  Fisher exact p = {pv:.4f}    {'★ 偏好性显著' if pv < 0.05 else '不显著'}")

    surv = build[pk_non | pk_dom]
    n_all = len(build)
    rate_new = 100 * kd / n_all
    print(f"\n  全样本 {n_all} 条可建成设计中，同时通过三目标支配 + 可及性面板的 = {kd}")
    print(f"  → 支配率 {rate_new:.3f}%   (原三目标支配率 {100*nd/n_all:.2f}%)")
    print(f"  → 可及性把支配者从 {nd} 条降到 {kd} 条")

    print("\n  诚实读法：")
    print("    · 支配率下降本身不是证据（加目标必然如此）")
    print(f"    · 证据是【偏好性】：支配者通过率 {100*kd/nd:.1f}% vs 非支配者 {100*kn/nn:.1f}%，p={pv:.4f}")
    print("    · 且本构建自身在可及性上比中位差——发现的是'对手更差'，不是'我们更好'")

    OUT_JSON.write_text(json.dumps({
        "what": "accessibility for all buildable designs; 4th objective group added",
        "accessibility_columns": ACCESS_COLS,
        "construct_accessibility": {c: float(o[c]) for c in ACCESS_COLS},
        "n_buildable": int(n_all),
        "dominators_three_objective": nd,
        "dominators_three_objective_passing_accessibility": kd,
        "non_dominators_passing_accessibility": kn,
        "n_non_dominators": nn,
        "fisher_p_preferential": float(pv),
        "domination_rate_three_objective_percent": 100 * nd / n_all,
        "domination_rate_with_accessibility_percent": rate_new,
        "caveats": [
            "adding objectives mechanically reduces dominators; the preferential Fisher "
            "test is the evidence, not the rate drop",
            "fourth objective group added after the three-objective front was known; "
            "requires out-of-sample replication",
            "the construct is worse than the background median on 4 of 5 declared "
            "accessibility metrics, including the discriminating one",
        ],
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT_JSON.name}")


if __name__ == "__main__":
    main()
