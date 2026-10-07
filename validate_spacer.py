#!/usr/bin/env python3
"""Out-of-sample validation of the IRES-CDS spacer accessibility result.

THE CLAIM UNDER TEST
--------------------
Report section 5.3 reports the construct's IRES-CDS spacer accessibility at the 91.2nd
percentile of buildable synonymous designs, and it is the one positive structural
result in the report. It was measured against a single 250-design background drawn
from the pooled seed-771+4242 samples. Everything else in this project was replicated
out of sample (the domination rate across three seeds, the fold, the sampler); this
was not.

WHAT MAKES THIS A REAL TEST
---------------------------
The background is drawn from an UNUSED seed (offset 9999), so it is a fresh draw from
the population, not a re-draw from the same 4,000 sequences. Two things are reported:

  1. the construct's percentile on each of the four declared accessibility metrics in
     the new background, against the 250-design background used in the report;
  2. a bootstrap interval on each percentile, so "is 91.2 reproducible" is answered
     with a number rather than by eyeballing two point estimates.

PRE-REGISTERED EXPECTATION
--------------------------
spacer accessibility stays in the top decile (percentile >= 90).

The splice happens in the parent; workers receive finished sequences, because
multiprocessing spawn on Windows gives them an empty template (this failure mode
reached a published number once).

Outputs validate_spacer.json.

Run:  python validate_spacer.py [--n-background 500] [--workers 28]
"""

from __future__ import annotations

import argparse
import json
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np

import accessibility as A
from advantage_space import makeup_shuffles, qc_fast

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "validate_spacer.json"

NEW_SEED_OFFSET = 9999          # unused by any other analysis
METRICS = [("acc_SPACER", "IRES–CDS 间隔区可及性", True),
           ("acc_TIR", "TIR 翻译起始区可及性", True),
           ("acc_CDS", "CDS 全局可及性", True),
           ("longest_helix", "最长螺旋长度", False)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-background", type=int, default=500)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    A._MATURE_RNA = "".join(
        l.strip() for l in A.MATURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("T", "U")
    recs = A.read_fasta(A.CDS_FASTA)

    # ---- the fresh background: a new population draw, QC only (no folding) --------
    print(f"drawing a fresh background from UNUSED seed offset {NEW_SEED_OFFSET}...")
    pool = makeup_shuffles(recs["original"], 2000, 20260804 + NEW_SEED_OFFSET)
    qc = [qc_fast(s) for s in pool]
    buildable = [s for s, m in zip(pool, qc)
                 if m["splice_donors"] == 0 and m["polya_signals"] == 0]
    rng = np.random.default_rng(20260804)
    if len(buildable) > args.n_background:
        idx = rng.choice(len(buildable), size=args.n_background, replace=False)
        buildable = [buildable[i] for i in idx]
    print(f"  {len(buildable)} buildable designs in the fresh background")

    # ---- fold: parent splices, workers only fold ---------------------------------
    todo = [("original", A.build_mature(recs["original"]))]
    todo += [(f"bg_{i:04d}", A.build_mature(s)) for i, s in enumerate(buildable)]
    if {len(s) for _, s in todo} != {len(A._MATURE_RNA)}:
        raise SystemExit("guard failed: spliced lengths wrong")
    print(f"  folding {len(todo)} sequences on {args.workers} workers...")
    with Pool(args.workers) as p:
        rows = p.map(A.worker, todo, chunksize=4)

    o = rows[0]
    if abs(o["mfe"] - -779.4000244140625) > 1e-6:
        raise SystemExit(f"guard failed: construct folds to {o['mfe']}, expected -779.4000")
    print(f"  guard passed: construct MFE {o['mfe']:.4f}")
    B = rows[1:]

    # the report's background, for comparison
    ref = ROOT / "accessibility.csv"
    old = {}
    if ref.exists():
        import pandas as pd
        a = pd.read_csv(ref, encoding="utf-8-sig")
        a = a[a.candidate_id.str.startswith("bg_")]
        for col, _, _ in METRICS:
            if col in a.columns:
                old[col] = a[col].to_numpy(float)

    rng2 = np.random.default_rng(7)
    print("\n" + "=" * 96)
    print("SPACER 结果的样本外验证（背景 = 全新种子的可建成设计）")
    print("=" * 96)
    print(f"  {'指标':<26}{'构造':>10}{'新背景中位':>12}{'新样本位次':>12}"
          f"{'95% 区间':>18}{'报告值':>10}")
    res = {}
    for col, lab, hi in METRICS:
        v = np.array([r[col] for r in B], dtype=float)
        obs = float(o[col])
        # lower percentile = construct better (for hi=True, "better" means fewer above)
        pctl = 100.0 * ((v < obs).mean() if hi else (v > obs).mean())
        boots = []
        n = len(v)
        for _ in range(4000):
            s = v[rng2.integers(0, n, n)]
            boots.append(100.0 * ((s < obs).mean() if hi else (s > obs).mean()))
        lo, hii = np.percentile(boots, [2.5, 97.5])
        oldp = None
        if col in old:
            ov = old[col]
            oldp = 100.0 * ((ov < obs).mean() if hi else (ov > obs).mean())
        res[col] = {"construct": obs, "fresh_median": float(np.median(v)),
                    "fresh_percentile": float(pctl),
                    "bootstrap_ci": [float(lo), float(hii)],
                    "report_percentile": float(oldp) if oldp is not None else None,
                    "n_background": int(n)}
        old_s = f"{oldp:.1f}%" if oldp is not None else "—"
        print(f"  {lab:<26}{obs:>10.4f}{np.median(v):>12.4f}{pctl:>11.1f}%"
              f"{f'[{lo:.1f}, {hii:.1f}]':>18}{old_s:>10}")

    sp = res["acc_SPACER"]
    ok = sp["fresh_percentile"] >= 90.0
    print(f"\n  预注册预期：SPACER 保持在最高十分位（>= 90）")
    print(f"  → {'REPRODUCED' if ok else 'NOT REPRODUCED'}  "
          f"(新样本 {sp['fresh_percentile']:.1f}%，95% 区间 "
          f"[{sp['bootstrap_ci'][0]:.1f}, {sp['bootstrap_ci'][1]:.1f}]，"
          f"报告值 {sp['report_percentile']:.1f}%)")

    OUT.write_text(json.dumps({
        "what": "out-of-sample validation of the IRES-CDS spacer accessibility result",
        "fresh_seed_offset": NEW_SEED_OFFSET,
        "n_background": len(B),
        "metrics": res,
        "prediction": "spacer accessibility stays in the top decile (>= 90)",
        "reproduced": bool(ok),
        "note": ("the report's spacer figure came from a single 250-design background "
                 "with no out-of-sample check; the rest of the project's results were "
                 "replicated across seeds"),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
