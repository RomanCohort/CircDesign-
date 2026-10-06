#!/usr/bin/env python3
"""Summarise the accessibility panel into a single artefact for the report.

Reads accessibility.csv (construct + 58 dominators + 250 background) and the two
advantage_space_seed*.csv files, and writes accessibility_summary.json with:

  * the construct's position on each pre-registered accessibility metric
  * the domination rate under the three-objective set and under the full objective
    set (three objectives + the four declared accessibility metrics)
  * the PREFERENTIAL test, which is the honest bound on the result: does the added
    objective group remove dominators at a higher rate than it removes equally
    buildable non-dominators? If not, the drop in domination rate is the mechanical
    consequence of adding objectives rather than a finding about the construct.

Both numbers and the test are written, because the report needs all three.

Run:  python accessibility_summary.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "accessibility_summary.json"

CON_MFE = -779.4000244140625
CON_REC = 0.7956989247311828
N_ALL, N_BUILD = 4000, 1553

# (column, label, higher_is_better)
PANEL = [
    ("acc_TIR", "TIR 翻译起始区可及性", True),
    ("acc_SPACER", "间隔区可及性", True),
    ("acc_CDS", "CDS 全局可及性", True),
    ("longest_helix", "最长螺旋长度", False),
]


def main() -> None:
    acc = pd.read_csv(ROOT / "accessibility.csv", encoding="utf-8-sig")
    o = acc[acc.candidate_id == "original"].iloc[0]
    B = acc[acc.candidate_id.str.startswith("bg_")]
    Dm = acc[acc.candidate_id.str.startswith("dom_")].copy()

    frames = []
    for s in (771, 4242):
        d = pd.read_csv(ROOT / f"advantage_space_seed{s}.csv", encoding="utf-8-sig")
        d["seed"] = s
        d["pos"] = range(len(d))
        frames.append(d)
    pool = pd.concat(frames, ignore_index=True)
    dom = pool[(pool.mfe <= CON_MFE) & (pool.recall >= CON_REC)].copy()

    def all_pass(d):
        m = pd.Series(True, index=d.index)
        for c, _, hi in PANEL:
            m &= (d[c] >= o[c]) if hi else (d[c] <= o[c])
        return m

    # align the dom_* rows back to their generating rows by (seed, order within seed)
    Dm["seed"] = Dm.candidate_id.str.extract(r"dom_seed(\d+)_").astype("Int64")
    Dm["ord"] = Dm.candidate_id.str.extract(r"_(\d+)$").astype("Int64")
    aligned = []
    for s in (771, 4242):
        sub = dom[dom.seed == s].reset_index(drop=True)
        dsub = Dm[Dm.seed == s].sort_values("ord").reset_index(drop=True)
        if len(sub) != len(dsub):
            raise SystemExit(f"alignment failed for seed {s}: {len(sub)} vs {len(dsub)}")
        sub = sub.copy()
        sub["acc_pass"] = all_pass(dsub).to_numpy()
        aligned.append(sub)
    doma = pd.concat(aligned, ignore_index=True)

    dom_ok = doma[doma.correct]
    kd, nd = int(dom_ok.acc_pass.sum()), len(dom_ok)
    kb, nb = int(all_pass(B).sum()), len(B)
    _, pv = stats.fisher_exact([[kd, nd - kd], [kb, nb - kb]])
    kd_all = int(doma.acc_pass.sum())

    positions = {}
    print("本构建在各可及性指标上的位置（背景 = 250 条随机可建成设计）")
    for c, lab, hi in PANEL:
        b = B[c].to_numpy(float)
        obs = float(o[c])
        better = 100.0 * ((b > obs).mean() if hi else (b < obs).mean())
        positions[c] = {
            "label": lab, "construct": obs, "background_median": float(np.median(b)),
            "background_sd": float(b.std(ddof=1)),
            "percent_better_than_construct": better,
            "direction": "higher_is_better" if hi else "lower_is_better",
        }
        print(f"  {lab:<22} 本构建={obs:.4f}  背景中位={np.median(b):.4f}  "
              f"比本构建更好的={better:.1f}%")

    print(f"\n支配率  三目标: {len(dom_ok)}/{N_BUILD} = {100*len(dom_ok)/N_BUILD:.2f}%")
    print(f"        +可及性: {kd}/{N_BUILD} = {100*kd/N_BUILD:.2f}%")
    print(f"  可建成支配者通过率 {kd}/{nd} = {100*kd/nd:.1f}%")
    print(f"  可建成背景通过率   {kb}/{nb} = {100*kb/nb:.1f}%")
    print(f"  偏好性 Fisher p = {pv:.4f}  "
          f"{'显著' if pv < 0.05 else '不显著 —— 降幅是力学性的，不是发现'}")

    OUT.write_text(json.dumps({
        "what": "accessibility panel as a fourth objective group",
        "objective_set": [c for c, _, _ in PANEL],
        "directions": {c: ("higher" if hi else "lower") for c, _, hi in PANEL},
        "construct_accessibility": {c: float(o[c]) for c, _, _ in PANEL},
        "construct_position": positions,
        "n_buildable": N_BUILD,
        "n_all_candidates": N_ALL,
        "dominators_three_objective": len(dom_ok),
        "dominators_full_objective": kd,
        "domination_rate_three_objective_percent": 100 * len(dom_ok) / N_BUILD,
        "domination_rate_full_objective_percent": 100 * kd / N_BUILD,
        "dominators_all_three_objective": len(doma),
        "dominators_all_full_objective": kd_all,
        "domination_rate_all_three_percent": 100 * len(doma) / N_ALL,
        "domination_rate_all_full_percent": 100 * kd_all / N_ALL,
        "preferential_test": {
            "dominators_pass_percent": 100 * kd / nd,
            "background_pass_percent": 100 * kb / nb,
            "fisher_p": float(pv),
            "significant": bool(pv < 0.05),
            "interpretation": (
                "not significant -- the reduction in dominance is the mechanical "
                "consequence of enlarging the objective set, not evidence that the "
                "added objectives single out the competitors"),
        },
        "disclosures_required": [
            "the accessibility group is a fourth objective group added after the "
            "three-objective front was known; it needs out-of-sample replication",
            "the construct ranks behind the background median on three of the four "
            "accessibility metrics, so the objective does not favour it in absolute terms",
            "Pareto optimality does not require being best on any single objective; "
            "the construct is not top-ranked on any one of the seven",
        ],
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
