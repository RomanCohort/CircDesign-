#!/usr/bin/env python3
"""Evaluate the new IRES metrics under the correct null models.

Two questions this answers, in order:

1. DOES THE NEW METRIC SEE ANYTHING THE OLD ONE DIDN'T?
   Spearman correlation between each new metric and the L_IRES it replaces.  A
   replacement that correlates ~1.0 with the thing it replaces has changed the
   vocabulary but not the information, and that should be said out loud rather
   than implied by the upgrade in naming.

2. IS THE CONSTRUCT'S IRES ACTUALLY DISTINCTIVE?
   The only null that isolates the arrangement axis is composition_shuffle: it
   holds the codon multiset exactly fixed and randomises only the order.  Under
   that null the construct's codon arrangement is exchangeable with the 64
   shuffles, so its rank among the 65 is uniform and the one-sided p-value is
   (1 + #{shuffles at least as good}) / 65.  The pooled / weighted / uniform
   ensembles co-vary composition and cannot isolate arrangement, so they are
   reported for context only.

Outputs ``ires_analysis.json``.

Run:  python ires_analysis.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
IRES = ROOT / "ires_metrics.csv"
CAND = ROOT / "candidate_metrics.csv"
OUT = ROOT / "ires_analysis.json"

# (base column name, human label, higher_is_better or None when ambiguous)
BASES = [
    ("ires_native_recall_hard", "IRES 自身结构保持率(硬)", True),
    ("ires_native_recall_soft", "IRES 自身结构保持率(软)", True),
    ("ires_cross_occupancy", "IRES 被环外侵占率", False),
    ("ires_mean_pair_prob", "IRES 平均配对概率", True),
    ("ires_unpaired_frac", "IRES 非配对残基比例", None),
]
# window tag -> the L_IRES column this metric is meant to replace
TAGS = {"666": "lires", "669": "lires_ires669"}


def rank_p(obs, sample, higher_is_better):
    """One-sided exchangeability p-value: rank of `obs` among {obs} U sample."""
    ge = int(np.sum(sample >= obs)) if higher_is_better else int(np.sum(sample <= obs))
    n = len(sample) + 1
    return (1 + ge) / n, ge, n


def main() -> None:
    ires = pd.read_csv(IRES, encoding="utf-8-sig")
    ens = pd.read_csv(CAND, encoding="utf-8-sig")[
        ["candidate_id", "ensemble", "lires", "lires_ires669"]
    ]
    df = ires.merge(ens, on="candidate_id", how="left", validate="one_to_one")
    if df["ensemble"].isna().any():
        raise SystemExit("some candidates have no ensemble label")
    o = df[df.candidate_id == "original"].iloc[0]
    sh = df[df.ensemble == "composition_shuffle"]
    print(f"{len(df)} candidates; shuffle n={len(sh)}")

    subsets = {
        "pooled 256": df[df.ensemble.isin(["human_weighted", "uniform", "composition_shuffle"])],
        "human_weighted 96": df[df.ensemble == "human_weighted"],
        "uniform 96": df[df.ensemble == "uniform"],
        "composition_shuffle 64": sh,
        "weighted+uniform 192": df[df.ensemble.isin(["human_weighted", "uniform"])],
    }

    percentiles: dict = {}
    exchange: dict = {}
    for tag in TAGS:
        print("\n" + "=" * 96)
        print(f"IRES 窗口 {tag} nt  ——  原序列实测值")
        print("=" * 96)
        for base, label, _ in BASES:
            print(f"    {label:<30} {o[f'{base}_{tag}']:.4f}")
        print(f"\n    {'指标':<30}" + "".join(f"{k:>19}" for k in subsets))
        for base, label, hib in BASES:
            col = f"{base}_{tag}"
            cells = []
            for name, sub in subsets.items():
                v = sub[col].dropna().to_numpy(float)
                obs = float(o[col])
                if hib is None or len(v) == 0:
                    cells.append(f"{'—':>19}")
                    continue
                worse = (v < obs).sum() if hib else (v > obs).sum()
                pct = 100.0 * (worse + 0.5 * (v == obs).sum()) / len(v)
                cells.append(f"{pct:>18.1f}%")
                percentiles[f"{col}|{name}"] = float(pct)
            print(f"    {label:<30}" + "".join(cells))

        print(f"\n    【关键】可交换性检验（组成匹配 shuffle，秩 1..65）")
        for base, label, hib in BASES:
            col = f"{base}_{tag}"
            v = sh[col].dropna().to_numpy(float)
            obs = float(o[col])
            if hib is None:
                print(f"    {label:<30} 原={obs:8.4f}  组内均值={v.mean():8.4f}   "
                      f"(方向不明,仅描述)")
                exchange[col] = {"original": obs, "shuffle_mean": float(v.mean()),
                                 "note": "direction ambiguous -- descriptive only"}
                continue
            p, ge, n = rank_p(obs, v, hib)
            sd = v.std(ddof=1)
            z = (obs - v.mean()) / sd if sd > 0 else float("nan")
            mark = "★ 显著" if p < 0.05 else "  不显著"
            print(f"    {label:<30} 原={obs:8.4f}  均值={v.mean():8.4f} SD={sd:6.4f}  "
                  f"z={z:+5.2f}  秩={ge + 1:2d}/65  p={p:.3f}  {mark}")
            exchange[col] = {"original": obs, "shuffle_mean": float(v.mean()),
                             "shuffle_sd": float(sd), "z": float(z), "rank": ge + 1,
                             "n": n, "p_exchangeability": float(p),
                             "significant_at_0.05": bool(p < 0.05)}

    print("\n" + "=" * 96)
    print("新指标 vs 被替换的 L_IRES（Spearman，n=261）")
    print("=" * 96)
    corr = {}
    for tag, old in TAGS.items():
        for base, label, _ in BASES[:3]:
            c = df[[f"{base}_{tag}", old]].dropna()
            rho, pv = stats.spearmanr(c.iloc[:, 0], c.iloc[:, 1])
            corr[f"{base}_{tag}_vs_{old}"] = {"spearman_rho": float(rho), "p": float(pv)}
            print(f"  {label:<30} vs {old:<14} rho={rho:+.3f}  p={pv:.2e}")
    print("  （|rho| 接近 1 = 换了说法没换信息；接近 0 = 确实测到了不同的东西）")

    OUT.write_text(json.dumps({
        "what": "new IRES metrics under correct null models",
        "percentiles_by_null": percentiles,
        "exchangeability_test": exchange,
        "vs_replaced_lires": corr,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
