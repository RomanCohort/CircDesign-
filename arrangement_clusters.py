#!/usr/bin/env python3
"""Do the designs that beat the construct cluster, or are they scattered tail draws?

THE QUESTION
------------
Across two pre-registered samples, 30 of 1,553 buildable synonymous rearrangements
dominate the construct, and they are more folded than the background (Fisher
p = 0.0197). Two things follow that are worth separating:

  * if the dominators form COHERENT GROUPS in arrangement space, there is a
    describable class of arrangements that does better -- a mechanism, and possibly
    a designable principle;
  * if they are SCATTERED, they are simply the tail of a smooth distribution and
    there is nothing to design against.

THE EMBEDDING, AND WHY NOT A LANGUAGE MODEL
-------------------------------------------
No learned RNA model is used, deliberately. The arrangement of a coding sequence is
already represented exactly by its codon-pair spectrum, and the shuffle null has an
exact expectation for that spectrum under a uniformly random permutation of the same
codon multiset:

    E[C_ij] = n_i * (n_j - [i == j]) / N

so the residual R = C - E is, by construction, "how far this arrangement sits from a
random arrangement of exactly these codons" -- the linearisation of the same null the
whole report uses. It is 3,721-dimensional, sequence-only, and reproducible by anyone
in seconds, which a black-box embedding is not. A learned embedding would be added
only if it answered a question this does not; it does not.

ANALYSES
--------
  1. PCA of the residual; the construct's position and its neighbourhood density.
  2. k-means over a range of k; for each k, a contingency test of whether dominators
     are ENRICHED in particular clusters (chi-square), reported for every k so the
     range cannot be cherry-picked.
  3. Correlation of the leading components with CDS folding degree, to check whether
     this geometry is the same finding as the folding result or a different one.
  4. The dominators' mutual spread, compared with the spread of a random subset of the
     same size -- if the dominators are tighter than random, they are a group.

Outputs arrangement_clusters.json and arrangement_clusters.csv.

Run:  python arrangement_clusters.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

import accessibility as A
from advantage_space import makeup_shuffles
from axis_decomposition import (
    SENSE,
    SENSE_IDX,
    pair_expectation,
    read_fasta,
    sense_codons,
)

ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "arrangement_clusters.json"
OUT_CSV = ROOT / "arrangement_clusters.csv"

CON_MFE = -779.4000244140625
CON_REC = 0.7956989247311828
N_PCS = 20
K_RANGE = (2, 3, 4, 5, 6, 8, 10, 12)


def residual(cds: str) -> np.ndarray:
    cod = sense_codons(cds)
    cnt = np.zeros(len(SENSE))
    pair = np.zeros((len(SENSE), len(SENSE)))
    for c in cod:
        cnt[SENSE_IDX[c]] += 1
    for k in range(len(cod) - 1):
        pair[SENSE_IDX[cod[k]], SENSE_IDX[cod[k + 1]]] += 1
    return (pair - pair_expectation(cnt)).ravel()


def main() -> None:
    recs = dict(read_fasta(ROOT / "synonymous_candidate_CDS.fasta"))
    orig = recs["original"]

    frames = []
    for s in (771, 4242):
        d = pd.read_csv(ROOT / f"advantage_space_seed{s}.csv", encoding="utf-8-sig")
        d["seed"] = s
        d["pos"] = range(len(d))
        frames.append(d)
    pool = pd.concat(frames, ignore_index=True)
    pool["domin"] = (pool.mfe <= CON_MFE) & (pool.recall >= CON_REC)

    # accessibility / folding degree, where computed
    acc_path = ROOT / "accessibility.csv"
    fold = {}
    if acc_path.exists():
        acc = pd.read_csv(acc_path, encoding="utf-8-sig")
        for _, r in acc.iterrows():
            fold[r.candidate_id] = 1.0 - r["acc_CDS"] if "acc_CDS" in acc.columns else np.nan

    print("computing arrangement residuals (sequence only, no folding)...")
    seqs = {s: makeup_shuffles(orig, 2000, 20260804 + int(s))
            for s in sorted(pool["seed"].unique())}
    rows, meta = [], []
    for _, r in pool.iterrows():
        cds = seqs[int(r.seed)][int(r.pos)]
        rows.append(residual(cds))
        meta.append(r)
    R = np.vstack(rows)
    X_orig = residual(orig)
    meta = pd.DataFrame(meta).reset_index(drop=True)
    print(f"  residual matrix {R.shape}; construct residual computed")

    Z = PCA(n_components=N_PCS, random_state=0).fit_transform(np.vstack([R, X_orig]))
    Z_all, z_orig = Z[:-1], Z[-1]
    ev = PCA(n_components=N_PCS, random_state=0).fit(np.vstack([R, X_orig])
                                                     ).explained_variance_ratio_

    dom_idx = meta.index[meta.domin].to_numpy()
    build_idx = meta.index[meta.correct.astype(bool)].to_numpy()
    dom_build = np.array([i for i in dom_idx if i in set(build_idx.tolist())])
    print(f"  dominators {len(dom_idx)} (buildable {len(dom_build)}), buildable total {len(build_idx)}")

    # ---- 1. is the construct an outlier? ---------------------------------------
    d_con = np.linalg.norm(Z_all - z_orig, axis=1)
    pct_near = 100.0 * (d_con < np.median(d_con)).mean()
    print(f"\n[1] 构造在排列空间的位置")
    print(f"    到全体候选的距离：中位 {np.median(d_con):.2f}，构造处于第 {pct_near:.1f} 百分位")
    nn = np.sort(d_con)[:20]
    print(f"    最近 20 个邻居的距离：{nn.min():.2f} – {nn.max():.2f}")

    # ---- 2. k-means enrichment scan ---------------------------------------------
    print(f"\n[2] k-means 富集扫描（每一 k 都报告，不做挑选）")
    print(f"    {'k':>3}{'卡方 p':>12}{'支配者最富集簇':>16}{'该簇支配率':>12}{'整体支配率':>12}")
    scan = {}
    for k in K_RANGE:
        km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(Z_all)
        lab = km.labels_
        tab = pd.crosstab(pd.Series(lab, name="cluster"),
                          pd.Series(meta.domin.values, name="dom"))
        if tab.shape[1] < 2:
            continue
        chi2, pv, _, _ = stats.chi2_contingency(tab)
        rates = {}
        for c in tab.index:
            n_c = int(tab.loc[c].sum())
            k_c = int(tab.loc[c].get(True, 0))
            rates[c] = (k_c, n_c, 100.0 * k_c / n_c if n_c else 0.0)
        best = max(rates.items(), key=lambda kv: kv[1][2])
        overall = 100.0 * meta.domin.mean()
        scan[k] = {"chi2_p": float(pv), "overall_rate_percent": overall,
                   "clusters": {str(c): {"dominators": v[0], "n": v[1], "rate_percent": v[2]}
                                for c, v in rates.items()}}
        print(f"    {k:>3}{pv:>12.4f}{best[0]:>16}{best[1][2]:>11.1f}%{overall:>11.2f}%")

    # ---- 3. which component carries folding? ------------------------------------
    print(f"\n[3] 主成分与折叠程度的关系")
    corr = {}
    fvals = np.array([fold.get(f"s{int(m.seed)}p{int(m.pos)}", np.nan)
                      for _, m in meta.iterrows()], dtype=float)
    if np.isfinite(fvals).sum() > 30:
        for i in range(4):
            m = np.isfinite(fvals)
            rho, pv = stats.spearmanr(Z_all[m, i], fvals[m])
            corr[f"PC{i+1}"] = {"spearman_rho": float(rho), "p": float(pv)}
            print(f"    PC{i+1}  vs CDS 折叠程度  rho={rho:+.3f}  p={pv:.2e}")
    else:
        print("    (accessibility.csv 覆盖的候选太少，跳过)")

    # ---- 4. are the dominators tighter than random? -----------------------------
    print(f"\n[4] 支配者彼此是否比随机更聚集？")
    if len(dom_idx) >= 5:
        def spread(idx):
            sub = Z_all[idx]
            c = sub.mean(0)
            return float(np.mean(np.linalg.norm(sub - c, axis=1)))
        obs = spread(dom_idx)
        rng = np.random.default_rng(20260804)
        null = np.array([spread(rng.choice(len(Z_all), size=len(dom_idx), replace=False))
                         for _ in range(2000)])
        pv = float((null <= obs).mean())
        z = (obs - null.mean()) / null.std(ddof=1)
        print(f"    支配者平均离散度 {obs:.2f}；随机同尺寸子集均值 {null.mean():.2f} "
              f"(SD {null.std(ddof=1):.2f})")
        print(f"    z = {z:+.2f}   单尾 p(支配者更紧) = {pv:.4f}  "
              f"{'★ 显著更聚集' if pv < 0.05 else '不显著'}")
        spread_res = {"observed": obs, "null_mean": float(null.mean()),
                      "null_sd": float(null.std(ddof=1)), "z": float(z), "p_tighter": pv}
    else:
        spread_res = None
        print("    支配者太少")

    scores = meta[["seed", "pos", "mfe", "recall", "correct", "domin"]].copy()
    for i in range(6):
        scores[f"PC{i+1}"] = Z_all[:, i]
    scores["dist_to_construct"] = d_con
    scores.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

    OUT_JSON.write_text(json.dumps({
        "what": "arrangement-space geometry of the dominators",
        "embedding": "codon-pair count matrix minus its exact expectation under a "
                     "uniformly random permutation of the same codon multiset",
        "n_candidates": int(len(Z_all)),
        "n_dominators": int(len(dom_idx)),
        "construct_position": {
            "median_distance": float(np.median(d_con)),
            "percentile_of_distance": float(pct_near),
            "nearest_20_range": [float(nn.min()), float(nn.max())],
        },
        "explained_variance_first_20": [float(x) for x in ev],
        "cluster_scan": scan,
        "component_vs_folding": corr,
        "dominator_spread": spread_res,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT_CSV.name} and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
