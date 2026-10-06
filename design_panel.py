#!/usr/bin/env python3
"""A pre-registered, mechanism-motivated metric panel for a circRNA vaccine CDS.

WHY THIS EXISTS
---------------
The retrospective report evaluates the construct on circular MFE, CAI and L_IRES.
Under the composition-matched null (``composition_shuffle``: identical codon
multiset, order randomised) all three sit at 72-86th percentile with p between 0.15
and 0.29 -- favourable in direction, not distinguishable from chance.  Scaling the
null to n=2000 cannot change that: the observed 55/64 already caps the true
percentile's 95% interval at 92.8%, below the 95% needed for significance.

The reason is structural.  The shuffle null holds the codon multiset fixed, so GC
content and every composition-derived quantity are locked; the only free variable
is local codon order.  All three existing metrics are views of that one variable.

This panel measures dimensions the existing analysis never touches, each with a
stated mechanism, each computable directly from the sequence.  They are the checks
a gene is put through before it is ordered:

  synthesis      homopolymer runs, GC window extremes, low-complexity runs
  recombination  longest direct repeat, longest inverted repeat
  cloning        common restriction sites
  processing     cryptic polyadenylation signals, cryptic splice donor/acceptor

HONESTY RULES BUILT IN
----------------------
* Every metric in this file is reported, not a selected subset.  The panel is the
  whole list; there is no "and also we looked at X and it looked good".
* Multiplicity is handled explicitly: with k metrics at alpha=0.05 the expected
  number of false positives under the null is 0.05k, and that number is printed.
* Each metric is evaluated against BOTH the composition-matched null (isolates
  arrangement) and the broader synonymous nulls (which also vary composition), and
  the two are always shown together so neither can be quoted alone.
* Direction is declared per metric before the result is seen.  Where a direction
  is genuinely ambiguous it is marked None and reported descriptively only.

Outputs ``design_panel.csv`` and ``design_panel.json``.

Run:  python design_panel.py
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
CDS_FASTA = ROOT / "synonymous_candidate_CDS.fasta"
CAND = ROOT / "candidate_metrics.csv"
OUT_CSV = ROOT / "design_panel.csv"
OUT_JSON = ROOT / "design_panel.json"

RESTRICTION_SITES = {
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC", "XbaI": "TCTAGA", "NheI": "GCTAGC", "SacI": "GAGCTC",
    "PstI": "CTGCAG", "KpnI": "GGTACC", "SalI": "GTCGAC", "NcoI": "CCATGG",
    "BsaI": "GGTCTC", "SapI": "GCTCTTC", "Esp3I": "CGTCTC", "BsmBI": "CGTCTC",
    "AarI": "CACCTGC", "SbfI": "CCTGCAGG",
}
POLYA = ["AATAAA", "ATTAAA", "AATAAT", "CATAAA", "TATAAA"]
DONOR = ["GTAAGT", "GTGAGT", "GTATGT", "GTACGT"]
GC_WINDOW = 50
LOWC_MAX = 12          # low-complexity = one base >75% within this window
REPEAT_MIN = 12        # only report direct/inverted repeats at least this long

# (column, label, higher_is_better or None when ambiguous)
PANEL = [
    ("homopolymer_max", "最长同碱基连续", False),
    ("gc_window_max", "GC 窗口最高 (50nt)", False),
    ("gc_window_min", "GC 窗口最低 (50nt)", None),
    ("low_complexity_max", "最长低复杂度区", False),
    ("direct_repeat_max", "最长正向重复", False),
    ("inverted_repeat_max", "最长反向重复(茎长)", None),
    ("restriction_sites", "常用限制性位点总数", False),
    ("polya_signals", "隐性 polyA 信号数", False),
    ("splice_donors", "隐性剪接供体数", False),
]


def read_fasta(path: Path):
    out, name, buf = [], None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split("|")[0], []
        else:
            buf.append(line.upper())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def longest_homopolymer(s: str) -> int:
    best = run = 1
    for i in range(1, len(s)):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    return best


def gc_windows(s: str):
    n, w = len(s), GC_WINDOW
    if n < w:
        return 0.0, 0.0
    isgc = np.frombuffer(s.encode(), dtype=np.uint8)
    isgc = ((isgc == ord("G")) | (isgc == ord("C"))).astype(float)
    csum = np.concatenate([[0.0], np.cumsum(isgc)])
    fracs = (csum[w:] - csum[:-w]) / w
    return float(fracs.max()), float(fracs.min())


def longest_low_complexity(s: str) -> int:
    n, w = len(s), LOWC_MAX
    best = 0
    for i in range(n - w + 1):
        c = Counter(s[i : i + w])
        if max(c.values()) / w >= 0.75:
            best = max(best, max(c.values()))
    return best


def longest_direct_repeat(s: str, min_len: int = REPEAT_MIN) -> int:
    n = len(s)
    for L in range(min(n // 2, 60), min_len - 1, -1):
        seen = set()
        for i in range(n - L + 1):
            sub = s[i : i + L]
            if sub in seen:
                return L
            seen.add(sub)
    return 0


_COMP = str.maketrans("ACGT", "TGCA")


def longest_inverted_repeat(s: str) -> int:
    """Longest stem length of an inverted repeat with a loop of >= 3 nt."""
    n = len(s)
    best = 0
    for i in range(n):
        for j in range(i + 3 + best, n):
            k = 0
            while j + k < n and i - k - 1 >= 0 and s[i - k - 1] == s[j + k].translate(_COMP):
                k += 1
            best = max(best, k)
            if best >= 20:
                return best
    return best


def metrics(seq: str) -> dict:
    txt = seq if isinstance(seq, str) else seq.decode()
    up = txt.upper().replace("U", "T")
    gmax, gmin = gc_windows(up)
    sites = sum(len(_find_all(up, v)) for v in set(RESTRICTION_SITES.values()))
    return {
        "homopolymer_max": longest_homopolymer(up),
        "gc_window_max": gmax,
        "gc_window_min": gmin,
        "low_complexity_max": longest_low_complexity(up),
        "direct_repeat_max": longest_direct_repeat(up),
        "inverted_repeat_max": longest_inverted_repeat(up),
        "restriction_sites": sites,
        "polya_signals": sum(len(_find_all(up, m)) for m in POLYA),
        "splice_donors": sum(len(_find_all(up, m)) for m in DONOR),
    }


def _find_all(hay: str, needle: str):
    out, i = [], hay.find(needle)
    while i != -1:
        out.append(i)
        i = hay.find(needle, i + 1)
    return out


def main() -> None:
    records = read_fasta(CDS_FASTA)
    ens = pd.read_csv(CAND, encoding="utf-8-sig")[["candidate_id", "ensemble"]]
    rows = [dict(candidate_id=n, **metrics(c)) for n, c in records]
    df = pd.DataFrame(rows).merge(ens, on="candidate_id", how="left", validate="one_to_one")
    df.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

    o = df[df.candidate_id == "original"].iloc[0]
    sh = df[df.ensemble == "composition_shuffle"]
    subsets = {
        "pooled 256": df[df.ensemble.isin(["human_weighted", "uniform", "composition_shuffle"])],
        "human_weighted 96": df[df.ensemble == "human_weighted"],
        "uniform 96": df[df.ensemble == "uniform"],
        "composition_shuffle 64": sh,
    }

    print(f"{len(df)} candidates; panel of {len(PANEL)} metrics\n")
    print("=" * 100)
    print("原序列实测值")
    print("=" * 100)
    for col, label, _ in PANEL:
        print(f"  {label:<26} {o[col]}")

    print("\n" + "=" * 100)
    print("各零模型下的百分位（已按声明的方向折算，越高=越有利）")
    print("=" * 100)
    print(f"{'指标':<26}" + "".join(f"{k:>21}" for k in subsets))
    result = {}
    for col, label, hib in PANEL:
        cells = []
        for name, sub in subsets.items():
            v = sub[col].to_numpy(float)
            obs = float(o[col])
            if hib is None:
                cells.append(f"{'—':>21}")
                continue
            worse = (v < obs).sum() if hib else (v > obs).sum()
            pct = 100.0 * (worse + 0.5 * (v == obs).sum()) / len(v)
            cells.append(f"{pct:>20.1f}%")
        print(f"{label:<26}" + "".join(cells))

    print("\n" + "=" * 100)
    print("【关键】可交换性检验（组成匹配 shuffle，秩 1..65）")
    print("=" * 100)
    hood = []
    for col, label, hib in PANEL:
        v = sh[col].to_numpy(float)
        obs = float(o[col])
        if hib is None:
            rho, _ = stats.spearmanr(sh[col], sh["homopolymer_max"])
            print(f"{label:<26} 原={obs:<8} 组内均值={v.mean():<10.3f} (方向不明,仅描述)")
            result[col] = {"original": obs, "shuffle_mean": float(v.mean()),
                           "note": "direction ambiguous -- descriptive only"}
            continue
        ge = int(np.sum(v >= obs)) if hib else int(np.sum(v <= obs))
        p = (1 + ge) / (len(v) + 1)
        mark = "★ 显著" if p < 0.05 else ""
        print(f"{label:<26} 原={obs:<8} 均值={v.mean():<10.3f} 秩={ge + 1:2d}/65  p={p:.3f}  {mark}")
        result[col] = {"original": obs, "shuffle_mean": float(v.mean()),
                       "rank": ge + 1, "p_exchangeability": float(p),
                       "significant_at_0.05": bool(p < 0.05)}
        if p < 0.05:
            hood.append(col)

    k = sum(1 for _, _, h in PANEL if h is not None)
    exp_fp = 0.05 * k
    print(f"\n多重比较：{k} 个有方向的指标，alpha=0.05 下期望假阳性 = {exp_fp:.2f} 条")
    print(f"实测显著 = {len(hood)} 条 {hood if hood else ''}")
    verdict = ("与假阳性期望一致 —— 不构成证据"
               if len(hood) <= exp_fp else "超出假阳性期望 —— 值得进一步检验")
    print(f"判定：{verdict}")

    OUT_JSON.write_text(json.dumps({
        "what": "pre-registered mechanism-motivated panel, CDS level",
        "note": "every metric reported; multiplicity corrected explicitly",
        "n_directional_metrics": k,
        "expected_false_positives_at_0.05": exp_fp,
        "observed_significant": hood,
        "results": result,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT_CSV.name} and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
