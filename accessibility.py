#!/usr/bin/env python3
"""RNA accessibility of the construct versus the designs that dominated it.

WHY THIS EXISTS
---------------
The pre-registered multi-objective test used three objectives (circular MFE, human
CAI, IRES preservation). Under it, 30 of 1,553 buildable synonymous rearrangements
dominate the construct. Two candidate explanations for discarding them were tested
and both failed: an effect-size threshold is self-defeating (the construct's own
advantage over a random arrangement, ~1.4 SD, is the same size as the dominators'
advantage over it, up to 1.19 SD), and the dominators show no systematic difference
on any quality-control metric (all Mann-Whitney p > 0.45).

A fourth objective was then specified by the project BEFORE this analysis was run:
RNA ACCESSIBILITY. It is a-priori justified for a circRNA vaccine -- the CDS has to be
readable by the ribosome, and long duplexes are what trigger PKR. It is therefore not
a post-hoc filter; it is an objective the original three-objective front omitted.

PRE-REGISTERED DEFINITIONS (fixed before any result was seen)
------------------------------------------------------------
per-residue accessibility   a_i = 1 - sum_j P(i, j), from the circular partition
                            function at the same ViennaRNA settings as the rest of
                            the project (circular, 37 C, dangles=2, noLP off)
window accessibility        mean of a_i over the window
  TIR      mature 782-811   first 10 codons of the CDS          HIGHER is better
  SPACER   mature 767-781   between the IRES and the CDS        HIGHER is better
  CDS      mature 782-1864  the whole coding sequence           HIGHER is better
  IRES     mature 101-766   descriptive only (an IRES must fold; accessibility
                            there has no declared direction)
longest helix              longest contiguous stack in the circular MFE dot-bracket
                            LOWER is better (long duplexes drive PKR)
structured fraction        fraction of CDS residues with a_i < 0.5

The paired statement of the two directions is not decorative: if the construct wins
on accessibility it also wins on structured fraction by construction, and both are
reported so that cannot be presented as two independent findings.

SAMPLE
------
construct + all dominators from both pre-registered runs + a random sample of
buildable (correctness-passing) designs to locate the construct in the population.
The random sample is drawn with a fixed seed so the draw is reproducible.

Outputs accessibility.csv and accessibility.json.

Run:  python accessibility.py [--n-background 250] [--workers 28]
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

import combined_circdesign_analysis as base
from advantage_space import makeup_shuffles

ROOT = Path(__file__).resolve().parent
MATURE_FASTA = ROOT / "combined_mature_circRNA.fasta"
CDS_FASTA = ROOT / "synonymous_candidate_CDS.fasta"
OUT_CSV = ROOT / "accessibility.csv"
OUT_JSON = ROOT / "accessibility.json"

CUTOFF = base.DEFAULT_BPP_CUTOFF
CON_MFE = -779.4000244140625
CON_REC = 0.7956989247311828
CDS_MATURE = base.relative_interval(base.CDS_START, base.CDS_END)          # (782, 1864)
IRES_MATURE = base.relative_interval(base.IRES_START, base.IRES_END_PRIMARY)  # (101, 766)

WINDOWS = {
    "TIR": (CDS_MATURE[0], CDS_MATURE[0] + 29),
    "SPACER": (IRES_MATURE[1] + 1, CDS_MATURE[0] - 1),
    "CDS": CDS_MATURE,
    "IRES": IRES_MATURE,
}
DIRECTION = {"acc_TIR": +1, "acc_SPACER": +1, "acc_CDS": +1, "acc_IRES": 0,
             "longest_helix": -1, "structured_frac_CDS": -1}

_MATURE_RNA = ""
_CDS_POS = None


def read_fasta(path: Path):
    recs, name, buf = {}, None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                recs[name] = "".join(buf)
            name, buf = line[1:].split("|")[0], []
        else:
            buf.append(line)
    if name is not None:
        recs[name] = "".join(buf)
    return recs


def build_mature(cds_dna: str) -> str:
    # Hard refusal rather than a silent fallback. On Windows, multiprocessing spawns
    # fresh interpreters, so a worker that calls this without the parent having set
    # _MATURE_RNA gets an EMPTY template and the function returns the bare 1083 nt CDS.
    # The fold still succeeds and returns a plausible-looking energy (-431 instead of
    # -779), which has silently corrupted one run of advantage_space.py and one run of
    # validate_seed.py. Refusing here catches the mistake at the source, in every
    # caller, instead of relying on each script to remember the rule.
    if len(_MATURE_RNA) < 1000:
        raise RuntimeError(
            "build_mature called with an unset/short _MATURE_RNA template "
            f"(len={len(_MATURE_RNA)}). In a multiprocessing worker this means the "
            "parent did not splice the sequence -- the worker would fold the bare CDS. "
            "Splice in the parent and pass the finished sequence to the worker."
        )
    lo, hi = CDS_MATURE
    out = _MATURE_RNA[: lo - 1] + cds_dna.replace("T", "U") + _MATURE_RNA[hi:]
    if len(out) != len(_MATURE_RNA):
        raise RuntimeError("build_mature produced a sequence of the wrong length")
    return out


def longest_stack(dbn: str) -> int:
    """Longest run of nested, contiguously stacked '(' ... ')' pairs."""
    best = 0
    stack: list[int] = []
    depth_pairs: list[tuple[int, int]] = []
    for idx, ch in enumerate(dbn):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if stack:
                depth_pairs.append((stack.pop(), idx))
    pairs = sorted(depth_pairs)
    run = 1
    for (i1, j1), (i2, j2) in zip(pairs, pairs[1:]):
        if i2 == i1 + 1 and j2 == j1 - 1:
            run += 1
            best = max(best, run)
        else:
            run = 1
    return max(best, 1 if pairs else 0)


def worker(item):
    """Measure one already-spliced mature-circle RNA.

    The splice is done in main(), NOT here. This function used to call build_mature()
    itself, which reads the module-level _MATURE_RNA template; under Windows spawn each
    worker is a fresh interpreter where that template is empty, so every worker folded
    the bare 1083 nt CDS instead of the 2013 nt circle. The fold succeeded and returned
    plausible numbers (MFE -450.6 rather than -779.4), so the result looked usable and
    was not caught until the construct's own row was checked against its recorded MFE.
    """
    sid, seq = item
    fc = base.RNA.fold_compound(seq, base.vienna_model(True),
                                base.RNA.OPTION_MFE | base.RNA.OPTION_PF)
    dbn, mfe = fc.mfe()
    fc.pf()
    plist = {(int(it.i), int(it.j)): float(it.p) for it in fc.plist_from_probs(CUTOFF)}

    n = len(seq)
    paired = np.zeros(n + 1)
    for (i, j), p in plist.items():
        paired[i] += p
        paired[j] += p
    acc = 1.0 - paired  # 1-based indexing; acc[0] unused

    row = {"candidate_id": sid, "mfe": float(mfe), "longest_helix": longest_stack(dbn)}
    for tag, (lo, hi) in WINDOWS.items():
        row[f"acc_{tag}"] = float(acc[lo : hi + 1].mean())
    cds_lo, cds_hi = CDS_MATURE
    row["structured_frac_CDS"] = float((acc[cds_lo : cds_hi + 1] < 0.5).mean())
    return row


def main() -> None:
    global _MATURE_RNA
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-background", type=int, default=250)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    _MATURE_RNA = "".join(
        l.strip() for l in MATURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("T", "U")
    recs = read_fasta(CDS_FASTA)

    # --- assemble the sample -----------------------------------------------------
    frames = []
    for s in (771, 4242):
        path = ROOT / f"advantage_space_seed{s}.csv"
        if path.exists():
            d = pd.read_csv(path, encoding="utf-8-sig")
            d["seed"] = s
            d["pos"] = range(len(d))   # positional index == sample order within this seed
            frames.append(d)
    if not frames:
        raise SystemExit("no advantage_space_seed*.csv found -- run advantage_space.py first")
    pool_df = pd.concat(frames, ignore_index=True)

    dom = pool_df[(pool_df.mfe <= CON_MFE) & (pool_df.recall >= CON_REC)].copy()
    print(f"dominators available: {len(dom)}")

    buildable = pool_df[pool_df.correct]
    bg = buildable.sample(n=min(args.n_background, len(buildable)), random_state=20260804)
    print(f"background sample: {len(bg)} buildable designs (seed 20260804)")

    # advantage_space.py drops the sequence column before writing, and no identifier is
    # kept. The row order is however the sample order, and the sampler is deterministic
    # given its seed, so the sequences are regenerated by replaying it. Recorded here
    # rather than silently patched, because it is a real gap in that file's output.
    def seq_for(seed_offset, n: int = 2000):
        # int() is required: pandas hands back numpy int64, which random.Random rejects.
        return makeup_shuffles(recs["original"], n, 20260804 + int(seed_offset))

    cache = {s: seq_for(s) for s in sorted(pool_df["seed"].unique())}

    def seq_at(row):
        return cache[int(row.seed)][int(row.pos)]

    todo = [("original", build_mature(recs["original"]))]
    todo += [(f"dom_seed{int(r.seed)}_{i:03d}", build_mature(seq_at(r)))
             for i, (_, r) in enumerate(dom.iterrows())]
    todo += [(f"bg_{i:04d}", build_mature(seq_at(r))) for i, (_, r) in enumerate(bg.iterrows())]
    lens = {len(s) for _, s in todo}
    if lens != {len(_MATURE_RNA)}:
        raise SystemExit(f"guard failed: spliced lengths {lens}, expected {len(_MATURE_RNA)}")
    probe = todo[0][1]
    fc0 = base.RNA.fold_compound(probe, base.vienna_model(True),
                                 base.RNA.OPTION_MFE | base.RNA.OPTION_PF)
    pm = float(fc0.mfe()[1])
    if abs(pm - CON_MFE) > 1e-6:
        raise SystemExit(f"guard failed: spliced construct folds to {pm}, expected "
                         f"{CON_MFE} -- wrong molecule, refusing to run")
    print(f"  guard passed: spliced construct reproduces MFE {pm:.4f}")
    print(f"folding {len(todo)} sequences on {args.workers} workers...")

    with Pool(args.workers) as pool:
        rows = pool.map(worker, todo, chunksize=2)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

    o = df[df.candidate_id == "original"].iloc[0]
    D = df[df.candidate_id.str.startswith("dom_")]
    B = df[df.candidate_id.str.startswith("bg_")]

    print("\n" + "=" * 84)
    print("构造的可及性 vs 背景样本（可建成设计）")
    print("=" * 84)
    print(f"{'指标':<20}{'构造':>10}{'背景中位':>11}{'构造百分位':>12}{'方向':>8}")
    res = {}
    for col, lab in (("acc_TIR", "TIR"), ("acc_SPACER", "SPACER"), ("acc_CDS", "CDS"),
                     ("acc_IRES", "IRES"), ("longest_helix", "最长螺旋")):
        b = B[col].to_numpy(float)
        obs = float(o[col])
        d = DIRECTION[col]
        worse = int((b > obs).sum()) if d > 0 else (int((b < obs).sum()) if d < 0 else 0)
        pctl = 100.0 * worse / len(b) if d else float("nan")
        res[col] = {"construct": obs, "background_median": float(np.median(b)),
                    "background_mean": float(b.mean()), "background_sd": float(b.std(ddof=1)),
                    "percentile": pctl, "direction": d}
        dirstr = "越高越好" if d > 0 else ("越短越好" if d < 0 else "描述性")
        print(f"{lab:<20}{obs:>10.4f}{np.median(b):>11.4f}{pctl:>11.1f}%{dirstr:>8}")

    print("\n" + "=" * 84)
    print("支配者 vs 背景：可及性上有没有系统性差异？")
    print("=" * 84)
    print(f"{'指标':<20}{'支配者中位':>12}{'背景中位':>11}{'Mann-Whitney p':>16}  {'支配者更差?':>12}")
    tests = {}
    for col, lab in (("acc_TIR", "TIR"), ("acc_SPACER", "SPACER"), ("acc_CDS", "CDS"),
                     ("longest_helix", "最长螺旋")):
        a = D[col].to_numpy(float)
        b = B[col].to_numpy(float)
        try:
            _, p2 = stats.mannwhitneyu(a, b, alternative="two-sided")
            alt = "less" if DIRECTION[col] > 0 else "greater"
            _, p1 = stats.mannwhitneyu(a, b, alternative=alt)
        except ValueError:
            p2 = p1 = float("nan")
        tests[col] = {"p_two_sided": float(p2), "p_dominator_worse": float(p1)}
        print(f"{lab:<20}{np.median(a):>12.4f}{np.median(b):>11.4f}{p2:>16.4f}"
              f"{('是' if p1 < 0.05 else '否'):>12}")

    k = len(tests)
    print(f"\n多重比较：{k} 个指标 × α=0.05 → 期望假阳性 {0.05 * k:.2f} 条")

    OUT_JSON.write_text(json.dumps({
        "what": "RNA accessibility, pre-registered windows",
        "pre_registration": {
            "windows_mature_1based": {k: list(v) for k, v in WINDOWS.items()},
            "direction": DIRECTION,
            "background_sample_seed": 20260804,
            "n_background": int(len(B)),
            "n_dominators": int(len(D)),
            "mechanism": "ribosome loading and elongation; long duplexes drive PKR",
        },
        "construct_vs_background": res,
        "dominators_vs_background": tests,
        "expected_false_positives_at_0.05": 0.05 * k,
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT_CSV.name} and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
