#!/usr/bin/env python3
"""Out-of-sample validation of the folding-degree finding, on a THIRD seed.

WHAT IS BEING TESTED
--------------------
On the two pre-registered samples (seed offsets 771 and 4242), adding CDS folding
degree (1 - acc_CDS) to the objective set reduced the domination rate from 1.93% to
0.97%, and the reduction PASSED the preferential test -- dominators passed at 50.0%
against a background rate of 72.0%, Fisher p = 0.0197. Neither IRES occlusion nor the
accessibility panel reached that; two of five tested combinations were significant
against an expected 0.25 false positives.

That result is therefore suggestive but not established: it came from a search over
combinations, and five tests were run. This script draws a COMPLETELY NEW sample with
an unused seed offset and re-runs the same preferential test, with the direction of
the prediction fixed in advance: if folding degree is a real discriminator, the
dominators of this new sample should also pass it at a lower rate than the new
sample's background.

STAGING
-------
MFE for all 2000 (cheap), then the partition function only for (a) MFE-competitive
candidates -- needed for the IRES recall that defines dominance -- and (b) a random
background of buildable designs, needed for the preferential comparison. The same
partition function yields the accessibility and folding quantities.

PRE-REGISTERED PREDICTION
-------------------------
folding degree (1 - acc_CDS) removes dominators preferentially at p < 0.05.

Outputs validate_seed_<offset>.json.

Run:  python validate_seed.py [--offset 9999] [--n 2000] [--workers 28]
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
import combined_circdesign_analysis as base
from advantage_space import makeup_shuffles

ROOT = Path(__file__).resolve().parent
CUTOFF = base.DEFAULT_BPP_CUTOFF
CON_MFE = -779.4000244140625
CON_REC = 0.7956989247311828
CDS_MATURE = base.relative_interval(base.CDS_START, base.CDS_END)
IRES_MATURE = base.relative_interval(base.IRES_START, base.IRES_END_PRIMARY)

_REF: dict = {}


def _init(ref_pairs):
    global _REF
    _REF = {"pairs": ref_pairs}


def fold_all(item):
    """MFE only -- cheap screen. Receives a SPLICED mature-circle RNA, not a CDS.

    The splice happens in the parent (see main) because multiprocessing spawn on
    Windows gives each worker a fresh interpreter with an empty _MATURE_RNA template.
    """
    idx, seq = item
    fc = base.RNA.fold_compound(seq, base.vienna_model(True), base.RNA.OPTION_MFE)
    return idx, float(fc.mfe()[1])


def fold_full(item):
    """MFE + partition function: IRES recall, accessibility, folding degree.

    Also receives a spliced mature-circle RNA.
    """
    idx, seq = item
    fc = base.RNA.fold_compound(seq, base.vienna_model(True),
                                base.RNA.OPTION_MFE | base.RNA.OPTION_PF)
    dbn, mfe = fc.mfe()
    fc.pf()
    pl = {(int(it.i), int(it.j)): float(it.p) for it in fc.plist_from_probs(CUTOFF)}

    ref = _REF["pairs"]
    recall = (sum(1 for p in ref if pl.get(p, 0.0) >= 0.5) / len(ref)) if ref else float("nan")

    n = len(seq)
    paired = np.zeros(n + 1)
    for (i, j), p in pl.items():
        paired[i] += p
        paired[j] += p
    acc = 1.0 - paired
    cds_lo, cds_hi = CDS_MATURE
    ires_lo, ires_hi = IRES_MATURE
    return {
        "idx": idx, "mfe": float(mfe), "recall": recall,
        "longest_helix": float(A.longest_stack(dbn)),
        "acc_TIR": float(acc[cds_lo : cds_lo + 29].mean()),
        "acc_SPACER": float(acc[ires_hi + 1 : cds_lo].mean()),
        "acc_CDS": float(acc[cds_lo : cds_hi + 1].mean()),
        "acc_IRES": float(acc[ires_lo : ires_hi + 1].mean()),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offset", type=int, default=9999, help="seed offset; must be unused")
    ap.add_argument("--n", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    if args.offset in (771, 4242):
        raise SystemExit("offset already used -- this must be an OUT-OF-SAMPLE draw")

    A._MATURE_RNA = "".join(
        l.strip() for l in A.MATURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("T", "U")
    recs = A.read_fasta(A.CDS_FASTA)

    ref, _info = base.constrained_ires_reference(
        A._MATURE_RNA.replace("T", "U"), IRES_MATURE, CUTOFF)
    ref_hard = sorted(p for p, v in ref.items() if v >= 0.5)
    if not ref_hard:
        raise SystemExit("reference IRES pair set is empty -- refusing to run")
    print(f"reference IRES pairs (p>=0.5): {len(ref_hard)}")

    print(f"[stage A] generating {args.n} shuffles at NEW seed offset {args.offset}...")
    seqs = makeup_shuffles(recs["original"], args.n, 20260804 + args.offset)
    print(f"  {len(seqs)} distinct sequences")

    print("[stage B] correctness QC...")
    qc = [A.qc_fast(s) if hasattr(A, "qc_fast") else None for s in seqs]
    if qc[0] is None:
        from advantage_space import qc_fast
        qc = [qc_fast(s) for s in seqs]
    df = pd.DataFrame(qc)
    df["correct"] = (df.splice_donors == 0) & (df.polya_signals == 0)
    print(f"  correctness-passing: {int(df.correct.sum())}/{len(df)}")

    # Splice in the PARENT. Workers must never rebuild the sequence from a module
    # global -- spawn gives them an empty template and they would fold the bare CDS,
    # which returns a plausible energy (-431 vs -779) and corrupts the whole run.
    mature = [A.build_mature(s) for s in seqs]
    if {len(m) for m in mature} != {len(A._MATURE_RNA)}:
        raise SystemExit(f"guard failed: spliced lengths {set(len(m) for m in mature)}")
    fc0 = base.RNA.fold_compound(A.build_mature(recs["original"]),
                                 base.vienna_model(True), base.RNA.OPTION_MFE)
    probe_mfe = float(fc0.mfe()[1])
    if abs(probe_mfe - CON_MFE) > 1e-6:
        raise SystemExit(f"guard failed: spliced construct folds to {probe_mfe}, "
                         f"recorded {CON_MFE} -- wrong molecule, refusing to run")
    print(f"  guard passed: spliced construct reproduces MFE {probe_mfe:.4f}")

    print(f"[stage C1] MFE for all {len(df)}...")
    with Pool(args.workers) as p:
        mfe = dict(p.imap_unordered(fold_all, list(enumerate(mature)), chunksize=4))
    df["mfe"] = [mfe[i] for i in range(len(df))]
    comp = df.index[df.mfe <= CON_MFE].tolist()
    print(f"  MFE-competitive: {len(comp)}/{len(df)}  "
          f"(median MFE {df.mfe.median():.2f})")

    buildable = df.index[df.correct].tolist()
    rng = np.random.default_rng(20260804)
    bg = list(rng.choice(buildable, size=min(250, len(buildable)), replace=False))
    targets = sorted(set(comp) | set(bg))
    print(f"[stage C2] partition function for {len(targets)} "
          f"({len(comp)} competitive + background, union)")
    with Pool(args.workers, initializer=_init, initargs=(ref_hard,)) as p:
        full = p.map(fold_full, [(i, mature[i]) for i in targets], chunksize=2)
    F = pd.DataFrame(full).set_index("idx")
    # fold_full also returns mfe, which the cheap screen already produced. Use the
    # overlap as a free consistency check before dropping it, rather than silently
    # suffixing the columns.
    both = df.join(F[["mfe"]], how="inner", rsuffix="_full")
    dev = float((both["mfe"] - both["mfe_full"]).abs().max())
    print(f"  guard: MFE from the screen vs from the full fold, max |dev| = {dev:.6f}")
    if dev > 1e-6:
        raise SystemExit("guard failed: the two fold passes disagree")
    F = F.drop(columns=["mfe"])
    df = df.join(F, how="left")

    df["folding_CDS"] = 1 - df["acc_CDS"]
    df["ires_occlusion"] = 1 - df["acc_IRES"]
    # worker() takes a SPLICED maturity-circle RNA, and it returns acc_* but not the
    # two derived columns, so both must be fixed here. Passing recs["original"] (a CDS)
    # would silently measure the wrong molecule again.
    o = A.worker(("original", A.build_mature(recs["original"])))
    o["folding_CDS"] = 1 - o["acc_CDS"]
    o["ires_occlusion"] = 1 - o["acc_IRES"]
    df.loc["ORIG"] = {"mfe": o["mfe"], "recall": CON_REC, "longest_helix": o["longest_helix"],
                      "acc_TIR": o["acc_TIR"], "acc_SPACER": o["acc_SPACER"],
                      "acc_CDS": o["acc_CDS"], "acc_IRES": o["acc_IRES"],
                      "folding_CDS": o["folding_CDS"], "ires_occlusion": o["ires_occlusion"],
                      "correct": True}

    dom = df[(df.index != "ORIG") & (df.mfe <= CON_MFE) & (df.recall >= CON_REC)]
    dom_ok = dom[dom.correct.astype(bool)]
    bgd = df.loc[[i for i in bg if i != "ORIG"]]
    bg_ok = bgd[bgd.correct.astype(bool)]

    print("\n" + "=" * 88)
    print(f"样本外验证（种子偏移 {args.offset}，全新样本）")
    print("=" * 88)
    print(f"  三目标支配者 {len(dom)} 条，其中通过正确性质控 {len(dom_ok)} 条")
    print(f"  背景样本 {len(bg_ok)} 条可建成设计")

    n_all = len(df) - 1
    print(f"  三目标支配率 {100 * len(dom_ok) / max(int(df.correct.sum()), 1):.2f}%")

    tests = {}
    print(f"\n  {'额外目标':<20}{'支配者通过':>12}{'背景通过':>11}{'Fisher p':>11}  预注册预测")
    for col, lab, hi in (("folding_CDS", "折叠程度", False),
                         ("ires_occlusion", "IRES 遮蔽", False)):
        kd = int((dom_ok[col] <= o[col] if not hi else dom_ok[col] >= o[col]).sum())
        kb = int((bg_ok[col] <= o[col] if not hi else bg_ok[col] >= o[col]).sum())
        _, pv = stats.fisher_exact([[kd, len(dom_ok) - kd], [kb, len(bg_ok) - kb]])
        tests[col] = {"dominators_pass": kd, "n_dominators": len(dom_ok),
                      "background_pass": kb, "n_background": len(bg_ok),
                      "fisher_p": float(pv), "significant": bool(pv < 0.05)}
        mark = "★ 复现" if (col == "folding_CDS" and pv < 0.05) else (
            "✗ 未复现" if col == "folding_CDS" else "")
        print(f"  {lab:<20}{kd:>7}/{len(dom_ok):<4}{kb:>6}/{len(bg_ok):<4}"
              f"{pv:>11.4f}  {mark}")

    fd = tests["folding_CDS"]
    verdict = ("REPRODUCED -- folding degree preferentially removes dominators in a "
               "fresh sample" if fd["significant"] else
               "NOT REPRODUCED -- the earlier significant result was a search artefact")
    print(f"\n  预注册预测：折叠程度在样本外仍显著偏好剔除支配者")
    print(f"  → {verdict}")

    (ROOT / f"validate_seed_{args.offset}.json").write_text(json.dumps({
        "what": "out-of-sample replication of the folding-degree preferential test",
        "seed_offset": args.offset,
        "n_generated": n_all,
        "n_correctness_passing": int(df.correct.astype(bool).sum()),
        "n_mfe_competitive": len(comp),
        "dominators_three_objective": int(len(dom)),
        "dominators_correctness_passing": int(len(dom_ok)),
        "n_background": int(len(bg_ok)),
        "domination_rate_percent": 100 * len(dom_ok) / max(int(df.correct.astype(bool).sum()), 1),
        "tests": tests,
        "prediction": "folding degree removes dominators preferentially at p < 0.05",
        "verdict": verdict,
        "note": ("the in-sample result came from a search over five combinations; this "
                 "run uses an unused seed and a prediction fixed before the draw"),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote validate_seed_{args.offset}.json")


if __name__ == "__main__":
    main()
