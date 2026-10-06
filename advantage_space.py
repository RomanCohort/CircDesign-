#!/usr/bin/env python3
"""Pre-registered test: is the construct Pareto-optimal among *buildable* designs?

THE QUESTION
------------
Within the composition-matched space (identical codon multiset, order permuted) the
construct sits at the ~86th percentile on circular MFE and ~72nd on IRES
preservation.  Neither is significant.  Scaling the null cannot change that: the
observed 55/64 caps the true percentile's 95% interval at 92.8%, short of the 95%
significance needs.

But the composition-matched space contains sequences that no designer would build.
A rearrangement that wins on MFE while introducing a cryptic splice donor or a
cryptic polyadenylation signal produces a WRONG transcript -- it is not a
competitor, it is a broken design.

This script asks the narrower, answerable question:

    Among synonymous rearrangements that are CORRECT (no cryptic splice donor, no
    cryptic polyA), how many still dominate the construct on the full objective set?

PRE-REGISTRATION (fixed before the sample was generated; do not edit after a run)
--------------------------------------------------------------------------------
sample            N = 2000 fresh composition-shuffle CDS, seed offset 771 so the
                  sample does not overlap the original 64-shuffle ensemble
correctness filter  splice_donors == 0 AND polya_signals == 0
                  Chosen because these are PRODUCT-CORRECTNESS criteria (they change
                  the transcript), not manufacturability criteria (which only make
                  synthesis harder).  Applied identically to the construct and to
                  every competitor.  The construct passes it.
objectives        (circular MFE, down) and (IRES native-pair recall, up).
                  human CAI is EXACTLY invariant across the shuffle null -- verified,
                  not assumed -- so it cannot enter a dominance test here and is
                  omitted rather than silently dropped.
dominator         MFE <= construct AND recall >= construct, at least one strict.
                  Because CAI is invariant, this is the full 2-objective front.
primary endpoint  domination rate among correctness-passing variants
decision rule     SUPPORTED if <= 2 dominators among >= 500 passing variants AND the
                  95% upper bound on the rate is < 2%.
                  Otherwise NOT SUPPORTED.
                  Both the filtered and unfiltered counts are always reported, so the
                  filter's effect is visible and cannot be quoted away.

STAGING (why this is affordable)
--------------------------------
Folding with the partition function is the expensive step (~35 s per 2013 nt
sequence).  A variant whose MFE is already worse than the construct's can never
dominate it, so the partition function is computed ONLY for MFE-competitive
variants, after a cheap MFE-only screen.

Outputs ``advantage_space.json`` and ``advantage_space.csv``.

Run:  python advantage_space.py [--n 2000] [--workers 28]
"""

from __future__ import annotations

import argparse
import json
import os
import random
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import combined_circdesign_analysis as base
from design_panel import (
    DONOR,
    POLYA,
    RESTRICTION_SITES,
    gc_windows,
    longest_homopolymer,
    read_fasta,
    _find_all,
)


def qc_fast(seq: str) -> dict:
    """The pre-registered filter, plus the cheap QC numbers worth carrying along.

    Deliberately EXCLUDES the direct/inverted repeat scans: those are O(n^2)-ish
    (~3.4 s per 1083 nt sequence, so ~110 minutes for 2000 of them) and the
    pre-registered filter does not read them. They are recomputed only for the
    construct and for whatever ends up dominating it.
    """
    up = seq.upper().replace("U", "T")
    gmax, gmin = gc_windows(up)
    return {
        "splice_donors": sum(len(_find_all(up, m)) for m in DONOR),
        "polya_signals": sum(len(_find_all(up, m)) for m in POLYA),
        "homopolymer_max": longest_homopolymer(up),
        "gc_window_max": gmax,
        "gc_window_min": gmin,
        "restriction_sites": sum(len(_find_all(up, v)) for v in set(RESTRICTION_SITES.values())),
    }

ROOT = Path(__file__).resolve().parent
CDS_FASTA = ROOT / "synonymous_candidate_CDS.fasta"
MATURE_FASTA = ROOT / "combined_mature_circRNA.fasta"
OUT_JSON = ROOT / "advantage_space.json"
OUT_CSV = ROOT / "advantage_space.json"  # rebound per-run in main() to a seed-specific name

SEED_OFFSET = 771          # pre-registered; keeps this sample disjoint from the original
N_DEFAULT = 2000
WINDOW = "666"
IRES = base.relative_interval(base.IRES_START, base.IRES_END_PRIMARY)
CDS_MATURE = base.relative_interval(base.CDS_START, base.CDS_END)
CUTOFF = base.DEFAULT_BPP_CUTOFF
PAIR_THRESHOLD = 0.5

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C",
    "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S",
    "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}

_REF = {}
_MATURE_RNA = ""


def makeup_shuffles(cds: str, n: int, seed: int):
    """n distinct per-amino-acid codon permutations of the original CDS."""
    sense = [cds[i : i + 3] for i in range(0, len(cds) - 3, 3)]
    stop = cds[-3:]
    groups = defaultdict(list)
    for idx, c in enumerate(sense):
        groups[CODON_TABLE[c]].append(idx)

    rng = random.Random(seed)
    seen, out = set(), []
    attempts = 0
    while len(out) < n and attempts < n * 50:
        attempts += 1
        codons = list(sense)
        for positions in groups.values():
            vals = [codons[p] for p in positions]
            rng.shuffle(vals)
            for p, v in zip(positions, vals):
                codons[p] = v
        s = "".join(codons) + stop
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def build_mature(cds_dna: str) -> str:
    lo, hi = CDS_MATURE
    return _MATURE_RNA[: lo - 1] + cds_dna.replace("T", "U") + _MATURE_RNA[hi:]


def _mfe_worker(item):
    """Fold one already-spliced mature-circle RNA. Index in, index out.

    The splice deliberately does NOT happen here. On Windows, multiprocessing uses
    spawn: each worker is a fresh interpreter that re-imports this module, so a
    module global assigned in the parent (_MATURE_RNA) is EMPTY in the child, and
    build_mature() would silently return the bare 1083 nt CDS instead of the
    2013 nt circle. The fold still succeeds and returns a plausible number. Doing
    the splice in the parent makes that failure mode structurally impossible.
    """
    idx, seq = item
    fc = base.RNA.fold_compound(seq, base.vienna_model(True), base.RNA.OPTION_MFE)
    _, mfe = fc.mfe()
    return idx, float(mfe)


def _pf_worker(item):
    idx, seq = item
    fc = base.RNA.fold_compound(seq, base.vienna_model(True),
                                base.RNA.OPTION_MFE | base.RNA.OPTION_PF)
    fc.mfe()
    fc.pf()
    plist = {(int(it.i), int(it.j)): float(it.p) for it in fc.plist_from_probs(CUTOFF)}
    ref = _REF["pairs"]
    if not ref:
        return idx, float("nan")
    still = sum(1 for p in ref if plist.get(p, 0.0) >= PAIR_THRESHOLD)
    return idx, still / len(ref)


def _init(ref_pairs):
    global _REF
    _REF = {"pairs": ref_pairs}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=N_DEFAULT)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    # Seed-specific outputs. An earlier version wrote a single shared filename, so the
    # confirmation run silently OVERWROTE the first run's per-variant CSV and its raw
    # evidence could not be re-examined. One file per seed makes that impossible.
    global OUT_CSV, OUT_JSON
    OUT_CSV = ROOT / f"advantage_space_seed{SEED_OFFSET}.csv"
    OUT_JSON = ROOT / f"advantage_space_seed{SEED_OFFSET}.json"
    print(f"outputs: {OUT_CSV.name}, {OUT_JSON.name}")

    global _MATURE_RNA
    # NOTE: keep this in RNA (U). build_mature() converts only the CDS from DNA.
    # Converting one side and not the other folds a chimera that is neither, and
    # ViennaRNA reads T as an unknown character -- silently, with no warning.
    _MATURE_RNA = "".join(
        l.strip() for l in MATURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("T", "U")
    records = dict(read_fasta(CDS_FASTA))
    orig_cds = records["original"]

    construct = pd.read_csv(ROOT / "ires_metrics.csv", encoding="utf-8-sig")
    c = construct[construct.candidate_id == "original"].iloc[0]
    CON_MFE = float(c["circular_mfe_kcal_mol"])
    CON_RECALL = float(c[f"ires_native_recall_hard_{WINDOW}"])
    print(f"construct: MFE={CON_MFE:.4f}  IRES recall({WINDOW})={CON_RECALL:.4f}")

    # --- guards: splice round-trips, and the reproduced fold matches the record ---
    probe = build_mature(orig_cds)
    if probe != _MATURE_RNA:
        raise SystemExit("guard failed: original CDS does not round-trip to the mature circle")
    fc = base.RNA.fold_compound(probe, base.vienna_model(True), base.RNA.OPTION_MFE)
    probe_mfe = float(fc.mfe()[1])
    if abs(probe_mfe - CON_MFE) > 1e-6:
        raise SystemExit(
            f"guard failed: reproduced MFE {probe_mfe:.4f} != recorded {CON_MFE:.4f}. "
            "The spliced sequence or its alphabet is wrong; refusing to run."
        )
    print(f"guards passed: splice round-trips; reproduced MFE {probe_mfe:.4f} matches the record")

    # ---- reference IRES pair set (identical construction to L_IRES) -------------
    ref, info = base.constrained_ires_reference(_MATURE_RNA.replace("T", "U"), IRES, CUTOFF)
    ref_hard = {p: v for p, v in ref.items() if v >= PAIR_THRESHOLD}
    print(f"reference IRES pairs with p>=0.5: {len(ref_hard)} (of {info['retained_pairs']} retained)")
    # An empty reference would make every variant's recall NaN, and _pf_worker's NaN
    # path combined with `dominates = False for NaN` would then report ZERO dominators
    # and declare the construct Pareto-optimal. A measurement failure must never be
    # able to produce evidence FOR the construct, so it is refused here.
    if not ref_hard:
        raise SystemExit(
            "reference IRES pair set is empty -- every recall would be NaN and the "
            "decision rule would read that as 'no dominators'. Refusing to run."
        )

    # ---- STAGE A: sample ---------------------------------------------------------
    print(f"\n[stage A] generating {args.n} fresh shuffles (seed offset {SEED_OFFSET})...")
    shuffles = makeup_shuffles(orig_cds, args.n, 20260804 + SEED_OFFSET)
    print(f"  {len(shuffles)} distinct sequences")

    # verify CAI invariance is a property of the sampling, not an assumption
    print("  verifying codon multiset invariance...")
    from collections import Counter as _C
    base_cnt = _C(orig_cds[i : i + 3] for i in range(0, len(orig_cds) - 3, 3))
    bad = sum(1 for s in shuffles if _C(s[i : i + 3] for i in range(0, len(s) - 3, 3)) != base_cnt)
    print(f"  sequences whose sense-codon multiset differs from the construct: {bad}")
    if bad:
        raise SystemExit("shuffle construction broken -- aborting")

    # ---- STAGE B: QC panel (free) ------------------------------------------------
    print("[stage B] QC panel...")
    rows = []
    for s in shuffles:
        m = qc_fast(s)
        m["cds"] = s
        rows.append(m)
    df = pd.DataFrame(rows)
    df["correct"] = (df["splice_donors"] == 0) & (df["polya_signals"] == 0)
    print(f"  correctness-passing: {int(df.correct.sum())}/{len(df)} "
          f"({100 * df.correct.mean():.1f}%)")

    # ---- STAGE C: MFE screen -----------------------------------------------------
    # splice in the parent, once, and assert the result before anything is folded
    mature = [build_mature(s) for s in df["cds"]]
    lengths = {len(m) for m in mature}
    if lengths != {len(_MATURE_RNA)}:
        raise SystemExit(
            f"guard failed: spliced lengths {lengths}, expected {len(_MATURE_RNA)}. "
            "The worker would have folded the wrong molecule."
        )
    if mature[0] != _MATURE_RNA and shuffles[0] == orig_cds:
        raise SystemExit("guard failed: spliced original is not the mature circle")
    print(f"[stage C] MFE for all {len(df)} on {args.workers} workers "
          f"(spliced in parent; all {len(_MATURE_RNA)} nt)...")
    t0 = time.time()
    with Pool(args.workers) as pool:
        mfes = dict(pool.imap_unordered(_mfe_worker, list(enumerate(mature)), chunksize=4))
    df["mfe"] = [mfes[i] for i in range(len(df))]
    print(f"  done in {time.time() - t0:.0f}s")
    print(f"  shuffle MFE: min={df.mfe.min():.2f} median={df.mfe.median():.2f} "
          f"max={df.mfe.max():.2f}   (construct {CON_MFE:.2f})")
    comp = df[df.mfe <= CON_MFE]
    comp_idx = comp.index.tolist()
    print(f"  MFE-competitive (mfe <= construct): {len(comp)}/{len(df)}")
    print(f"  competitive AND correctness-passing: {int(comp.correct.sum())}")
    if len(comp) == 0:
        raise SystemExit(
            "guard failed: zero MFE-competitive shuffles. With n=2000 and an expected "
            "~14% rate this cannot happen -- treat it as a bug, not a result."
        )

    # ---- STAGE D: partition function only where it can matter --------------------
    print(f"[stage D] IRES recall for {len(comp)} competitive variants...")
    t1 = time.time()
    with Pool(args.workers, initializer=_init, initargs=(sorted(ref_hard),)) as pool:
        rec = dict(pool.imap_unordered(
            _pf_worker, [(i, mature[i]) for i in comp_idx], chunksize=2))
    df["recall"] = np.nan
    for i, r in rec.items():
        df.at[i, "recall"] = r
    print(f"  done in {time.time() - t1:.0f}s")

    df["dominates"] = (df.mfe <= CON_MFE) & (df.recall >= CON_RECALL)
    df.loc[df.recall.isna(), "dominates"] = False

    # ---- analysis -----------------------------------------------------------------
    def wilson(k, n, z=1.96):
        if n == 0:
            return 0.0, 0.0
        p = k / n
        d = 1 + z * z / n
        ctr = (p + z * z / (2 * n)) / d
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
        return max(0.0, ctr - half) * 100, min(1.0, ctr + half) * 100

    print("\n" + "=" * 84)
    print("PRE-REGISTERED ENDPOINT")
    print("=" * 84)
    total = len(df)
    k_all = int(df.dominates.sum())
    feas = df[df.correct]
    k_f = int(feas.dominates.sum())
    lo_a, hi_a = wilson(k_all, total)
    lo_f, hi_f = wilson(k_f, len(feas))
    print(f"  all variants        : {k_all:4d}/{total:<5d} = {100*k_all/total:6.2f}%  95%CI [{lo_a:.2f}, {hi_a:.2f}]")
    print(f"  correctness-passing : {k_f:4d}/{len(feas):<5d} = {100*k_f/max(len(feas),1):6.2f}%  95%CI [{lo_f:.2f}, {hi_f:.2f}]")

    supported = (k_f <= 2) and (len(feas) >= 500) and (hi_f < 2.0)
    print(f"\n  decision rule: <=2 dominators AND >=500 passing AND upper bound <2%")
    print(f"  -> {'SUPPORTED' if supported else 'NOT SUPPORTED'}")
    print(f"\n  construct MFE percentile among passing variants: "
          f"{100*(feas.mfe < CON_MFE).mean():.1f}%")
    print(f"  construct recall percentile among passing variants: "
          f"{100*(feas.recall < CON_RECALL).mean():.1f}%")

    df.drop(columns=["cds"]).to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    OUT_JSON.write_text(json.dumps({
        "pre_registration": {
            "n_requested": args.n, "seed_offset": SEED_OFFSET,
            "correctness_filter": "splice_donors == 0 and polya_signals == 0",
            "objectives": ["circular MFE (down)", "IRES native recall (up)"],
            "decision_rule": "<=2 dominators among >=500 passing AND upper 95% bound < 2%",
        },
        "construct": {"mfe": CON_MFE, "ires_recall": CON_RECALL,
                      "reference_ires_pairs": len(ref_hard)},
        "n_generated": total,
        "n_correctness_passing": int(len(feas)),
        "n_mfe_competitive": int(len(comp)),
        "dominators_all": k_all, "dominators_correctness_passing": k_f,
        "rate_all_percent": 100 * k_all / total,
        "rate_passing_percent": 100 * k_f / max(len(feas), 1),
        "ci_passing_percent": [lo_f, hi_f],
        "supported": bool(supported),
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT_CSV.name} and {OUT_JSON.name}")


if __name__ == "__main__":
    main()
