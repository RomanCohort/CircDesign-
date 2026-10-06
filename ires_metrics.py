#!/usr/bin/env python3
"""Replace the L_IRES proxy with three directly interpretable IRES metrics.

Why
---
The delivered report measures IRES integrity as ``L_IRES``: the L2 distance
between the full-circle base-pair-probability matrix and an IRES-only constrained
reference matrix.  A composition-matched null analysis showed this quantity is not
predictable from any codon-level feature: cross-validated R^2 is about 0.02 for
codon usage and the arrangement block adds no measurable increment (see
axis_decomposition.json, whose fit EXCLUDES the two rows carrying the
partition-function sentinel -- with them included the value is -0.11, which is an
artefact of two broken rows, not a finding).  A quantity the sequence does not
determine cannot serve as an optimisation objective.  It also has no verbal
reading: "L_IRES = 11.9540" does not tell you what happened to the IRES.

This script keeps the *same physics* -- identical ViennaRNA model, identical
constrained reference, identical BPP cutoff, imported from
``combined_circdesign_analysis`` so the two are directly comparable -- and reports
three quantities that each answer a question in words:

  ires_native_recall_hard   of the base pairs the IRES forms on its own
                            (reference p >= 0.5), what fraction still form
                            (candidate p >= 0.5) in the full circle?
                            -> "the IRES keeps X% of its own structure"

  ires_native_recall_soft   the same, but probability-weighted instead of
                            thresholded: mean candidate probability over the
                            reference pair set.  Smooth, no threshold artefact.

  ires_cross_occupancy      of the IRES's own residues, what fraction acquire a
                            partner OUTSIDE the IRES (p >= 0.5)?
                            -> "how much of the IRES has the rest of the circle
                               taken over?"

  ires_unpaired_frac        of the IRES's own residues, what fraction are
                            effectively unpaired in the full-circle ensemble
                            (total pairing probability < 0.5)?
                            -> accessibility; an IRES must be readable

All four are ratios, so they have a null model and a CI, unlike L2 norms.

Outputs ``ires_metrics.csv`` (per candidate, both IRES windows) and
``ires_metrics.json`` (reference construction + definitions).

Run:  python ires_metrics.py [--workers 8]
"""

from __future__ import annotations

import argparse
import json
import os
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

import combined_circdesign_analysis as base

ROOT = Path(__file__).resolve().parent
MATURE_FASTA = ROOT / "combined_mature_circRNA.fasta"
CDS_FASTA = ROOT / "synonymous_candidate_CDS.fasta"
OUT_CSV = ROOT / "ires_metrics.csv"
OUT_JSON = ROOT / "ires_metrics.json"

CUTOFF = base.DEFAULT_BPP_CUTOFF
PAIR_THRESHOLD = 0.5

# mature-circle coordinates (1-based inclusive) -- from summary.json
WINDOWS = {
    "666": base.relative_interval(base.IRES_START, base.IRES_END_PRIMARY),
    "669": base.relative_interval(base.IRES_START, base.IRES_END_SENSITIVITY),
}
CDS_MATURE = base.relative_interval(base.CDS_START, base.CDS_END)

_REFS: dict = {}
_WINDOWS: dict = {}


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


def build_mature(mature_dna: str, cds_dna: str) -> str:
    """Splice a candidate CDS (DNA) into the fixed mature circle (DNA).

    Everything stays in DNA here and is converted to RNA only at fold time, matching
    combined_circdesign_analysis.evaluate_candidate (`full_rna = mature_dna.replace("T","U")`).
    Converting on one side only silently produces a sequence that is neither.
    """
    lo, hi = CDS_MATURE
    lo0, hi0 = lo - 1, hi  # 1-based inclusive -> 0-based half-open
    return mature_dna[:lo0] + cds_dna + mature_dna[hi0:]


def _init(refs, windows):
    global _REFS, _WINDOWS
    _REFS, _WINDOWS = refs, windows


def _metrics_for_window(cand_plist, ref, interval):
    start, end = interval
    span = end - start + 1

    # reference pairs the IRES forms on its own
    n_ref_hard = {p: v for p, v in ref.items() if v >= PAIR_THRESHOLD}
    if n_ref_hard:
        still = sum(1 for p in n_ref_hard if cand_plist.get(p, 0.0) >= PAIR_THRESHOLD)
        recall_hard = still / len(n_ref_hard)
        recall_soft = sum(cand_plist.get(p, 0.0) for p in n_ref_hard) / len(n_ref_hard)
    else:
        recall_hard = recall_soft = float("nan")

    # per-residue total pairing probability inside the IRES, from the candidate
    total = {i: 0.0 for i in range(start, end + 1)}
    cross_res = set()
    for (i, j), p in cand_plist.items():
        i_in = start <= i <= end
        j_in = start <= j <= end
        if i_in and j_in:
            total[i] += p
            total[j] += p
        elif i_in:
            total[i] += p
            if p >= PAIR_THRESHOLD:
                cross_res.add(i)
        elif j_in:
            total[j] += p
            if p >= PAIR_THRESHOLD:
                cross_res.add(j)

    unpaired = sum(1 for v in total.values() if v < PAIR_THRESHOLD)
    return {
        "ires_native_recall_hard": recall_hard,
        "ires_native_recall_soft": recall_soft,
        "ires_n_ref_pairs": len(n_ref_hard),
        "ires_cross_occupancy": len(cross_res) / span,
        "ires_unpaired_frac": unpaired / span,
        "ires_mean_pair_prob": sum(total.values()) / span,
    }


def _worker(rec):
    cid, mature_dna = rec
    fc = base.RNA.fold_compound(mature_dna.replace("T", "U"), base.vienna_model(True),
                                base.RNA.OPTION_MFE | base.RNA.OPTION_PF)
    structure, mfe = fc.mfe()
    fc.pf()
    plist = {(int(it.i), int(it.j)): float(it.p) for it in fc.plist_from_probs(CUTOFF)}

    row = {"candidate_id": cid, "circular_mfe_kcal_mol": float(mfe)}
    for tag, interval in _WINDOWS.items():
        m = _metrics_for_window(plist, _REFS[tag], interval)
        row.update({f"{k}_{tag}": v for k, v in m.items()})
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    mature_dna = "".join(
        l.strip() for l in MATURE_FASTA.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("U", "T")
    records = read_fasta(CDS_FASTA)
    print(f"mature circle {len(mature_dna)} nt (DNA); {len(records)} candidate CDS")

    # sanity: splicing the ORIGINAL cds back must reproduce the mature circle exactly
    orig_cds = next(c for n, c in records if n == "original")
    if build_mature(mature_dna, orig_cds) != mature_dna:
        raise SystemExit("splice check failed: original CDS does not round-trip")
    print("splice check: original CDS round-trips to the mature circle exactly")

    print("building IRES-only constrained references (same construction as L_IRES)...")
    mature_rna = mature_dna.replace("T", "U")
    refs = {}
    for tag, interval in WINDOWS.items():
        ref, info = base.constrained_ires_reference(mature_rna, interval, CUTOFF)
        refs[tag] = ref
        n_hard = sum(1 for v in ref.values() if v >= PAIR_THRESHOLD)
        print(f"  window {tag}: interval={interval}  ref pairs={info['retained_pairs']}  "
              f"p>=0.5: {n_hard}  dG={info['ensemble_free_energy_kcal_mol']:.2f}  "
              f"({info['elapsed_seconds']:.1f}s)")

    seqs = [(n, build_mature(mature_dna, c)) for n, c in records]
    print(f"folding {len(seqs)} candidates on {args.workers} workers...")
    t0 = time.time()
    with Pool(args.workers, initializer=_init, initargs=(refs, WINDOWS)) as pool:
        rows = []
        for k, row in enumerate(pool.imap_unordered(_worker, seqs), 1):
            rows.append(row)
            if k % 20 == 0 or k == len(seqs):
                el = time.time() - t0
                print(f"  {k}/{len(seqs)}  {el:.0f}s  eta {el / k * (len(seqs) - k):.0f}s",
                      flush=True)

    df = pd.DataFrame(rows).sort_values("candidate_id").reset_index(drop=True)
    df.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")

    payload = {
        "what": "interpretable IRES metrics, same physics as L_IRES",
        "replaces": "L_IRES (L2 distance between BPP matrices)",
        "reference_construction": (
            "circular full-length partition function with every non-IRES nucleotide "
            "forced unpaired; identical to constrained_ires_reference() used for L_IRES"
        ),
        "bpp_cutoff": CUTOFF,
        "pair_threshold": PAIR_THRESHOLD,
        "windows": {k: list(v) for k, v in WINDOWS.items()},
        "definitions": {
            "ires_native_recall_hard": "of ref pairs with p>=0.5, fraction still p>=0.5 in candidate",
            "ires_native_recall_soft": "mean candidate probability over the same ref pair set",
            "ires_cross_occupancy": "fraction of IRES residues that pair outside the IRES (p>=0.5)",
            "ires_unpaired_frac": "fraction of IRES residues with total pairing probability <0.5",
            "ires_mean_pair_prob": "mean total pairing probability per IRES residue",
        },
        "n_candidates": int(len(df)),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT_CSV.name} and {OUT_JSON.name}  ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
