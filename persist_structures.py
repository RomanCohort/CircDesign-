#!/usr/bin/env python3
"""Persist the secondary structure of all 261 candidates.

WHY THIS EXISTS, AND WHY IT DOES NOT JUST RE-RUN THE PIPELINE
------------------------------------------------------------
``combined_circdesign_analysis.py`` computes the MFE dot-bracket for every candidate
and then threw it away for all but the original, so any structural analysis of the
reference set needed a full refold. That line is fixed at the source.

Re-running the pipeline to regenerate the column would, however, OVERWRITE
``candidate_metrics.csv`` -- and that file, per the audit in the report, cannot be
reproduced by the delivered source (its mtime predates the script's, and the script
cannot emit the 1e5 sentinel the file contains). Overwriting it would silently
invalidate every downstream analysis. So this script computes the structures
SEPARATELY and writes a new file, leaving the original untouched.

WHAT IS STORED
--------------
  mfe_dot_bracket   the minimum-free-energy structure of the full 2013 nt circle,
                    ViennaRNA circular model, identical settings to the pipeline

The base-pair-probability matrices are NOT stored: a dense 2013x2013 float matrix is
~16 MB per candidate, i.e. ~4 GB for the set. The report documents that omission
rather than hiding it; the arrays are reproducible from the sequences in seconds.

Outputs structures_261.csv (candidate_id, mfe, mfe_dot_bracket).

Run:  python persist_structures.py [--workers 28]
"""

from __future__ import annotations

import argparse
import os
import time
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

import combined_circdesign_analysis as base

ROOT = Path(__file__).resolve().parent
MATURE = ROOT / "combined_mature_circRNA.fasta"
CDS = ROOT / "synonymous_candidate_CDS.fasta"
OUT = ROOT / "structures_261.csv"

CDS_MATURE = base.relative_interval(base.CDS_START, base.CDS_END)
_MATURE_RNA = ""


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
    if len(_MATURE_RNA) < 1000:
        raise RuntimeError("splice template unset -- workers must receive a finished "
                           "sequence from the parent, not rebuild it from a global")
    lo, hi = CDS_MATURE
    out = _MATURE_RNA[: lo - 1] + cds_dna.replace("T", "U") + _MATURE_RNA[hi:]
    if len(out) != len(_MATURE_RNA):
        raise RuntimeError("spliced length changed")
    return out


def worker(item):
    cid, seq = item                      # seq is already spliced, in the parent
    fc = base.RNA.fold_compound(seq, base.vienna_model(True), base.RNA.OPTION_MFE)
    dbn, mfe = fc.mfe()
    return {"candidate_id": cid, "mfe": float(mfe), "mfe_dot_bracket": dbn}


def main() -> None:
    global _MATURE_RNA
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 4))
    args = ap.parse_args()

    _MATURE_RNA = "".join(
        l.strip() for l in MATURE.read_text(encoding="utf-8").splitlines()
        if not l.startswith(">")
    ).upper().replace("T", "U")
    recs = read_fasta(CDS)

    metrics = pd.read_csv(ROOT / "candidate_metrics.csv", encoding="utf-8-sig")
    ids = metrics.candidate_id.tolist()
    missing = [i for i in ids if i not in recs]
    if missing:
        raise SystemExit(f"{len(missing)} candidate ids absent from the FASTA: {missing[:4]}")

    todo = [(cid, build_mature(recs[cid])) for cid in ids]
    if {len(s) for _, s in todo} != {len(_MATURE_RNA)}:
        raise SystemExit("guard failed: spliced lengths differ from the mature circle")

    probe = build_mature(recs["original"])
    fc0 = base.RNA.fold_compound(probe, base.vienna_model(True), base.RNA.OPTION_MFE)
    pm = float(fc0.mfe()[1])
    recorded = float(metrics.loc[metrics.candidate_id == "original",
                                 "circular_mfe_kcal_mol"].iloc[0])
    if abs(pm - recorded) > 1e-6:
        raise SystemExit(f"guard failed: construct folds to {pm}, recorded {recorded}")
    print(f"guard passed: construct reproduces MFE {pm:.4f}")

    print(f"folding {len(todo)} candidates on {args.workers} workers "
          f"(spliced in parent)...")
    t0 = time.time()
    with Pool(args.workers) as pool:
        rows = pool.map(worker, todo, chunksize=4)

    df = pd.DataFrame(rows)
    # sanity: recomputed MFE must match the delivered column
    chk = df.merge(metrics[["candidate_id", "circular_mfe_kcal_mol"]], on="candidate_id")
    dev = float((chk.mfe - chk.circular_mfe_kcal_mol).abs().max())
    print(f"guard: recomputed MFE vs delivered column, max |dev| = {dev:.6f}")
    if dev > 1e-6:
        raise SystemExit("guard failed: this is not the same reference set")
    print("guard passed: all 261 structures belong to the delivered reference set")

    if (df.mfe_dot_bracket.str.len() != len(_MATURE_RNA)).any():
        raise SystemExit("guard failed: a dot-bracket has the wrong length")

    df.to_csv(OUT, index=False, encoding="utf-8-sig")
    paired = df.mfe_dot_bracket.str.count(r"\(")
    print(f"\nwrote {OUT.name}: {len(df)} structures, "
          f"base pairs {paired.min()}-{paired.max()} (median {paired.median():.0f})")
    print(f"elapsed {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
