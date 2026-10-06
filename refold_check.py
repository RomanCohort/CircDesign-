#!/usr/bin/env python3
"""Decisive check: does my fold reproduce the delivered MFE for the ORIGINAL 64
composition_shuffle sequences?

The pre-registered advantage-space run found a domination rate an order of magnitude
below the delivered one (26/2000 = 1.3% vs 5/64 = 7.8%). Under a binomial with
p = 0.013 the delivered 5/64 has probability ~4e-4, so the two samples cannot come
from the same distribution. Either my fold differs from the delivered fold, or my
shuffle sampler draws from a different distribution than the delivered one. This
script separates the two by refolding the DELIVERED shuffle sequences with MY code
path and comparing, sequence by sequence, against the delivered MFE column.

Writes refold_check.json.
"""

from __future__ import annotations

import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

import advantage_space as A
import combined_circdesign_analysis as base

ROOT = Path(__file__).resolve().parent
CON_MFE = -779.4000244140625


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


def fold(item):
    cid, seq = item
    fc = base.RNA.fold_compound(seq, base.vienna_model(True), base.RNA.OPTION_MFE)
    return cid, float(fc.mfe()[1]), len(seq)


def main() -> None:
    A._MATURE_RNA = "".join(
        l.strip() for l in (ROOT / "combined_mature_circRNA.fasta").read_text(
            encoding="utf-8").splitlines() if not l.startswith(">")
    ).upper().replace("T", "U")

    recs = read_fasta(ROOT / "synonymous_candidate_CDS.fasta")
    m = pd.read_csv(ROOT / "candidate_metrics.csv", encoding="utf-8-sig")
    sh = m[m.ensemble == "composition_shuffle"].candidate_id.tolist()
    print(f"delivered composition_shuffle sequences: {len(sh)}")

    items = [(cid, A.build_mature(recs[cid])) for cid in sh]
    lens = {len(s) for _, s in items}
    print(f"spliced length set: {lens} (expect {{2013}})")
    if lens != {2013}:
        raise SystemExit("splice produced the wrong molecule -- aborting")

    with Pool(28) as pool:
        res = pool.map(fold, items, chunksize=2)

    rec = m.set_index("candidate_id")["circular_mfe_kcal_mol"]
    mf = np.array([r[1] for r in res])
    dev = np.array([abs(r[1] - rec[r[0]]) for r in res])

    print(f"\nrefolded MFE: min={mf.min():.2f} median={np.median(mf):.2f} max={mf.max():.2f}")
    print(f"delivered   : min={rec[sh].min():.2f} median={rec[sh].median():.2f} "
          f"max={rec[sh].max():.2f}")
    print(f"max per-sequence |deviation| = {dev.max():.6f} kcal/mol")
    print(f"bit-identical on all 64: {bool(dev.max() < 1e-9)}")

    mine = int((mf <= CON_MFE).sum())
    print(f"\ncompetitive (mfe <= construct) with MY fold      : {mine}/64")
    print(f"competitive according to the delivered CSV       : "
          f"{int((rec[sh] <= CON_MFE).sum())}/64")

    verdict = ("fold path IDENTICAL -- the discrepancy is in SAMPLE GENERATION, "
               "not in the fold" if dev.max() < 1e-9 else
               "FOLD PATH DIFFERS -- every downstream MFE-derived number is suspect")
    print(f"\nVERDICT: {verdict}")

    (ROOT / "refold_check.json").write_text(json.dumps({
        "n_delivered_shuffles": len(sh),
        "max_abs_deviation_kcal_mol": float(dev.max()),
        "fold_identical": bool(dev.max() < 1e-9),
        "competitive_with_my_fold": mine,
        "competitive_delivered": int((rec[sh] <= CON_MFE).sum()),
        "refolded_mfe_min": float(mf.min()), "refolded_mfe_median": float(np.median(mf)),
        "refolded_mfe_max": float(mf.max()),
        "verdict": verdict,
    }, indent=2), encoding="utf-8")
    print(f"\nwrote refold_check.json")


if __name__ == "__main__":
    main()
