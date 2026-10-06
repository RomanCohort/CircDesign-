#!/usr/bin/env python3
"""Composition-axis / arrangement-axis decomposition of the synonymous reference space.

Purpose
-------
The retrospective report originally presented three null models (human-weighted,
uniform, composition-shuffle) as if they were interchangeable samples from one
"random synonymous" reference distribution.  They are not.  They differ in a
controlled way along two interpretable axes:

  * composition axis  -- WHICH synonymous codons are used (codon usage vector)
  * arrangement axis  -- WHERE those codons are placed (order along the CDS)

`composition_shuffle` holds the codon multiset exactly fixed and randomises only
the order, so it is the only null that isolates the arrangement axis.  The other
two co-vary composition and arrangement.

This script makes that decomposition explicit and quantitative:

  1. codon usage matrix U        (n x 61 sense codons)      -> composition space
  2. codon-pair residual matrix R (n x 3721 ordered pairs)  -> arrangement space
     where R = C - E and E is the exact expectation of the codon-pair count
     matrix under a uniformly random permutation of the SAME codon multiset
     (i.e. under the composition_shuffle null itself).
  3. PCA on each block; cross-validated variance decomposition of each metric
     into composition-explained / arrangement-incremental parts.

Outputs
-------
  axis_decomposition.json   machine-readable results (consumed by build_report.py)
  axis_scores.csv           per-candidate scores on the leading axes

Self-checks (fail loudly rather than silently produce a wrong decomposition):
  * every CDS must translate to the identical 360 aa protein
  * the 65 composition-shuffle + original candidates must share ONE codon usage
    vector exactly -- this is the algebraic fact that makes shuffle the
    arrangement-isolating null

Run:  python axis_decomposition.py
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import KFold, cross_val_score

ROOT = Path(__file__).resolve().parent
CANDIDATE_FASTA = ROOT / "synonymous_candidate_CDS.fasta"
METRICS_CSV = ROOT / "candidate_metrics.csv"
OUT_JSON = ROOT / "axis_decomposition.json"
OUT_SCORES = ROOT / "axis_scores.csv"

CODONS = [
    "TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG",
    "ATT", "ATC", "ATA", "ATG",
    "GTT", "GTC", "GTA", "GTG",
    "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG",
    "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
    "TAT", "TAC", "TAA", "TAG", "CAT", "CAC", "CAA", "CAG",
    "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
    "TGT", "TGC", "TGA", "TGG", "CGT", "CGC", "CGA", "CGG",
    "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG",
]
STOPS = {"TAA", "TAG", "TGA"}
SENSE = [c for c in CODONS if c not in STOPS]  # 61
SENSE_IDX = {c: i for i, c in enumerate(SENSE)}

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

METRIC_SPECS = [
    ("circular_mfe_kcal_mol", "Circular MFE", True),
    ("human_cai", "Human CAI", False),
    ("lires", "L- IRES deviation (666 nt)", True),
    ("lires_ires669", "L- IRES deviation (669 nt)", True),
]


def read_fasta(path: Path) -> "list[tuple[str, str]]":
    records: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = line[1:].split("|")[0]
            chunks = []
        else:
            chunks.append(line.upper())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def sense_codons(cds: str) -> "list[str]":
    """First 360 codons of the 1083 nt CDS (stop codon excluded)."""
    codons = [cds[i : i + 3] for i in range(0, len(cds) - 3, 3)]
    return codons


def translate(codons) -> str:
    return "".join(CODON_TABLE[c] for c in codons)


def pair_expectation(counts: np.ndarray) -> np.ndarray:
    """E[C_ij] for a uniformly random permutation of a multiset with counts `counts`.

    For a linear arrangement of N items with class counts n_i, the number of
    adjacent ordered pairs is N-1, and

        P(pos k = i, pos k+1 = j) = (n_i / N) * (n_j - [i==j]) / (N - 1)

    so E[C_ij] = n_i * (n_j - [i==j]) / N.
    """
    n = counts.astype(float)
    n_total = n.sum()
    expected = np.outer(n, n)
    expected[np.diag_indices_from(expected)] -= n  # n_i * (n_i - 1) / N
    return expected / n_total


def main() -> None:
    records = read_fasta(CANDIDATE_FASTA)
    metrics = pd.read_csv(METRICS_CSV, encoding="utf-8-sig").set_index("candidate_id")
    print(f"loaded {len(records)} CDS records, {len(metrics)} metric rows")

    ids: list[str] = []
    usage_rows: list[np.ndarray] = []
    resid_rows: list[np.ndarray] = []
    reference_protein = None

    for name, cds in records:
        codons = sense_codons(cds)
        if len(codons) != 360:
            raise SystemExit(f"{name}: expected 360 sense codons, got {len(codons)}")
        if any(c not in CODON_TABLE for c in codons):
            raise SystemExit(f"{name}: non-sense codon encountered")
        protein = translate(codons)
        if reference_protein is None:
            reference_protein = protein
        elif protein != reference_protein:
            raise SystemExit(f"{name}: protein differs from reference -- not synonymous")

        counts = np.zeros(len(SENSE))
        pair = np.zeros((len(SENSE), len(SENSE)))
        for k, c in enumerate(codons):
            counts[SENSE_IDX[c]] += 1
        for k in range(len(codons) - 1):
            pair[SENSE_IDX[codons[k]], SENSE_IDX[codons[k + 1]]] += 1

        ids.append(name)
        usage_rows.append(counts)
        resid_rows.append((pair - pair_expectation(counts)).ravel())

    U = np.vstack(usage_rows)          # n x 61
    R = np.vstack(resid_rows)          # n x 3721
    print(f"composition matrix {U.shape}, arrangement residual matrix {R.shape}")

    # --- self-check: shuffle + original must share ONE codon usage vector -------
    ens = metrics.loc[ids, "ensemble"].to_numpy()
    shuffle_mask = (ens == "composition_shuffle") | (ens == "current_construct")
    unique_usage = np.unique(U[shuffle_mask], axis=0)
    print(
        f"self-check  codon-usage uniqueness inside shuffle+original group: "
        f"{unique_usage.shape[0]} distinct vector(s) over {shuffle_mask.sum()} candidates"
    )
    if unique_usage.shape[0] != 1:
        raise SystemExit("shuffle group does NOT share one codon-usage vector -- assumption broken")
    # and confirm it differs from the other ensembles
    others = U[~shuffle_mask]
    min_dist = np.min(np.abs(others - unique_usage[0]).sum(axis=1))
    print(f"self-check  minimum L1 distance from that vector to any non-shuffle candidate: {min_dist:.1f}")

    # --- PCA on each block ------------------------------------------------------
    n_comp_pcs = 15
    n_arr_pcs = 30
    # centre within the pooled random ensembles only; used for regression below
    pca_comp = PCA(n_components=n_comp_pcs, random_state=0).fit(U)
    pca_arr = PCA(n_components=n_arr_pcs, random_state=0).fit(R)
    Zc = pca_comp.transform(U)
    Za = pca_arr.transform(R)
    ev_c = pca_comp.explained_variance_ratio_
    ev_a = pca_arr.explained_variance_ratio_
    print(f"composition block: PC1..{n_comp_pcs} explain {ev_c.sum()*100:.1f}% (PC1 {ev_c[0]*100:.1f}%)")
    print(f"arrangement block: PC1..{n_arr_pcs} explain {ev_a.sum()*100:.1f}% (PC1 {ev_a[0]*100:.1f}%)")

    # --- cross-validated variance decomposition per metric ----------------------
    # Two bounding controls (cai_max, gc3_high) carry a FAILED partition function
    # recorded as the sentinel ensemble_free_energy_kcal_mol == 1e5 with an empty
    # pair list, so their lires equals ||p_ref||_2 by construction. They are also the
    # extreme corners of the composition space, i.e. exactly the points with the most
    # leverage on a linear fit. Including them moves the lires regression across zero:
    # measured R^2 -0.1135 (contaminated) vs +0.0222 (clean). The sign is not a finding,
    # it is an artefact of two broken rows, so the clean fit is reported as primary and
    # the contaminated one is kept only so the difference is visible.
    sentinel = (
        metrics.loc[ids, "ensemble_free_energy_kcal_mol"].to_numpy(dtype=float) >= 1e4
    )
    if sentinel.any():
        bad = [i for i, s in zip(ids, sentinel) if s]
        print(f"\nEXCLUDED from the variance decomposition: {len(bad)} row(s) carrying the "
              f"partition-function sentinel -- {bad}")
        print("  (their lires is an artefact: empty pair list -> lires == ||p_ref||_2)")

    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    alphas = np.logspace(-3, 3, 25)
    keep = ~sentinel
    results = {}
    print()
    print(f"{'metric':<28}{'R2 comp':>10}{'R2 comp+arr':>13}{'arr incr':>10}"
          f"{'[contam]':>11}")
    for col, label, _ in METRIC_SPECS:
        y_all = metrics.loc[ids, col].to_numpy(dtype=float)
        y = y_all[keep]
        r2_c = cross_val_score(
            RidgeCV(alphas=alphas), Zc[keep], y, cv=cv, scoring="r2"
        ).mean()
        r2_ca = cross_val_score(
            RidgeCV(alphas=alphas), np.hstack([Zc, Za])[keep], y, cv=cv, scoring="r2"
        ).mean()
        r2_c_contam = cross_val_score(
            RidgeCV(alphas=alphas), Zc, y_all, cv=cv, scoring="r2"
        ).mean()
        results[col] = {
            "label": label,
            "r2_composition_only": float(r2_c),
            "r2_composition_plus_arrangement": float(r2_ca),
            "r2_arrangement_incremental": float(r2_ca - r2_c),
            "r2_composition_only_contaminated": float(r2_c_contam),
            "n_rows_used": int(keep.sum()),
            "n_rows_excluded_sentinel": int((~keep).sum()),
            "sign_flips_when_cleaned": bool((r2_c < 0) != (r2_c_contam < 0)),
        }
        flip = " SIGN FLIP" if results[col]["sign_flips_when_cleaned"] else ""
        print(f"{label:<28}{r2_c:>10.3f}{r2_ca:>13.3f}{r2_ca - r2_c:>10.3f}"
              f"{r2_c_contam:>11.3f}{flip}")

    # --- metric-space effective rank (restricted to the random pool) ------------
    pool = metrics.loc[ids, "ensemble"].isin(
        ["human_weighted", "uniform", "composition_shuffle"]
    ).to_numpy()
    M = np.column_stack(
        [metrics.loc[ids, c].to_numpy(dtype=float) for c, _, _ in METRIC_SPECS[:3]]
    )
    Ms = (M[pool] - M[pool].mean(0)) / M[pool].std(0, ddof=1)
    corr = np.corrcoef(Ms.T)
    eig = np.linalg.eigvalsh(corr)[::-1]
    eff_rank = float(eig.sum() ** 2 / (eig**2).sum())
    print()
    print(f"metric-space correlation (MFE, CAI, LIRES666):\n{np.round(corr, 3)}")
    print(f"metric-space eigenvalues {np.round(eig, 4)}  effective rank {eff_rank:.3f} / 3")

    # --- how much arrangement variance each null actually carries ---------------
    var_by_ensemble = {}
    for ens_name in ["human_weighted", "uniform", "composition_shuffle"]:
        m = (ens == ens_name)
        var_by_ensemble[ens_name] = {
            "composition_pc1_sd": float(Zc[m, 0].std(ddof=1)),
            "arrangement_pc1_sd": float(Za[m, 0].std(ddof=1)),
        }
    print("\nspread along the two axes, per ensemble (SD of leading PC score):")
    for k, v in var_by_ensemble.items():
        print(f"  {k:<22} composition {v['composition_pc1_sd']:9.3f}   arrangement {v['arrangement_pc1_sd']:9.3f}")

    payload = {
        "description": "composition / arrangement decomposition of the synonymous reference space",
        "definitions": {
            "composition_axis": "61-dim sense-codon usage vector of the 1083 nt CDS",
            "arrangement_axis": (
                "3721-dim codon-pair count matrix minus its exact expectation under a "
                "uniformly random permutation of the SAME codon multiset, i.e. under the "
                "composition_shuffle null"
            ),
            "n_composition_pcs": n_comp_pcs,
            "n_arrangement_pcs": n_arr_pcs,
            "n_sequences": len(ids),
        },
        "self_checks": {
            "protein_identical_across_all_candidates": True,
            "shuffle_and_original_share_one_codon_usage_vector": unique_usage.shape[0] == 1,
            "min_l1_distance_usage_to_non_shuffle": float(min_dist),
        },
        "explained_variance": {
            "composition_pcs_cumulative": [float(x) for x in np.cumsum(ev_c)],
            "arrangement_pcs_cumulative": [float(x) for x in np.cumsum(ev_a)],
        },
        "variance_decomposition": results,
        "metric_space": {
            "correlation_matrix": [[float(x) for x in row] for row in corr],
            "eigenvalues": [float(x) for x in eig],
            "effective_rank": eff_rank,
            "labels": ["circular_mfe", "human_cai", "lires_666"],
        },
        "axis_spread_by_ensemble": var_by_ensemble,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")

    scores = pd.DataFrame({"candidate_id": ids, "ensemble": ens})
    for i in range(4):
        scores[f"composition_pc{i + 1}"] = Zc[:, i]
    for i in range(4):
        scores[f"arrangement_pc{i + 1}"] = Za[:, i]
    scores.to_csv(OUT_SCORES, index=False, encoding="utf-8-sig")
    print(f"\nwrote {OUT_JSON.name} and {OUT_SCORES.name}")


if __name__ == "__main__":
    main()
