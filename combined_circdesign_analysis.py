#!/usr/bin/env python3
"""Retrospective circDesign-inspired audit of the JLU-FBH Combined circRNA.

This program deliberately does not claim to reproduce the unpublished circDesign
production code.  It implements the three published post-hoc objectives:

1. circular minimum free energy (MFE),
2. human codon adaptation index (CAI), and
3. IRES structural deviation (L2 distance between base-pairing-probability
   matrices for the unconstrained full circle and an IRES-only constrained
   full circle).

The exact mature circle and CDS coordinates are reconstructed from the uploaded
SnapGene file.  Synonymous benchmarking uses deterministic, reproducible random
ensembles that all encode the identical 360-aa antigen.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import shutil
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import RNA
from Bio import SeqIO
from Bio.Seq import Seq


# Uploaded construct coordinates (1-based, inclusive) determined from SnapGene
# annotations and verified by ORF continuity.
MATURE_START = 723
MATURE_END = 2735
CDS_START = 1504
CDS_END = 2586
IRES_START = 823
IRES_END_PRIMARY = 1488  # explicit synIRES-RC25 feature
IRES_END_SENSITIVITY = 1491  # broader HRV-B3 IRES annotation

DEFAULT_SEED = 20260804
DEFAULT_BPP_CUTOFF = 1e-5


# Kazusa Homo sapiens codon-usage table (species 9606; 93,487 CDSs,
# 40,662,582 codons).  Counts are used rather than rounded fractions.
_HUMAN_USAGE_TEXT = """
TTT F 714298 TCT S 618711 TAT Y 495699 TGT C 430311
TTC F 824692 TCC S 718892 TAC Y 622407 TGC C 513028
TTA L 311881 TCA S 496448 TAA * 40285 TGA * 63237
TTG L 525688 TCG S 179419 TAG * 32109 TGG W 535595
CTT L 536515 CCT P 713233 CAT H 441711 CGT R 184609
CTC L 796638 CCC P 804620 CAC H 613713 CGC R 423516
CTA L 290751 CCA P 688038 CAA Q 501911 CGA R 250760
CTG L 1611801 CCG P 281570 CAG Q 1391973 CGG R 464485
ATT I 650473 ACT T 533609 AAT N 689701 AGT S 493429
ATC I 846466 ACC T 768147 AAC N 776603 AGC S 791383
ATA I 304565 ACA T 614523 AAA K 993621 AGA R 494682
ATG M 896005 ACG T 246105 AAG K 1295568 AGG R 486463
GTT V 448607 GCT A 750096 GAT D 885429 GGT G 437126
GTC V 588138 GCC A 1127679 GAC D 1020595 GGC G 903565
GTA V 287712 GCA A 643471 GAA E 1177632 GGA G 669873
GTG V 1143534 GCG A 299495 GAG E 1609975 GGG G 669768
"""


def parse_human_usage() -> Tuple[Dict[str, Tuple[str, int]], Dict[str, List[Tuple[str, int]]]]:
    toks = _HUMAN_USAGE_TEXT.split()
    codon_usage: Dict[str, Tuple[str, int]] = {}
    by_aa: Dict[str, List[Tuple[str, int]]] = defaultdict(list)
    for i in range(0, len(toks), 3):
        codon, aa, count = toks[i], toks[i + 1], int(toks[i + 2])
        codon_usage[codon] = (aa, count)
        by_aa[aa].append((codon, count))
    return codon_usage, dict(by_aa)


CODON_USAGE, CODONS_BY_AA = parse_human_usage()
CODON_WEIGHTS = {
    codon: count / max(n for _, n in CODONS_BY_AA[aa])
    for codon, (aa, count) in CODON_USAGE.items()
}


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    ensemble: str
    cds_dna: str
    mature_dna: str
    generation_index: int


def wrap_fasta(sequence: str, width: int = 80) -> str:
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def md5(sequence: str) -> str:
    return hashlib.md5(sequence.encode("ascii")).hexdigest()


def human_cai(cds_dna: str) -> float:
    """Sharp-Li style CAI from Kazusa human relative synonymous usage.

    The terminal stop is excluded.  Single-codon amino acids have weight 1 and
    therefore do not alter the geometric mean.
    """

    if len(cds_dna) % 3:
        raise ValueError("CDS length is not divisible by three")
    codons = [cds_dna[i : i + 3] for i in range(0, len(cds_dna) - 3, 3)]
    if not codons:
        raise ValueError("CDS has no sense codons")
    weights = []
    for codon in codons:
        try:
            weights.append(CODON_WEIGHTS[codon])
        except KeyError as exc:
            raise ValueError(f"Unknown codon {codon}") from exc
    return math.exp(sum(math.log(w) for w in weights) / len(weights))


def gc_fraction(sequence: str) -> float:
    return sum(base in "GC" for base in sequence) / len(sequence)


def gc3_fraction(cds_dna: str) -> float:
    sense = [cds_dna[i : i + 3] for i in range(0, len(cds_dna) - 3, 3)]
    return sum(codon[2] in "GC" for codon in sense) / len(sense)


def translate_cds(cds_dna: str) -> str:
    return str(Seq(cds_dna).translate())


def relative_interval(global_start: int, global_end: int) -> Tuple[int, int]:
    """Return 1-based inclusive positions within the reconstructed mature circle."""

    return global_start - MATURE_START + 1, global_end - MATURE_START + 1


def extract_construct(snapgene_path: Path) -> Dict[str, object]:
    record = SeqIO.read(str(snapgene_path), "snapgene")
    plasmid = str(record.seq).upper()
    mature = plasmid[MATURE_START - 1 : MATURE_END]
    cds = plasmid[CDS_START - 1 : CDS_END]
    ires_primary = plasmid[IRES_START - 1 : IRES_END_PRIMARY]
    ires_sensitivity = plasmid[IRES_START - 1 : IRES_END_SENSITIVITY]
    protein = translate_cds(cds)

    if len(plasmid) != 5267:
        raise ValueError(f"Unexpected plasmid length: {len(plasmid)} bp")
    if len(mature) != 2013:
        raise ValueError(f"Unexpected mature circle length: {len(mature)} nt")
    if len(cds) != 1083:
        raise ValueError(f"Unexpected CDS length: {len(cds)} nt")
    if not cds.startswith("ATG") or cds[-3:] not in {"TAA", "TAG", "TGA"}:
        raise ValueError("Expected canonical start and stop codons")
    if protein.count("*") != 1 or not protein.endswith("*"):
        raise ValueError("CDS contains an internal stop or lacks the terminal stop")
    cds_rel_start, cds_rel_end = relative_interval(CDS_START, CDS_END)
    if mature[cds_rel_start - 1 : cds_rel_end] != cds:
        raise AssertionError("CDS does not map into mature circle")

    feature_rows = []
    wanted = {
        "T4td US intron retained",
        "Intron Scar 2",
        "PABPspacer50v3",
        "synIRES-RC25",
        "Human rhinovirus 3 5' UTR IRES",
        "Kozak sequence",
        "tPA signal/pro sequence",
        "TROP2 CTL1",
        "TROP2 CTL2",
        "B7H4 CTL1",
        "B7H4 CTL2",
        "Nectin4 CTL1",
        "Nectin4 CTL2",
        "LIV-1 CTL2",
        "LIV-1 CTL3",
        "mC3dP28-1",
        "mC3dP28-2",
        "mC3dP28-3",
        "mC3dP28-4",
        "MITD",
        "FLAG",
        "HBA1 full 3' UTR",
        "T4td DS intron retained",
        "Intron Scar 1",
    }
    seen = set()
    for feature in record.features:
        labels = feature.qualifiers.get("label", [])
        if not labels:
            continue
        label = str(labels[0])
        if label not in wanted:
            continue
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        key = (label, start, end)
        if key in seen:
            continue
        seen.add(key)
        feature_rows.append(
            {
                "feature": label,
                "plasmid_start": start,
                "plasmid_end": end,
                "mature_start": start - MATURE_START + 1 if MATURE_START <= start <= MATURE_END else None,
                "mature_end": end - MATURE_START + 1 if MATURE_START <= end <= MATURE_END else None,
                "length_nt": end - start + 1,
                "type": feature.type,
            }
        )

    return {
        "record_id": record.id,
        "record_name": record.name,
        "plasmid": plasmid,
        "mature": mature,
        "cds": cds,
        "protein": protein[:-1],
        "protein_with_stop": protein,
        "ires_primary": ires_primary,
        "ires_sensitivity": ires_sensitivity,
        "features": feature_rows,
        "junction": mature[-32:] + "|" + mature[:32],
    }


def synonymous_cds_from_choices(protein: str, chooser) -> str:
    codons = [chooser(aa, index) for index, aa in enumerate(protein)]
    return "".join(codons) + "TGA"


def weighted_choice(rng: random.Random, options: Sequence[Tuple[str, int]]) -> str:
    total = sum(weight for _, weight in options)
    draw = rng.randrange(total)
    running = 0
    for codon, weight in options:
        running += weight
        if draw < running:
            return codon
    return options[-1][0]


def replace_cds(mature: str, cds_dna: str) -> str:
    cds_rel_start, cds_rel_end = relative_interval(CDS_START, CDS_END)
    out = mature[: cds_rel_start - 1] + cds_dna + mature[cds_rel_end:]
    if len(out) != len(mature):
        raise AssertionError("Synonymous replacement changed mature circle length")
    return out


def build_candidates(
    mature: str,
    original_cds: str,
    protein: str,
    n_human: int,
    n_uniform: int,
    n_shuffle: int,
    seed: int,
) -> List[Candidate]:
    candidates: List[Candidate] = []
    seen = {original_cds}

    def add(candidate_id: str, ensemble: str, cds: str, generation_index: int) -> None:
        if translate_cds(cds) != protein + "*":
            raise AssertionError(f"{candidate_id} does not preserve the protein")
        candidates.append(
            Candidate(candidate_id, ensemble, cds, replace_cds(mature, cds), generation_index)
        )

    add("original", "current_construct", original_cds, 0)

    def generate_unique(ensemble: str, count: int, start_index: int, generator) -> None:
        made = 0
        attempts = 0
        while made < count:
            attempts += 1
            cds = generator(made, attempts)
            if cds in seen:
                continue
            seen.add(cds)
            made += 1
            add(f"{ensemble}_{made:03d}", ensemble, cds, start_index + made)

    human_rng = random.Random(seed + 101)

    def human_generator(_made: int, _attempt: int) -> str:
        return synonymous_cds_from_choices(
            protein,
            lambda aa, _idx: weighted_choice(human_rng, CODONS_BY_AA[aa]),
        )

    generate_unique("human_weighted", n_human, 1, human_generator)

    uniform_rng = random.Random(seed + 202)

    def uniform_generator(_made: int, _attempt: int) -> str:
        return synonymous_cds_from_choices(
            protein,
            lambda aa, _idx: uniform_rng.choice(CODONS_BY_AA[aa])[0],
        )

    generate_unique("uniform", n_uniform, 1 + n_human, uniform_generator)

    original_sense_codons = [original_cds[i : i + 3] for i in range(0, len(original_cds) - 3, 3)]
    aa_positions: Dict[str, List[int]] = defaultdict(list)
    for idx, aa in enumerate(protein):
        aa_positions[aa].append(idx)
    shuffle_rng = random.Random(seed + 303)

    def shuffle_generator(_made: int, _attempt: int) -> str:
        codons = list(original_sense_codons)
        for positions in aa_positions.values():
            values = [codons[pos] for pos in positions]
            shuffle_rng.shuffle(values)
            for pos, value in zip(positions, values):
                codons[pos] = value
        return "".join(codons) + original_cds[-3:]

    generate_unique("composition_shuffle", n_shuffle, 1 + n_human + n_uniform, shuffle_generator)

    def deterministic(strategy: str) -> str:
        def choose(aa: str, _idx: int) -> str:
            options = CODONS_BY_AA[aa]
            if strategy == "cai_max":
                return max(options, key=lambda x: (x[1], x[0]))[0]
            if strategy == "cai_min":
                return min(options, key=lambda x: (x[1], x[0]))[0]
            if strategy == "gc3_high":
                return max(options, key=lambda x: (x[0][2] in "GC", x[1], x[0]))[0]
            if strategy == "gc3_low":
                return max(options, key=lambda x: (x[0][2] in "AT", x[1], x[0]))[0]
            raise ValueError(strategy)

        return synonymous_cds_from_choices(protein, choose)

    for offset, strategy in enumerate(("cai_max", "cai_min", "gc3_high", "gc3_low"), 1):
        cds = deterministic(strategy)
        if cds in seen:
            continue
        seen.add(cds)
        add(strategy, "bounding_control", cds, 10_000 + offset)

    return candidates


def vienna_model(circular: bool) -> RNA.md:
    md = RNA.md()
    md.circ = 1 if circular else 0
    md.temperature = 37.0
    md.dangles = 2
    md.noLP = 0
    return md


def constrained_ires_reference(
    full_rna: str,
    ires_interval: Tuple[int, int],
    cutoff: float,
) -> Tuple[Dict[Tuple[int, int], float], Dict[str, float]]:
    """Full-circle partition function with all non-IRES nucleotides forced unpaired."""

    start, end = ires_interval
    fc = RNA.fold_compound(full_rna, vienna_model(True), RNA.OPTION_PF)
    for pos in range(1, len(full_rna) + 1):
        if not (start <= pos <= end):
            fc.hc_add_up(pos, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    t0 = time.time()
    ensemble_free_energy = float(fc.pf()[1])
    plist = fc.plist_from_probs(cutoff)
    ref = {
        (int(item.i), int(item.j)): float(item.p)
        for item in plist
        if start <= item.i <= end and start <= item.j <= end
    }
    return ref, {
        "ensemble_free_energy_kcal_mol": ensemble_free_energy,
        "retained_pairs": len(ref),
        "elapsed_seconds": time.time() - t0,
    }


_WORKER_REF_PRIMARY: Dict[Tuple[int, int], float] = {}
_WORKER_REF_SENSITIVITY: Dict[Tuple[int, int], float] = {}
_WORKER_IRES_PRIMARY: Tuple[int, int] = (0, 0)
_WORKER_IRES_SENSITIVITY: Tuple[int, int] = (0, 0)
_WORKER_CUTOFF = DEFAULT_BPP_CUTOFF


def init_worker(
    ref_primary: Dict[Tuple[int, int], float],
    ref_sensitivity: Dict[Tuple[int, int], float],
    ires_primary: Tuple[int, int],
    ires_sensitivity: Tuple[int, int],
    cutoff: float,
) -> None:
    global _WORKER_REF_PRIMARY, _WORKER_REF_SENSITIVITY
    global _WORKER_IRES_PRIMARY, _WORKER_IRES_SENSITIVITY, _WORKER_CUTOFF
    _WORKER_REF_PRIMARY = ref_primary
    _WORKER_REF_SENSITIVITY = ref_sensitivity
    _WORKER_IRES_PRIMARY = ires_primary
    _WORKER_IRES_SENSITIVITY = ires_sensitivity
    _WORKER_CUTOFF = cutoff


def ires_deviation(
    plist,
    reference: Mapping[Tuple[int, int], float],
    interval: Tuple[int, int],
) -> Dict[str, float]:
    start, end = interval
    candidate: Dict[Tuple[int, int], float] = {
        (int(item.i), int(item.j)): float(item.p)
        for item in plist
        if start <= item.i <= end or start <= item.j <= end
    }
    internal_sq = 0.0
    cross_sq = 0.0
    cross_mass = 0.0
    for pair in set(candidate).union(reference):
        p_cand = candidate.get(pair, 0.0)
        p_ref = reference.get(pair, 0.0)
        delta_sq = (p_cand - p_ref) ** 2
        if start <= pair[0] <= end and start <= pair[1] <= end:
            internal_sq += delta_sq
        else:
            cross_sq += delta_sq
            cross_mass += p_cand
    # ViennaRNA returns the unique upper triangle.  Multiplication by two
    # reports the L2 norm of the full symmetric probability matrix in Eq. 5.
    return {
        "lires": math.sqrt(2.0 * (internal_sq + cross_sq)),
        "lires_internal": math.sqrt(2.0 * internal_sq),
        "lires_crosstalk": math.sqrt(2.0 * cross_sq),
        "ires_cross_pair_probability_mass": cross_mass,
        "ires_candidate_pairs_retained": len(candidate),
    }


def evaluate_candidate(candidate: Candidate) -> Dict[str, object]:
    t0 = time.time()
    full_rna = candidate.mature_dna.replace("T", "U")
    fc = RNA.fold_compound(full_rna, vienna_model(True), RNA.OPTION_MFE | RNA.OPTION_PF)
    structure, mfe = fc.mfe()
    ensemble_free_energy = float(fc.pf()[1])
    plist = fc.plist_from_probs(_WORKER_CUTOFF)
    primary = ires_deviation(plist, _WORKER_REF_PRIMARY, _WORKER_IRES_PRIMARY)
    sensitivity = ires_deviation(plist, _WORKER_REF_SENSITIVITY, _WORKER_IRES_SENSITIVITY)
    result: Dict[str, object] = {
        "candidate_id": candidate.candidate_id,
        "ensemble": candidate.ensemble,
        "generation_index": candidate.generation_index,
        "sequence_md5": md5(candidate.mature_dna),
        "cds_md5": md5(candidate.cds_dna),
        "mature_length_nt": len(candidate.mature_dna),
        "cds_length_nt": len(candidate.cds_dna),
        "circular_mfe_kcal_mol": float(mfe),
        "circular_mfe_per_100nt": float(mfe) / len(candidate.mature_dna) * 100.0,
        "ensemble_free_energy_kcal_mol": ensemble_free_energy,
        "human_cai": human_cai(candidate.cds_dna),
        "mature_gc_fraction": gc_fraction(candidate.mature_dna),
        "cds_gc_fraction": gc_fraction(candidate.cds_dna),
        "cds_gc3_fraction": gc3_fraction(candidate.cds_dna),
        "bpp_cutoff": _WORKER_CUTOFF,
        "fold_elapsed_seconds": time.time() - t0,
        # `structure` is computed for every candidate; this used to keep it only for
        # the original, so any later structural analysis needed a full 261-sequence
        # refold. 2 KB per row is negligible -- persist it for all of them.
        "mfe_dot_bracket": structure,
    }
    result.update(primary)
    result.update({f"{key}_ires669": value for key, value in sensitivity.items()})
    return result


def write_checkpoint(rows: List[Mapping[str, object]], path: Path) -> None:
    df = pd.DataFrame(rows).sort_values(["generation_index", "candidate_id"])
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> Tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def empirical_percentile(values: np.ndarray, current: float, higher_better: bool) -> Dict[str, float]:
    if higher_better:
        worse = int(np.sum(values < current))
        ties = int(np.sum(values == current))
    else:
        worse = int(np.sum(values > current))
        ties = int(np.sum(values == current))
    # Mid-rank for ties; Wilson interval uses conservative integer worse count.
    percentile = (worse + 0.5 * ties) / len(values)
    # Use the same mid-rank convention for the uncertainty interval.  This is
    # important for the composition-shuffle ensemble, whose CAI values are
    # intentionally identical to the current construct.
    low, high = wilson_interval(worse + 0.5 * ties, len(values))
    return {
        "percentile": percentile * 100.0,
        "ci95_low": low * 100.0,
        "ci95_high": high * 100.0,
        "n": int(len(values)),
        "worse": worse,
        "ties": ties,
    }


def dominance_layers(df: pd.DataFrame) -> Tuple[Dict[str, int], Dict[str, List[str]]]:
    ids = list(df["candidate_id"])
    benefit = np.column_stack(
        [
            -df["circular_mfe_kcal_mol"].to_numpy(float),
            df["human_cai"].to_numpy(float),
            -df["lires"].to_numpy(float),
        ]
    )
    dominates: Dict[int, List[int]] = {i: [] for i in range(len(ids))}
    dominated_count = np.zeros(len(ids), dtype=int)
    dominators: Dict[str, List[str]] = {candidate_id: [] for candidate_id in ids}
    for i in range(len(ids)):
        ge = np.all(benefit >= benefit[i], axis=1)
        gt = np.any(benefit > benefit[i], axis=1)
        dom_i = np.where(ge & gt)[0]
        dominated_count[i] = len(dom_i)
        dominators[ids[i]] = [ids[j] for j in dom_i]
        for j in np.where(np.all(benefit[i] >= benefit, axis=1) & np.any(benefit[i] > benefit, axis=1))[0]:
            dominates[i].append(int(j))

    remaining = set(range(len(ids)))
    layers: Dict[str, int] = {}
    layer = 1
    while remaining:
        front = []
        for i in remaining:
            is_dominated = False
            for j in remaining:
                if i == j:
                    continue
                if np.all(benefit[j] >= benefit[i]) and np.any(benefit[j] > benefit[i]):
                    is_dominated = True
                    break
            if not is_dominated:
                front.append(i)
        if not front:
            raise RuntimeError("Pareto layer calculation stalled")
        for i in front:
            layers[ids[i]] = layer
            remaining.remove(i)
        layer += 1
    return layers, dominators


def metric_summary(df: pd.DataFrame, column: str) -> Dict[str, float]:
    arr = df[column].to_numpy(float)
    return {
        "min": float(np.min(arr)),
        "q1": float(np.quantile(arr, 0.25)),
        "median": float(np.median(arr)),
        "q3": float(np.quantile(arr, 0.75)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "sd": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
    }


def make_figures(df: pd.DataFrame, summary: Mapping[str, object], output_dir: Path) -> None:
    random_df = df[df["ensemble"].isin(["human_weighted", "uniform", "composition_shuffle"])]
    original = df.loc[df["candidate_id"] == "original"].iloc[0]
    colors = {
        "human_weighted": "#2A6FBB",
        "uniform": "#E08B2C",
        "composition_shuffle": "#4B9B69",
    }
    labels = {
        "human_weighted": "Human-usage weighted",
        "uniform": "Uniform synonymous",
        "composition_shuffle": "Original-composition shuffle",
    }

    metrics = [
        ("circular_mfe_kcal_mol", "Circular MFE (kcal/mol)", "Lower is better"),
        ("human_cai", "Human CAI", "Higher is better"),
        ("lires", r"IRES deviation $L_{IRES}$", "Lower is better"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.1), constrained_layout=True)
    for ax, (column, label, direction) in zip(axes, metrics):
        for ensemble in colors:
            values = random_df.loc[random_df["ensemble"] == ensemble, column]
            ax.hist(
                values,
                bins=18,
                alpha=0.42,
                color=colors[ensemble],
                label=labels[ensemble],
                edgecolor="white",
                linewidth=0.4,
            )
        ax.axvline(original[column], color="#B00020", linewidth=2.2, label="Current Combined")
        ax.set_xlabel(label)
        ax.set_ylabel("Candidate count")
        ax.set_title(direction, fontsize=10, color="#555555")
        ax.grid(axis="y", alpha=0.18)
    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.06))
    for ext in ("png", "svg"):
        fig.savefig(output_dir / f"figure_1_metric_distributions.{ext}", dpi=260, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.7, 5.7), constrained_layout=True)
    scatter = ax.scatter(
        random_df["circular_mfe_kcal_mol"],
        random_df["human_cai"],
        c=random_df["lires"],
        cmap="viridis_r",
        s=38,
        alpha=0.76,
        linewidths=0.25,
        edgecolors="white",
    )
    controls = df[df["ensemble"] == "bounding_control"]
    ax.scatter(
        controls["circular_mfe_kcal_mol"],
        controls["human_cai"],
        marker="x",
        color="#333333",
        s=70,
        linewidths=1.3,
        label="Bounding controls",
    )
    ax.scatter(
        [original["circular_mfe_kcal_mol"]],
        [original["human_cai"]],
        marker="*",
        color="#D00000",
        s=250,
        edgecolors="white",
        linewidths=0.8,
        zorder=5,
        label="Current Combined",
    )
    ax.set_xlabel("Circular MFE (kcal/mol)  ← more stable")
    ax.set_ylabel("Human CAI  → better adapted")
    ax.grid(alpha=0.18)
    cb = fig.colorbar(scatter, ax=ax)
    cb.set_label(r"IRES deviation $L_{IRES}$  (lower is better)")
    ax.legend(frameon=False, loc="best")
    for ext in ("png", "svg"):
        fig.savefig(output_dir / f"figure_2_pareto_overview.{ext}", dpi=260, bbox_inches="tight")
    plt.close(fig)

    ensemble_order = ["all_random", "human_weighted", "uniform", "composition_shuffle"]
    display = ["All random", "Human weighted", "Uniform", "Composition shuffle"]
    percentile_map = summary["percentiles"]
    matrix = np.array(
        [
            [
                percentile_map[name]["circular_mfe"]["percentile"],
                percentile_map[name]["human_cai"]["percentile"],
                percentile_map[name]["ires_deviation"]["percentile"],
            ]
            for name in ensemble_order
        ]
    )
    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    x = np.arange(len(display))
    width = 0.24
    bar_colors = ["#315E9A", "#4B9B69", "#9A5A9A"]
    names = ["Circular MFE", "Human CAI", "IRES integrity"]
    for idx in range(3):
        bars = ax.bar(x + (idx - 1) * width, matrix[:, idx], width, label=names[idx], color=bar_colors[idx])
        ax.bar_label(bars, fmt="%.0f", padding=2, fontsize=8)
    ax.axhline(50, color="#666666", linewidth=0.8, linestyle="--")
    ax.set_ylim(0, 108)
    ax.set_ylabel("Empirical desirability percentile")
    ax.set_xticks(x, display)
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.01))
    ax.grid(axis="y", alpha=0.18)
    for ext in ("png", "svg"):
        fig.savefig(output_dir / f"figure_3_percentiles.{ext}", dpi=260, bbox_inches="tight")
    plt.close(fig)


def build_summary(
    df: pd.DataFrame,
    construct: Mapping[str, object],
    args: argparse.Namespace,
    reference_info: Mapping[str, object],
) -> Dict[str, object]:
    original = df.loc[df["candidate_id"] == "original"].iloc[0]
    ensemble_frames = {
        "all_random": df[df["ensemble"].isin(["human_weighted", "uniform", "composition_shuffle"])],
        "human_weighted": df[df["ensemble"] == "human_weighted"],
        "uniform": df[df["ensemble"] == "uniform"],
        "composition_shuffle": df[df["ensemble"] == "composition_shuffle"],
    }
    percentiles: Dict[str, Dict[str, Dict[str, float]]] = {}
    ensemble_descriptives: Dict[str, Dict[str, Dict[str, float]]] = {}
    for name, subset in ensemble_frames.items():
        percentiles[name] = {
            "circular_mfe": empirical_percentile(
                subset["circular_mfe_kcal_mol"].to_numpy(float),
                float(original["circular_mfe_kcal_mol"]),
                higher_better=False,
            ),
            "human_cai": empirical_percentile(
                subset["human_cai"].to_numpy(float),
                float(original["human_cai"]),
                higher_better=True,
            ),
            "ires_deviation": empirical_percentile(
                subset["lires"].to_numpy(float),
                float(original["lires"]),
                higher_better=False,
            ),
            "ires_deviation_669": empirical_percentile(
                subset["lires_ires669"].to_numpy(float),
                float(original["lires_ires669"]),
                higher_better=False,
            ),
        }
        ensemble_descriptives[name] = {
            column: metric_summary(subset, column)
            for column in ("circular_mfe_kcal_mol", "human_cai", "lires", "lires_ires669")
        }

    pareto_input = df.copy()
    layers, dominators = dominance_layers(pareto_input)
    pareto_input_669 = df.copy()
    pareto_input_669["lires"] = pareto_input_669["lires_ires669"]
    layers_669, dominators_669 = dominance_layers(pareto_input_669)
    df["pareto_layer"] = df["candidate_id"].map(layers)
    df["dominators_count"] = df["candidate_id"].map(lambda x: len(dominators[x]))
    original_dominators = dominators["original"]
    original_dominators_669 = dominators_669["original"]
    original_layer = layers["original"]

    p_all = percentiles["all_random"]
    core_percentiles = [
        p_all["circular_mfe"]["percentile"],
        p_all["human_cai"]["percentile"],
        p_all["ires_deviation"]["percentile"],
    ]
    random_n = len(ensemble_frames["all_random"])
    random_dominators = [
        candidate_id
        for candidate_id in original_dominators
        if str(df.loc[df["candidate_id"] == candidate_id, "ensemble"].iloc[0])
        in {"human_weighted", "uniform", "composition_shuffle"}
    ]
    random_dominators_669 = [
        candidate_id
        for candidate_id in original_dominators_669
        if str(df.loc[df["candidate_id"] == candidate_id, "ensemble"].iloc[0])
        in {"human_weighted", "uniform", "composition_shuffle"}
    ]
    domination_fraction = len(random_dominators) / random_n if random_n else float("nan")
    domination_fraction_669 = len(random_dominators_669) / random_n if random_n else float("nan")

    dominator_details = []
    for candidate_id in original_dominators:
        row = df.loc[df["candidate_id"] == candidate_id].iloc[0]
        dominator_details.append(
            {
                "candidate_id": candidate_id,
                "ensemble": str(row["ensemble"]),
                "circular_mfe_kcal_mol": float(row["circular_mfe_kcal_mol"]),
                "human_cai": float(row["human_cai"]),
                "lires": float(row["lires"]),
                "lires_ires669": float(row["lires_ires669"]),
                "mfe_improvement_percent": (
                    float(original["circular_mfe_kcal_mol"]) - float(row["circular_mfe_kcal_mol"])
                )
                / abs(float(original["circular_mfe_kcal_mol"]))
                * 100.0,
                "lires_improvement_percent": (
                    float(original["lires"]) - float(row["lires"])
                )
                / float(original["lires"])
                * 100.0,
            }
        )

    # Decision rubric is explicit and project-defined; it is not an official
    # circDesign threshold.  It prevents post-hoc wording from being adjusted to
    # force a favorable answer.
    if (
        float(original["human_cai"]) >= 0.80
        and min(core_percentiles) >= 25.0
        and sum(p >= 50.0 for p in core_percentiles) >= 2
        and domination_fraction <= 0.05
    ):
        decision = "retain_supported"
        decision_cn = "计算结果支持保留现有 Combined 构建"
    elif (
        float(original["human_cai"]) >= 0.75
        and min(core_percentiles) >= 15.0
        and domination_fraction <= 0.10
    ):
        decision = "retain_with_tradeoff"
        decision_cn = "可在明确多目标权衡的前提下保留现有 Combined 构建"
    else:
        decision = "redesign_or_targeted_optimization"
        decision_cn = "现有计算证据不足以直接支持保留，建议定向同义优化后复评"

    ires_boundary_rank_delta = abs(
        p_all["ires_deviation"]["percentile"] - p_all["ires_deviation_669"]["percentile"]
    )
    boundary_robust = ires_boundary_rank_delta <= 10.0

    return {
        "analysis_title": "Combined mature circRNA retrospective circDesign-inspired evaluation",
        "analysis_timestamp_local": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "input_file": str(args.input.name),
        "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
        "software": {
            "python": sys.version.split()[0],
            "ViennaRNA": RNA.__version__,
            "Biopython": __import__("Bio").__version__,
            "random_seed": args.seed,
            "bpp_cutoff": args.bpp_cutoff,
            "workers": args.workers,
        },
        "coordinate_model": {
            "plasmid_length_bp": len(construct["plasmid"]),
            "mature_circle_plasmid_coordinates_1based_inclusive": [MATURE_START, MATURE_END],
            "mature_circle_length_nt": len(construct["mature"]),
            "junction_32nt_each_side": construct["junction"],
            "cds_plasmid_coordinates_1based_inclusive": [CDS_START, CDS_END],
            "cds_mature_coordinates_1based_inclusive": list(relative_interval(CDS_START, CDS_END)),
            "cds_length_nt": len(construct["cds"]),
            "protein_length_aa": len(construct["protein"]),
            "internal_stop_count": 0,
            "ires_primary_plasmid_coordinates": [IRES_START, IRES_END_PRIMARY],
            "ires_primary_mature_coordinates": list(relative_interval(IRES_START, IRES_END_PRIMARY)),
            "ires_primary_length_nt": len(construct["ires_primary"]),
            "ires_sensitivity_plasmid_coordinates": [IRES_START, IRES_END_SENSITIVITY],
            "ires_sensitivity_mature_coordinates": list(relative_interval(IRES_START, IRES_END_SENSITIVITY)),
            "ires_sensitivity_length_nt": len(construct["ires_sensitivity"]),
        },
        "sequence_checksums": {
            "plasmid_md5": md5(str(construct["plasmid"])),
            "mature_circle_md5": md5(str(construct["mature"])),
            "cds_md5": md5(str(construct["cds"])),
            "protein_md5": md5(str(construct["protein"])),
        },
        "benchmark_design": {
            "human_weighted_n": args.n_human,
            "uniform_synonymous_n": args.n_uniform,
            "composition_shuffle_n": args.n_shuffle,
            "random_total_n": random_n,
            "bounding_controls": list(df.loc[df["ensemble"] == "bounding_control", "candidate_id"]),
            "all_candidates_with_original_n": len(df),
        },
        "ires_references": reference_info,
        "current_metrics": {
            key: (float(original[key]) if isinstance(original[key], (float, np.floating, int, np.integer)) else original[key])
            for key in (
                "circular_mfe_kcal_mol",
                "circular_mfe_per_100nt",
                "ensemble_free_energy_kcal_mol",
                "human_cai",
                "mature_gc_fraction",
                "cds_gc_fraction",
                "cds_gc3_fraction",
                "lires",
                "lires_internal",
                "lires_crosstalk",
                "ires_cross_pair_probability_mass",
                "lires_ires669",
                "lires_internal_ires669",
                "lires_crosstalk_ires669",
            )
        },
        "percentiles": percentiles,
        "ensemble_descriptives": ensemble_descriptives,
        "pareto": {
            "current_layer": original_layer,
            "current_on_front": original_layer == 1,
            "all_dominators_count": len(original_dominators),
            "all_dominators": original_dominators,
            "random_dominators_count": len(random_dominators),
            "random_dominators": random_dominators,
            "random_domination_fraction": domination_fraction,
            "dominator_details": dominator_details,
        },
        "pareto_sensitivity_669nt": {
            "current_layer": layers_669["original"],
            "current_on_front": layers_669["original"] == 1,
            "all_dominators_count": len(original_dominators_669),
            "all_dominators": original_dominators_669,
            "random_dominators_count": len(random_dominators_669),
            "random_dominators": random_dominators_669,
            "random_domination_fraction": domination_fraction_669,
        },
        "boundary_sensitivity": {
            "ires_integrity_percentile_666nt": p_all["ires_deviation"]["percentile"],
            "ires_integrity_percentile_669nt": p_all["ires_deviation_669"]["percentile"],
            "absolute_percentile_change": ires_boundary_rank_delta,
            "robust_within_10_percentile_points": boundary_robust,
        },
        "decision": {
            "code": decision,
            "conclusion_cn": decision_cn,
            "rubric": {
                "strong_retain": "CAI >= 0.80; every core percentile >= 25; at least two >= 50; <=5% random candidates dominate",
                "retain_with_tradeoff": "CAI >= 0.75; every core percentile >= 15; <=10% random candidates dominate",
                "status": "project-defined retrospective rule, not an official circDesign threshold",
            },
        },
        "limitations": [
            "This is an independent circDesign-inspired retrospective implementation, not the unpublished original production code.",
            "MFE and base-pairing probabilities are thermodynamic predictions and do not directly prove cellular half-life or translation.",
            "CAI uses the Kazusa general Homo sapiens codon-usage table, not a tissue-specific tRNA model.",
            "Base-pairing probabilities below the recorded cutoff are omitted from the sparse L2 calculation; their squared contribution is negligible but nonzero.",
            "The mature-circle boundary is inferred from retained-scar annotations and should be confirmed by junction sequencing for experimental identity.",
        ],
    }


def write_sequence_outputs(
    output_dir: Path,
    construct: Mapping[str, object],
    candidates: Sequence[Candidate],
    df: pd.DataFrame,
) -> None:
    (output_dir / "combined_mature_circRNA.fasta").write_text(
        ">Combined_mature_circRNA|plasmid_723-2735|junction_2735_to_723|length_2013_nt\n"
        + wrap_fasta(str(construct["mature"]))
        + "\n",
        encoding="utf-8",
    )
    (output_dir / "combined_CDS.fasta").write_text(
        ">Combined_CDS|plasmid_1504-2586|length_1083_nt\n"
        + wrap_fasta(str(construct["cds"]))
        + "\n",
        encoding="utf-8",
    )
    (output_dir / "combined_protein.fasta").write_text(
        ">Combined_antigen|length_360_aa\n" + wrap_fasta(str(construct["protein"])) + "\n",
        encoding="utf-8",
    )
    candidate_by_id = {candidate.candidate_id: candidate for candidate in candidates}
    with (output_dir / "synonymous_candidate_CDS.fasta").open("w", encoding="utf-8") as handle:
        for candidate_id in df.sort_values(["generation_index", "candidate_id"])["candidate_id"]:
            candidate = candidate_by_id[candidate_id]
            handle.write(
                f">{candidate.candidate_id}|ensemble={candidate.ensemble}|CAI={human_cai(candidate.cds_dna):.6f}\n"
            )
            handle.write(wrap_fasta(candidate.cds_dna) + "\n")

    original_row = df.loc[df["candidate_id"] == "original"].iloc[0]
    dbn = (
        ">Combined_mature_circRNA circular MFE structure\n"
        + str(construct["mature"]).replace("T", "U")
        + "\n"
        + str(original_row["mfe_dot_bracket"])
        + f" ({float(original_row['circular_mfe_kcal_mol']):.2f} kcal/mol)\n"
    )
    (output_dir / "combined_circular_MFE_structure.dbn").write_text(dbn, encoding="utf-8")


def write_readme(output_dir: Path, summary: Mapping[str, object]) -> None:
    current = summary["current_metrics"]
    p = summary["percentiles"]["all_random"]
    pareto = summary["pareto"]
    decision = summary["decision"]
    text = f"""# Combined circRNA circDesign-inspired retrospective audit

## Headline result

{decision['conclusion_cn']}

The reconstructed mature circle is 2,013 nt (plasmid positions 723-2735;
junction 2735->723).  Its 1,083-nt CDS encodes the intended 360-aa protein
without internal stops.

| Metric | Current Combined | Empirical desirability percentile (all random synonymous candidates) |
|---|---:|---:|
| Circular MFE | {current['circular_mfe_kcal_mol']:.2f} kcal/mol | {p['circular_mfe']['percentile']:.1f}% |
| Human CAI | {current['human_cai']:.4f} | {p['human_cai']['percentile']:.1f}% |
| IRES deviation, 666-nt boundary | {current['lires']:.4f} | {p['ires_deviation']['percentile']:.1f}% |

Pareto layer: {pareto['current_layer']}; synonymous random candidates that
simultaneously dominate the current construct in all three objectives:
{pareto['random_dominators_count']} ({pareto['random_domination_fraction'] * 100:.2f}%).

## Reproduce

Install Python 3.12, ViennaRNA 2.7.0, Biopython 1.85, NumPy, pandas,
SciPy, and Matplotlib, then run:

```bash
python combined_circdesign_analysis.py \
  --input "source_vector/111 YRH pUC57-ALL4-4xmC3dP28 vector 1101(1).dna" \
  --output . --workers 8
```

The random seed, coordinate model, software versions, cutoffs, sequence
checksums, percentile confidence intervals, and limitations are recorded in
`summary.json`.

## Interpretation boundary

This analysis supports or challenges a sequence-retention decision under the
published circDesign objectives.  It is not an official circDesign
certification and does not replace experimental measurements of circularization,
RNA half-life, or protein output.
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Combined SnapGene .dna file")
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--n-human", type=int, default=96)
    parser.add_argument("--n-uniform", type=int, default=96)
    parser.add_argument("--n-shuffle", type=int, default=64)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--bpp-cutoff", type=float, default=DEFAULT_BPP_CUTOFF)
    parser.add_argument("--fresh", action="store_true", help="Discard checkpoint and recompute")
    args = parser.parse_args()

    args.input = args.input.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    checkpoint = args.output / "candidate_metrics_checkpoint.csv"
    final_metrics = args.output / "candidate_metrics.csv"
    if args.fresh:
        checkpoint.unlink(missing_ok=True)
        final_metrics.unlink(missing_ok=True)

    print("STAGE extract_construct", flush=True)
    construct = extract_construct(args.input)
    pd.DataFrame(construct["features"]).to_csv(args.output / "mature_circle_feature_map.csv", index=False)
    candidates = build_candidates(
        str(construct["mature"]),
        str(construct["cds"]),
        str(construct["protein"]),
        args.n_human,
        args.n_uniform,
        args.n_shuffle,
        args.seed,
    )
    print(
        f"STAGE candidates total={len(candidates)} human={args.n_human} uniform={args.n_uniform} shuffle={args.n_shuffle}",
        flush=True,
    )

    full_rna = str(construct["mature"]).replace("T", "U")
    ires_primary = relative_interval(IRES_START, IRES_END_PRIMARY)
    ires_sensitivity = relative_interval(IRES_START, IRES_END_SENSITIVITY)
    print(f"STAGE reference_primary interval={ires_primary}", flush=True)
    ref_primary, ref_primary_info = constrained_ires_reference(full_rna, ires_primary, args.bpp_cutoff)
    print(
        f"REFERENCE primary pairs={len(ref_primary)} elapsed={ref_primary_info['elapsed_seconds']:.1f}s",
        flush=True,
    )
    print(f"STAGE reference_sensitivity interval={ires_sensitivity}", flush=True)
    ref_sensitivity, ref_sensitivity_info = constrained_ires_reference(
        full_rna, ires_sensitivity, args.bpp_cutoff
    )
    print(
        f"REFERENCE sensitivity pairs={len(ref_sensitivity)} elapsed={ref_sensitivity_info['elapsed_seconds']:.1f}s",
        flush=True,
    )

    rows: List[Dict[str, object]] = []
    if checkpoint.exists():
        existing = pd.read_csv(checkpoint)
        rows = existing.to_dict(orient="records")
        print(f"RESUME checkpoint_rows={len(rows)}", flush=True)
    completed_ids = {str(row["candidate_id"]) for row in rows}
    pending = [candidate for candidate in candidates if candidate.candidate_id not in completed_ids]
    print(f"STAGE fold pending={len(pending)} workers={args.workers}", flush=True)
    t_fold = time.time()
    if pending:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            initializer=init_worker,
            initargs=(
                ref_primary,
                ref_sensitivity,
                ires_primary,
                ires_sensitivity,
                args.bpp_cutoff,
            ),
        ) as executor:
            future_map = {executor.submit(evaluate_candidate, candidate): candidate for candidate in pending}
            finished_this_run = 0
            for future in as_completed(future_map):
                candidate = future_map[future]
                try:
                    row = future.result()
                except Exception as exc:
                    print(f"ERROR candidate={candidate.candidate_id} error={exc!r}", flush=True)
                    for other in future_map:
                        other.cancel()
                    raise
                rows.append(row)
                finished_this_run += 1
                if finished_this_run % args.workers == 0 or finished_this_run == len(pending):
                    write_checkpoint(rows, checkpoint)
                    print(
                        f"PROGRESS completed={len(rows)}/{len(candidates)} elapsed={time.time() - t_fold:.1f}s last={candidate.candidate_id}",
                        flush=True,
                    )

    df = pd.DataFrame(rows)
    expected = {candidate.candidate_id for candidate in candidates}
    actual = set(df["candidate_id"].astype(str))
    if expected != actual:
        raise RuntimeError(f"Candidate results mismatch: missing={sorted(expected-actual)} extra={sorted(actual-expected)}")
    df = df.drop_duplicates("candidate_id", keep="last").sort_values(["generation_index", "candidate_id"])

    reference_info = {
        "primary_666nt": ref_primary_info,
        "sensitivity_669nt": ref_sensitivity_info,
        "definition": "circular full-length partition function with every non-IRES nucleotide forced unpaired",
    }
    summary = build_summary(df, construct, args, reference_info)
    layers, dominators = dominance_layers(df)
    df["pareto_layer"] = df["candidate_id"].map(layers)
    df["dominators_count"] = df["candidate_id"].map(lambda x: len(dominators[x]))
    df.to_csv(final_metrics, index=False)
    checkpoint.unlink(missing_ok=True)

    write_sequence_outputs(args.output, construct, candidates, df)
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    make_figures(df, summary, args.output)
    write_readme(args.output, summary)
    print(
        "RESULT "
        + json.dumps(
            {
                "decision": summary["decision"],
                "current_metrics": summary["current_metrics"],
                "percentiles_all_random": summary["percentiles"]["all_random"],
                "pareto": summary["pareto"],
                "boundary_sensitivity": summary["boundary_sensitivity"],
            },
            ensure_ascii=False,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
