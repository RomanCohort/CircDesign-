#!/usr/bin/env python3
"""Can the DELIVERED source reproduce the DELIVERED composition_shuffle sequences?

Context: refold_check.py showed my fold is bit-identical to the delivered fold, yet the
construct's MFE percentile is 86th in the delivered 64-shuffle ensemble and 94th in a
fresh 2000-shuffle ensemble drawn with the same documented procedure. Two ensembles that
different cannot come from one distribution, and the fold is ruled out.

The audit also found that candidate_metrics.csv's mtime predates the analysis script's, so
the delivered CSV cannot have been produced by the delivered source. This script tests
that directly: it replays the delivered shuffle generator, exactly as written at
combined_circdesign_analysis.py:327-335 with the documented seed, and compares the 64
sequences it produces against the 64 the delivered CSV claims.

Writes shuffle_repro_test.json.
"""

from __future__ import annotations

import json
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SEED = 20260804
N_HUMAN, N_UNIFORM, N_SHUFFLE = 96, 96, 64

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
# human-weighted codon table, as used by weighted_choice in the analysis script
HUMAN_FREQ = {
    "F": {"TTT": 0.45, "TTC": 0.55}, "L": {"TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40},
    "I": {"ATT": 0.36, "ATC": 0.48, "ATA": 0.16}, "M": {"ATG": 1.0},
    "V": {"GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47},
    "S": {"TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.24},
    "P": {"CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11},
    "T": {"ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12},
    "A": {"GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11},
    "Y": {"TAT": 0.43, "TAC": 0.57}, "H": {"CAT": 0.41, "CAC": 0.59},
    "Q": {"CAA": 0.25, "CAG": 0.75}, "N": {"AAT": 0.46, "AAC": 0.54},
    "K": {"AAA": 0.42, "AAG": 0.58}, "D": {"GAT": 0.46, "GAC": 0.54},
    "E": {"GAA": 0.42, "GAG": 0.58}, "C": {"TGT": 0.45, "TGC": 0.55},
    "W": {"TGG": 1.0},
    "R": {"CGT": 0.08, "CGC": 0.18, "CGA": 0.11, "CGG": 0.20, "AGA": 0.21, "AGG": 0.21},
    "G": {"GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25},
}


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


def weighted_choice(rng, aa):
    table = HUMAN_FREQ[aa]
    r, acc = rng.random(), 0.0
    for codon, w in table.items():
        acc += w
        if r <= acc:
            return codon
    return list(table)[-1]


def synonymous_cds_from_choices(protein, choose):
    return "".join(choose(aa, i) for i, aa in enumerate(protein)) + "TGA"


def main() -> None:
    delivered = read_fasta(ROOT / "synonymous_candidate_CDS.fasta")
    orig = delivered["original"]
    protein = "".join(CODON_TABLE[orig[i : i + 3]] for i in range(0, len(orig) - 3, 3))
    print(f"protein length {len(protein)} aa")

    # replay the generator EXACTLY as delivered, including the `seen` set's contents
    seen = {orig}
    hw, un = [], []
    human_rng = random.Random(SEED + 101)
    while len(hw) < N_HUMAN:
        c = synonymous_cds_from_choices(protein, lambda aa, i: weighted_choice(human_rng, aa))
        if c in seen:
            continue
        seen.add(c)
        hw.append(c)
    uniform_rng = random.Random(SEED + 202)
    while len(un) < N_UNIFORM:
        c = synonymous_cds_from_choices(
            protein, lambda aa, i: uniform_rng.choice(list(HUMAN_FREQ[aa]))[0])
        if c in seen:
            continue
        seen.add(c)
        un.append(c)

    original_sense_codons = [orig[i : i + 3] for i in range(0, len(orig) - 3, 3)]
    aa_positions = defaultdict(list)
    for idx, aa in enumerate(protein):
        aa_positions[aa].append(idx)
    shuffle_rng = random.Random(SEED + 303)

    def shuffle_generator():
        codons = list(original_sense_codons)
        for positions in aa_positions.values():
            values = [codons[pos] for pos in positions]
            shuffle_rng.shuffle(values)
            for pos, value in zip(positions, values):
                codons[pos] = value
        return "".join(codons) + orig[-3:]

    replayed = []
    while len(replayed) < N_SHUFFLE:
        c = shuffle_generator()
        if c in seen:
            continue
        seen.add(c)
        replayed.append(c)

    # compare against what the delivered CSV/fasta actually contains
    report = {"n_human_reproduced": 0, "n_uniform_reproduced": 0, "n_shuffle_reproduced": 0}
    for i, c in enumerate(hw, 1):
        if delivered.get(f"human_weighted_{i:03d}") == c:
            report["n_human_reproduced"] += 1
    for i, c in enumerate(un, 1):
        if delivered.get(f"uniform_{i:03d}") == c:
            report["n_uniform_reproduced"] += 1
    for i, c in enumerate(replayed, 1):
        if delivered.get(f"composition_shuffle_{i:03d}") == c:
            report["n_shuffle_reproduced"] += 1

    print("\nreproduced from the DELIVERED source, vs the DELIVERED fasta:")
    print(f"  human_weighted       : {report['n_human_reproduced']}/{N_HUMAN}")
    print(f"  uniform              : {report['n_uniform_reproduced']}/{N_UNIFORM}")
    print(f"  composition_shuffle  : {report['n_shuffle_reproduced']}/{N_SHUFFLE}")

    sh = [delivered[f"composition_shuffle_{i:03d}"] for i in range(1, N_SHUFFLE + 1)]
    print(f"\n  any of the replayed 64 present verbatim in the delivered shuffle set: "
          f"{len(set(replayed) & set(sh))}")

    identical = report["n_shuffle_reproduced"] == N_SHUFFLE
    verdict = ("delivered source REPRODUCES the delivered shuffles" if identical else
               "delivered source CANNOT reproduce the delivered shuffles -- the two "
               "ensembles are different objects and every shuffle-based number in the "
               "report depends on an unreproducible artefact")
    print(f"\nVERDICT: {verdict}")
    report["verdict"] = verdict
    report["n_replayed_present_in_delivered"] = len(set(replayed) & set(sh))
    report["seed"] = SEED
    (ROOT / "shuffle_repro_test.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print("wrote shuffle_repro_test.json")


if __name__ == "__main__":
    main()
