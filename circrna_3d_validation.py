#!/usr/bin/env python3
"""Validate the 2-D circDesign-style prediction against the delivered 3-D model.

Why this exists
---------------
The retrospective report evaluates the construct with three 2-D proxies:
circular MFE, CAI, and an L2 distance between base-pair-probability matrices
(``L_IRES``).  A composition-matched null analysis showed that ``L_IRES`` is not
predictable from any codon-level feature: cross-validated R^2 is about 0.02 for
codon usage and the arrangement block adds no measurable increment (the fit excludes
the two rows carrying the partition-function sentinel; including them moves the value
to -0.11 for reasons unrelated to the biology).  A quantity the sequence does not
determine cannot serve as an optimisation objective.

Meanwhile a complete all-atom 3-D model of the SAME 2013 nt circle already exists
(TorusFold-Hybrid, ``artifacts/2013nt/isrnaclong_final.pdb``, OpenMM 8.5.2).  This
script replaces the proxy with a measurement:

  * derive the base pairs actually realised in the 3-D coordinates
  * compare them, pair by pair, against the 2-D circular MFE dot-bracket
  * report the agreement separately for the IRES region, the CDS, and the rest

Outputs ``circrna_3d_validation.json``.  Read-only with respect to TorusFold.

H-bond criteria (heavy-atom donor/acceptor distances, standard geometries):
    A-U   A:N1-U:N3, A:N6-U:O4
    G-C   G:N1-C:N3, G:O6-C:N4, G:N2-C:O2
    G-U   G:N1-U:O2, G:O6-U:N3
A pair counts as a base pair when |i-j| >= 4 (no self/adjacent stacking), the
base reference atoms (N9 for purines, N1 for pyrimidines) are within 12 A, and at
least two of the listed H-bonds are within HBOND_MAX.

Run:  python circrna_3d_validation.py [--pdb PATH]
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
DEFAULT_PDB = Path("C:/baidunetdiskdownload/torusfold-hybrid/artifacts/2013nt/isrnaclong_final.pdb")
DBN = ROOT / "combined_circular_MFE_structure.dbn"
METRICS = ROOT / "candidate_metrics.csv"
OUT = ROOT / "circrna_3d_validation.json"

HBOND_MAX = 3.6          # A, donor-acceptor heavy-atom distance
BASE_REF_MAX = 12.0      # A, base-reference atom separation
MIN_HBONDS = 2

PURINES = {"A", "G"}
BASE_REF_ATOM = {"A": "N9", "G": "N9", "C": "N1", "U": "N1"}

HBOND_PAIRS = {
    ("A", "U"): [("N1", "N3"), ("N6", "O4")],
    ("U", "A"): [("N3", "N1"), ("O4", "N6")],
    ("G", "C"): [("N1", "N3"), ("O6", "N4"), ("N2", "O2")],
    ("C", "G"): [("N3", "N1"), ("N4", "O6"), ("O2", "N2")],
    ("G", "U"): [("N1", "O2"), ("O6", "N3")],
    ("U", "G"): [("O2", "N1"), ("N3", "O6")],
}

# mature-circle coordinates from summary.json
IRES_PRIMARY = (101, 766)     # 666 nt
IRES_SENS = (101, 769)        # 669 nt
CDS = (782, 1864)             # 1083 nt


def parse_pdb(path: Path):
    """Return per-residue base identity and a {residue_index: {atom_name: xyz}} map."""
    atoms: dict[int, dict[str, np.ndarray]] = defaultdict(dict)
    bases: dict[int, str] = {}
    with path.open() as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                resname = line[17:20].strip()
                resseq = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            _ = serial
            base = resname[-1] if resname != "RA" else "A"
            base = {"RA": "A", "RC": "C", "RG": "G", "RU": "U"}.get(resname, base)
            bases[resseq] = base
            atoms[resseq][name] = xyz
    return bases, atoms


def load_dbn(path: Path):
    """Return (sequence, dot_bracket) from a plain or FASTA-ish .dbn file."""
    seq = None
    struct = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        if seq is None:
            seq = line.upper()
        elif struct is None:
            # The writer appends " (<mfe> kcal/mol)" to the structure line. Feeding the
            # whole line to the bracket parser lets the annotation's parentheses count:
            # observed effect is a phantom 658th pair at positions (2015, 2032), outside
            # a 2013 nt molecule. Truncate to the sequence length.
            struct = line[: len(seq)]
            break
    if seq is None or struct is None:
        raise SystemExit(f"could not parse .dbn: {path}")
    return seq, struct


def dbn_pairs(struct: str) -> set[tuple[int, int]]:
    """1-based (i, j) pairs from a dot-bracket string."""
    stack: list[int] = []
    pairs: set[tuple[int, int]] = set()
    for idx, ch in enumerate(struct, start=1):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                raise SystemExit("unbalanced dot-bracket")
            pairs.add((stack.pop(), idx))
    if stack:
        raise SystemExit("unbalanced dot-bracket")
    return pairs


def hbond_count(atoms_i, atoms_j, bi, bj) -> int:
    n = 0
    for a1, a2 in HBOND_PAIRS[(bi, bj)]:
        p = atoms_i.get(a1)
        q = atoms_j.get(a2)
        if p is None or q is None:
            continue
        if np.linalg.norm(p - q) <= HBOND_MAX:
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", type=Path, default=DEFAULT_PDB)
    args = ap.parse_args()

    bases, atoms = parse_pdb(args.pdb)
    n_res = len(bases)
    print(f"3-D model: {args.pdb.name}  residues={n_res}")

    seq3d = "".join(bases[i] for i in sorted(bases))
    seq2d, struct2d = load_dbn(DBN)
    print(f"2-D model: {DBN.name}  length={len(seq2d)}")

    if seq3d != seq2d:
        diff = [k for k, (a, b) in enumerate(zip(seq3d, seq2d)) if a != b]
        raise SystemExit(
            f"sequence mismatch between 3-D model and 2-D file: {len(diff)} positions "
            f"(first {diff[:5]}) -- indices are not comparable, refusing to proceed"
        )
    print("sequence check: 3-D and 2-D are the same 2013 nt sequence -- indices comparable")

    # base reference atoms, indexed 1..n
    ref = np.zeros((n_res + 1, 3))
    for i in range(1, n_res + 1):
        ref[i] = atoms[i][BASE_REF_ATOM[bases[i]]]
    D = np.linalg.norm(ref[1:, None, :] - ref[None, 1:, :], axis=-1)

    print("deriving base pairs from 3-D geometry...")
    pairs3d: set[tuple[int, int]] = set()
    for i in range(1, n_res + 1):
        for j in range(i + 4, n_res + 1):
            if D[i - 1, j - 1] > BASE_REF_MAX:
                continue
            bi, bj = bases[i], bases[j]
            if (bi, bj) not in HBOND_PAIRS:
                continue
            if hbond_count(atoms[i], atoms[j], bi, bj) >= MIN_HBONDS:
                pairs3d.add((i, j))
    print(f"  {len(pairs3d)} base pairs realised in 3-D")
    print(f"  {len(struct2d.replace('.', '')) // 2} base pairs in the 2-D MFE structure")

    pairs2d = dbn_pairs(struct2d)
    shared = pairs2d & pairs3d
    only2d = pairs2d - pairs3d
    only3d = pairs3d - pairs2d

    def region(i: int, lo: int, hi: int) -> bool:
        return lo <= i <= hi

    regions = {
        "IRES 666 (101-766)": IRES_PRIMARY,
        "IRES 669 (101-769)": IRES_SENS,
        "CDS (782-1864)": CDS,
    }
    per_region = {}
    for label, (lo, hi) in regions.items():
        p2 = {(i, j) for i, j in pairs2d if region(i, lo, hi) and region(j, lo, hi)}
        p3 = {(i, j) for i, j in pairs3d if region(i, lo, hi) and region(j, lo, hi)}
        sh = p2 & p3
        per_region[label] = {
            "2d_pairs": len(p2),
            "3d_pairs": len(p3),
            "agreed": len(sh),
            "sensitivity_3d_recovers_2d": (len(sh) / len(p2)) if p2 else None,
            "ppv_2d_pairs_realised_in_3d": (len(sh) / len(p3)) if p3 else None,
        }
        if p2:
            print(
                f"  {label:<22} 2D={len(p2):4d}  3D={len(p3):4d}  一致={len(sh):4d}  "
                f"3D 找回 2D 的 {100 * len(sh) / len(p2):5.1f}%"
            )

    overall = {
        "2d_pairs_total": len(pairs2d),
        "3d_pairs_total": len(pairs3d),
        "agreed": len(shared),
        "only_in_2d": len(only2d),
        "only_in_3d": len(only3d),
        "agreement_over_union_jaccard": len(shared) / (len(pairs2d | pairs3d) or 1),
    }
    print(
        f"\noverall: 2D={len(pairs2d)} 3D={len(pairs3d)} shared={len(shared)} "
        f"Jaccard={overall['agreement_over_union_jaccard']:.3f}"
    )

    # BSJ closure: distance between residue 1 and residue 2013 as modelled
    closed = float(np.linalg.norm(ref[1] - ref[n_res]))
    print(f"BSJ closure |P(1)-P(last)| over base-reference atoms: {closed:.3f} A")

    payload = {
        "what": "3-D realised base pairs vs 2-D circular MFE prediction, same 2013 nt construct",
        "source_3d": str(args.pdb),
        "source_2d": DBN.name,
        "sequence_identical": True,
        "criteria": {
            "hbond_max_angstrom": HBOND_MAX,
            "base_ref_max_angstrom": BASE_REF_MAX,
            "min_hbonds": MIN_HBONDS,
            "base_ref_atom": BASE_REF_ATOM,
        },
        "overall": overall,
        "per_region": per_region,
        "bsj_closure_angstrom": closed,
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {OUT.name}")


if __name__ == "__main__":
    main()
