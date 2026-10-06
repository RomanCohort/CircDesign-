#!/usr/bin/env python3
"""Independent verification pass. Computes raw quantities from raw data files only.

Deliberately does NOT read any conclusion fields (no ires_analysis.json,
design_panel.json, advantage_space.json, ires_metrics.json, no .md, and only the
coordinate_model + sequence_checksums fields of summary.json).
"""
import io
import json
import numpy as np
import pandas as pd

BASE = "."
OUT = {}

# ---------------------------------------------------------------- load table
df = pd.read_csv(f"{BASE}/candidate_metrics.csv")
random_ensembles = ["human_weighted", "uniform", "composition_shuffle"]

# ---------------------------------------------------------------- 1. counts
counts = df["ensemble"].value_counts().to_dict()
OUT["q1_total_rows"] = int(len(df))
OUT["q1_ensemble_counts"] = {k: int(v) for k, v in counts.items()}

# ---------------------------------------------------------------- 2. CAI max
orig_row = df[df["candidate_id"] == "original"].iloc[0]
orig_cai = float(orig_row["human_cai"])
subset = df[df["ensemble"].isin(random_ensembles)]
max_cai = float(subset["human_cai"].max())
n_tie = int((subset["human_cai"] == orig_cai).sum())
OUT["q2_orig_human_cai"] = orig_cai
OUT["q2_max_human_cai_among_random"] = max_cai
OUT["q2_orig_is_max"] = bool(orig_cai >= max_cai)
OUT["q2_n_candidates_tie_original"] = n_tie

# ---------------------------------------------------------------- 3. Pareto
def n_dominating(frame, o_mfe, o_cai, o_lires):
    mfe = frame["circular_mfe_kcal_mol"].to_numpy()
    cai = frame["human_cai"].to_numpy()
    lires = frame["lires"].to_numpy()
    dom = (
        (mfe <= o_mfe) & (cai >= o_cai) & (lires <= o_lires)
        & ((mfe < o_mfe) | (cai > o_cai) | (lires < o_lires))
    )
    return dom

o_mfe = float(orig_row["circular_mfe_kcal_mol"])
o_lires = float(orig_row["lires"])

dom_all = n_dominating(df, o_mfe, orig_cai, o_lires)
OUT["q3_dominators_all_261_count"] = int(dom_all.sum())
OUT["q3_dominators_all_261_ids"] = df.loc[dom_all, "candidate_id"].tolist()

sub3 = df[df["ensemble"].isin(["human_weighted", "uniform", "current_construct"])]
dom_sub = n_dominating(sub3, o_mfe, orig_cai, o_lires)
OUT["q3_dominators_restricted_count"] = int(dom_sub.sum())
OUT["q3_dominators_restricted_ids"] = sub3.loc[dom_sub, "candidate_id"].tolist()

# ---------------------------------------------------------------- 4. CAI >= orig
rand = df[df["ensemble"].isin(random_ensembles)]
pass_cai = rand["human_cai"] >= orig_cai
OUT["q4_random_rows"] = int(len(rand))
OUT["q4_pass_cai_ge_orig"] = int(pass_cai.sum())

# ---------------------------------------------------------------- 5. percentiles
shuf = df[df["ensemble"] == "composition_shuffle"]
OUT["q5_n_shuffle_rows"] = int(len(shuf))
OUT["q5_shuffles_worse_than_orig_mfe"] = int((shuf["circular_mfe_kcal_mol"] > o_mfe).sum())
OUT["q5_shuffles_worse_than_orig_lires"] = int((shuf["lires"] > o_lires).sum())

# ---------------------------------------------------------------- 6. effective rank
metrics = ["circular_mfe_kcal_mol", "human_cai", "lires"]
X = rand[metrics].to_numpy(dtype=float)
Z = (X - X.mean(axis=0)) / X.std(axis=0)
corr = np.corrcoef(Z, rowvar=False)
eig = np.linalg.eigvalsh(corr)[::-1]          # descending
pr = (eig.sum() ** 2) / (eig ** 2).sum()
OUT["q6_n_rows"] = int(len(rand))
OUT["q6_corr_matrix"] = corr.tolist()
OUT["q6_eigenvalues"] = eig.tolist()
OUT["q6_participation_ratio"] = float(pr)

# ---------------------------------------------------------------- 7. CDS fasta
def read_fasta(path):
    seqs = []
    header, buf = None, []
    with io.open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(buf)))
                header, buf = line[1:], []
            elif line:
                buf.append(line.strip())
    if header is not None:
        seqs.append((header, "".join(buf)))
    return seqs

CODON = {}
bases = "TCAG"
aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODON[a + b + c] = aas[i]
            i += 1

def translate(seq):
    return "".join(CODON.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))

syn = read_fasta(f"{BASE}/synonymous_candidate_CDS.fasta")
OUT["q7_n_sequences"] = len(syn)
OUT["q7_all_1083nt"] = bool(all(len(s) == 1083 for _, s in syn))
prot_set = set()
for _, s in syn:
    p = translate(s)
    if p.endswith("*"):
        p = p[:-1]
    prot_set.add(p)
OUT["q7_n_distinct_proteins"] = len(prot_set)
OUT["q7_all_translate_identical_360aa"] = bool(len(prot_set) == 1 and len(next(iter(prot_set))) == 360)
OUT["q7_protein_length"] = len(next(iter(prot_set)))

# ---------------------------------------------------------------- 8. repeats
cds = read_fasta(f"{BASE}/combined_CDS.fasta")[0][1]
n = len(cds)

def longest_direct_repeat(s):
    lo, hi = 1, len(s)
    best_len, best_pos = 0, None
    # binary search on length; track positions of best
    cache = {}
    def positions(L):
        if L not in cache:
            d = {}
            for i in range(len(s) - L + 1):
                sub = s[i:i + L]
                d.setdefault(sub, []).append(i)
            cache[L] = d
        return cache[L]

    # first find the max length
    def feasible(L):
        return any(len(v) >= 2 for v in positions(L).values())

    lo, hi = 1, n - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        if feasible(mid):
            best_len = mid
            lo = mid + 1
        else:
            hi = mid - 1
    if best_len == 0:
        return 0, None, None
    d = positions(best_len)
    # choose the pair with the smallest start positions
    best_sub = None
    for sub, plist in d.items():
        if len(plist) >= 2:
            plist = sorted(plist)
            if best_sub is None or plist[:2] < best_sub[1][:2]:
                best_sub = (sub, plist)
    sub, plist = best_sub
    return best_len, plist[0], plist[1]

rep_len, rep_p1, rep_p2 = longest_direct_repeat(cds)
OUT["q8_longest_direct_repeat_len"] = int(rep_len)
OUT["q8_longest_direct_repeat_positions"] = [int(rep_p1), int(rep_p2)]

def longest_homopolymer(s):
    best, bl, bp, cur, cp = 0, None, None, 0, 0
    for idx, ch in enumerate(s):
        if idx > 0 and ch == s[idx - 1]:
            cur += 1
        else:
            cur, cp = 1, idx
        if cur > best:
            best, bl, bp = cur, ch, cp
    return best, bl, bp

hp_len, hp_base, hp_pos = longest_homopolymer(cds)
OUT["q8_longest_homopolymer_len"] = int(hp_len)
OUT["q8_longest_homopolymer_base"] = hp_base
OUT["q8_longest_homopolymer_start"] = int(hp_pos)

# ---------------------------------------------------------------- 9. motifs
donor = ["GTAAGT", "GTGAGT", "GTATGT", "GTACGT"]
polya = ["AATAAA", "ATTAAA", "AATAAT", "CATAAA", "TATAAA"]

def count_motifs(s, motifs):
    return int(sum(s.count(m) for m in motifs))

OUT["q9_donor_motif_counts"] = {m: int(cds.count(m)) for m in donor}
OUT["q9_donor_motif_total"] = count_motifs(cds, donor)
OUT["q9_polya_motif_counts"] = {m: int(cds.count(m)) for m in polya}
OUT["q9_polya_motif_total"] = count_motifs(cds, polya)

# ---------------------------------------------------------------- 10. splice check
mature = read_fasta(f"{BASE}/combined_mature_circRNA.fasta")[0][1]
start, end = 782, 1864      # 1-based inclusive
spliced = mature[: start - 1] + cds + mature[end:]
OUT["q10_mature_length"] = len(mature)
OUT["q10_cds_length"] = len(cds)
OUT["q10_splice_reproduces_mature"] = bool(spliced == mature)
OUT["q10_region_equals_cds"] = bool(mature[start - 1:end] == cds)

# ---------------------------------------------------------------- emit
with io.open(f"{BASE}/verify_numbers.json", "w", encoding="utf-8") as fh:
    json.dump(OUT, fh, ensure_ascii=False, indent=2)

print(json.dumps(OUT, ensure_ascii=False, indent=2))
