# Combined circRNA — circDesign-inspired retrospective evaluation

A computational audit of a 2,013 nt mature circRNA construct (JLU-FBH CirCure),
plus the verification layer that audits the audit.

**This is an independent retrospective implementation of the published circDesign
objectives. It is not the unpublished original production code, and it is not an
official circDesign certification.**

---

## Headline result

The construct is retained for the current development chain. Three independently
supportable statements, each scoped to the space in which it holds:

| Statement | Space | Evidence |
|---|---|---|
| **Codon composition reaches the supremum** | 256 synonymous references | CAI = 0.8110 is the maximum; none of the 192 human-weighted and uniform candidates exceeds it |
| **Top 1.9% of buildable designs** | 1,553 correctness-passing synonymous rearrangements | 30 dominate it on the three objectives; 95% CI [1.36%, 2.74%] |
| **No systematic improvement exists** | all synonymous rearrangements | domination rate 1.45% (58/4,000); permutation test puts the observed count at chance |

And one statement that is **not** supportable:

> **The construct is not the absolute optimum in the synonymous arrangement space.**
> 5.50% of rearrangements have a lower MFE and 1.45% are not worse on all three
> objectives. These are enumerable, reproducible counterexamples. Any claim that the
> sequence is "optimal" holds only when the space and objective set are named.

---

## What was corrected

The original delivered report measured its three headline percentiles against a
**pooled** reference distribution of 96 human-weighted + 96 uniform + 64
composition-shuffle sequences. That is methodologically unsound: the first two
groups vary codon *composition*, the third holds composition exactly fixed and
varies only *arrangement*. A pooled distribution is dominated by composition
variance and pulls every metric toward values favourable to the construct.

| Issue | Original | Corrected |
|---|---|---|
| Reference distribution | one pooled n=256 | per-metric: composition axis vs arrangement axis |
| MFE percentile | 95.7% (pooled) | **94.50%** (composition-matched, n=4,000) |
| CAI percentile | 87.5% (pooled) | **100% of 192** composition candidates |
| Significance test | sign test (H0: at the median) | **exchangeability** (H0: exchangeable with its own codon permutations) |
| IRES metric | `L_IRES`, an L2 distance | four interpretable ratios; all non-significant, p = 0.29–0.71 |
| Domination denominator | 5/256 | **30/1,553** (192 of the 256 were algebraically ineligible) |
| "Dominators are a defect" | Pareto layer 4 | **at the chance expectation** (expected 4.64, observed 5, P = 0.54) |

Three defects in the delivered data are documented rather than silently worked around:

- Two rows of `candidate_metrics.csv` carry a partition-function failure recorded as
  the sentinel `ensemble_free_energy_kcal_mol = 1e5` with an empty pair list. They flip
  the sign of a regression R² from −0.11 to +0.02. All statistics here exclude them.
- `candidate_metrics.csv` **cannot be reproduced by the delivered source** (their file
  timestamps are contradictory, and the current script cannot emit that sentinel).
  Every key result was therefore recomputed independently.
- The delivered 64-shuffle baseline is a **1-in-104 unlucky draw**. Verified three ways:
  the fold is bit-identical on all 64, the sampler reproduces all 64 from the same seed,
  and the MFE distributions are KS-indistinguishable — yet 9/64 beat the construct
  against a true rate of 5.50%.

A `composition_shuffle_055` second-generation-candidate recommendation in the original
report is **withdrawn**: its IRES improvement is z = −1.51 within the shuffle group and
does not survive multiple-comparison correction.

---

## Reproduce

```bash
pip install -r requirements.txt          # Python 3.12+, ViennaRNA 2.7.0

# original implementation (unmodified)
python combined_circdesign_analysis.py \
  --input "source_vector/111 YRH pUC57-ALL4-4xmC3dP28 vector 1101(1).dna" \
  --output . --workers 8

# verification layer
python axis_decomposition.py             # composition / arrangement axes, effective rank
python ires_metrics.py                   # four interpretable IRES ratios
python ires_analysis.py
python design_panel.py                   # 9 mechanism QC metrics, multiplicity-corrected
python boundary_test.py                  # boundary robustness redone
python advantage_space.py --n 2000       # pre-registered multi-objective test
python accessibility.py                  # translation accessibility objectives
python pool_advantage.py                 # pool independent runs
python build_report.py                   # regenerate the DOCX report
```

`multiprocessing` on Windows uses `spawn`: scripts splice sequences in the **parent**
and pass the finished molecule to workers. `accessibility.build_mature` raises if the
template is unset rather than silently folding the wrong molecule — that failure mode
reached a published number once and is now structurally blocked.

---

## Files

| Path | What |
|---|---|
| `combined_circdesign_analysis.py` | original implementation, unmodified |
| `build_report.py` | DOCX generator; all numbers read live from the JSON artefacts |
| `axis_decomposition.py` | composition (61-d codon usage) vs arrangement (3,721-d codon-pair residual) |
| `ires_metrics.py` / `ires_analysis.py` | replacement IRES ratios |
| `design_panel.py` | synthesis / recombination / processing QC panel |
| `advantage_space.py` | pre-registered multi-objective test, seed-specific outputs |
| `accessibility.py` / `accessibility_full.py` | translation accessibility objectives |
| `validate_seed.py` | out-of-sample replication on an unused seed |
| `refold_check.py` / `shuffle_repro_test.py` | fold and sampler collision tests |
| `verify_numbers.py` | independent recomputation, written without reading conclusions |
| `sequence_optimality_report.md` | optimality argument: provable scope and its boundary |

---

## Interpretation boundary

- Retention is supported by **the absence of any executable systematic improvement**,
  by **the composition axis reaching its supremum**, and by the construct having
  already completed cloning and prior experimental validation — not by any claim of
  global optimality.
- MFE is not cellular half-life; CAI is not translation efficiency; IRES structural
  retention is not IRES function. All three are computational proxies.
- The IRES dimension carries **no detectable sequence-readable signal** in this design
  (four metrics, all p = 0.29–0.71) and should not be used as an optimisation objective.
- The translation-accessibility objective group was added **after** the three-objective
  front was known. Enlarging an objective set mechanically reduces the number of
  dominators, so the resulting zero is a statement about that objective set, not
  evidence that accessibility disproves the competitors (preferential test p = 0.38).
  A folding-degree sub-effect that appeared in-sample did **not** replicate out of
  sample (p = 0.76) and is not reported as a finding.
- The mature-circle boundary is inferred from retained-scar annotations and should be
  confirmed by back-splice junction RT-PCR and Sanger sequencing.
