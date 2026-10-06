# Combined circRNA circDesign-inspired retrospective audit

## Headline result

计算结果支持保留现有 Combined 构建

The reconstructed mature circle is 2,013 nt (plasmid positions 723-2735;
junction 2735->723).  Its 1,083-nt CDS encodes the intended 360-aa protein
without internal stops.

| Metric | Current Combined | Empirical desirability percentile (all random synonymous candidates) |
|---|---:|---:|
| Circular MFE | -779.40 kcal/mol | 95.7% |
| Human CAI | 0.8110 | 87.5% |
| IRES deviation, 666-nt boundary | 11.9540 | 69.1% |

Pareto layer: 4; synonymous random candidates that
simultaneously dominate the current construct in all three objectives:
5 (1.95%).

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
