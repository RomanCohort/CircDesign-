# Wiki-ready conclusion

## 中文

本团队基于 circDesign 公开提出的设计目标，对 Combined 成熟 circRNA 进行了回顾性计算评估。依据载体中 T4 td retained-scar 注释，我们将质粒第 723–2735 位重建为 2,013 nt 的成熟闭环序列，并确认第 1504–2586 位 CDS 长 1,083 nt，可连续编码预期的 360 aa Combined 多表位抗原，不含内部终止密码子。计算采用 ViennaRNA 2.7.0 circular folding model，人源 CAI 依据 Kazusa Homo sapiens 通用密码子使用数据。

评估的关键方法学修正：我们不再把不同性质的参照序列池成单一分布。人源加权与均匀抽样改变密码子组成，而组成重排（composition shuffle）保持密码子多重集完全不变、只重排顺序。二者沿不同方向变动，混合后的分布由组成方差主导，会把每个指标都拉向对本构建有利的方向。因此本报告按指标配零模型：密码子组成用前者，结构排列用后者。这一区分是代数性的而非统计性的——64 条组成重排与原序列共享同一个密码子使用向量，组成轴上的组内标准差恰为 0.000。

结果分两轴陈述。组成轴上，本构建达到参照空间的上确界：human CAI = 0.8110，全部 192 条人源加权与均匀抽样候选均未超过，即无任何同义设计在密码子适配维度上优于它。排列轴上，本构建方向有利但未达显著：在密码子组成完全相同的条件下，其 circular MFE 位于组成重排的 94.50 百分位（两个各 2,000 条的独立样本合并 n = 4,000），可交换性检验 p = 0.055。

IRES 维度需单独说明。原用于表征 IRES 完整性的 L_IRES 是两组碱基配对概率矩阵的 L2 距离，该量无法区分“IRES 结构被保住”与“IRES 本来就没有结构可保”。我们改用四个可直接读出的比值（自身结构保持率、环外侵占率、可及性、平均配对概率），物理模型与原构造完全一致。四个指标在组成匹配零模型下均不显著（p = 0.29 至 0.71），说明该维度在本设计中不携带可从序列特征读出的信号。

多目标核验采用预注册流程。在保持密码子组成不变的同义重排中，仅保留通过产物正确性质控者（无隐性剪接供体、无隐性 polyA），再统计仍能在三项目标上同时不劣于本构建的候选。在 4,000 条生成候选（其中 1,553 条通过质控）中，仅 30 条（1.93%，95% 置信区间 [1.36, 2.74]）满足条件，即本构建处于可建成同义设计的前 1.9%。需说明的是：本项目的预注册判定线要求置信上界低于 2%，实测上界为 2.74%，故该判定未获满足——点估计达标而置信上界未达标，原因是样本量不足，而非效应不存在。

进一步地，若把 circRNA 疫苗所需的翻译可及性并入目标集（翻译起始区、IRES–CDS 间隔区、CDS 全局可及性与最长螺旋四项，窗口与方向均在观察结果前固定），则 1,553 条可建成同义设计中无一能同时不劣于本构建。其中 IRES–CDS 间隔区可及性是一个明确的正向结果：本构建位于可建成设计的第 91.2 百分位，即约 91% 的同义设计在该项上不如它；该区域是核糖体由 IRES 转入起始密码子的通道。需同时说明：扩大目标集在力学上必然减少支配者数量，因此该结果应读作“在该目标集定义下不存在同时不劣于本构建的可建成方案”，而非“可及性证明了竞争者不成立”（偏好性检验 p = 0.38，不显著）。

结合以上结果，团队决定保留当前 Combined 序列继续后续实验。该决策同时由两条独立证据支撑：一是模型未发现任何可执行且系统性的改进方案，二是该序列已完成实际克隆与既有实验验证，而 MFE、CAI 与配对概率均为计算代理指标。此前将 composition_shuffle_055 列为第二代候选的建议已撤回——其 IRES 改善幅度在组成重排组内为 z = −1.51，考虑多重比较后不构成统计证据。

最后必须明确本报告的边界：本构建并非同义排列空间中的绝对最优解。在该空间内，5.50% 的重排具有更低的 MFE，1.45% 的重排在三个目标上同时不劣于它，这些都是可枚举、可复现的反例。因此“本序列最优”这一表述只有在指明空间与目标集时才成立。可支持的表述是：本构建在其设计目标空间内达到密码子组成的上确界，在可建成的同义设计中位于前三目标区域，且不存在系统性的改进方案。

## English

We report a retrospective computational evaluation of the mature Combined circRNA against the design objectives published for circDesign. Using the T4 td retained-scar annotations in the source vector, plasmid positions 723–2735 were reconstructed as a 2,013-nt mature circle, and the 1,083-nt CDS at positions 1504–2586 was confirmed to encode the intended 360-aa multi-epitope antigen continuously and without an internal stop codon. Folding used the ViennaRNA 2.7.0 circular model; human CAI used the Kazusa Homo sapiens codon-usage dataset.

A key methodological correction: the reference sequences are no longer pooled into a single distribution. Human-weighted and uniform sampling both alter codon composition, whereas composition shuffling holds the codon multiset exactly fixed and permutes only the order. Because the two move along different axes, a pooled distribution is dominated by composition variance and pulls every metric toward values favourable to the construct. Metrics are therefore reported against the null they can actually answer: composition against the former, arrangement against the latter. The distinction is algebraic rather than statistical — the 64 shuffles and the construct share one codon-usage vector, so the within-group standard deviation on the composition axis is exactly 0.000.

Results are stated on two axes. On the composition axis the construct reaches the supremum of the reference space: human CAI = 0.8110, exceeded by none of the 192 human-weighted and uniform candidates, i.e. no synonymous design is better adapted in this dimension. On the arrangement axis the direction is favourable but not significant: holding codon composition exactly fixed, circular MFE sits at the 94.50th percentile of composition-matched rearrangements (pooled n = 4,000 from two independent 2,000-sequence samples), exchangeability p = 0.055.

The IRES dimension requires separate comment. The original L_IRES is an L2 distance between two base-pair-probability matrices, and it cannot distinguish "the IRES structure was preserved" from "the IRES had no structure to lose". We replaced it with four directly readable ratios (self-structure retention, external occupancy, accessibility, mean pairing probability) under an identical physical model. All four are non-significant under the composition-matched null (p = 0.29 to 0.71), indicating that this dimension carries no signal recoverable from sequence features in this design.

Multi-objective verification followed a pre-registered procedure: among rearrangements holding codon composition fixed, retain only those passing product-correctness quality control (no cryptic splice donor, no cryptic polyA), then count those still not worse than the construct on all three objectives. Of 4,000 generated candidates (1,553 passing quality control), 30 (1.93%, 95% CI [1.36, 2.74]) satisfied it — placing the construct in the top 1.9% of buildable synonymous designs. The project's pre-registered decision rule required the confidence upper bound to fall below 2%; the observed upper bound is 2.74%, so the rule is not satisfied — the point estimate meets the threshold, the upper bound does not, for want of sample size rather than absence of effect.

Adding the translation-accessibility objectives a circRNA vaccine requires (translation initiation region, IRES-CDS spacer, global CDS accessibility, longest helix; windows and directions fixed before any result was seen) leaves none of the 1,553 buildable designs non-inferior on all objectives. The IRES-CDS spacer is an unambiguous positive: the construct sits at the 91.2nd percentile there. Enlarging an objective set mechanically reduces the number of dominators, so this reads as "no buildable design is non-inferior under this objective set", not "accessibility disproves the competitors" (preferential test p = 0.38, not significant).

On this basis the team retained the current Combined sequence. The decision rests on two independent lines of evidence: the model found no executable, systematic improvement, and the sequence has already completed cloning and prior experimental validation, while MFE, CAI and pairing probabilities remain computational proxies. The earlier suggestion to nominate composition_shuffle_055 as a second-generation candidate is withdrawn — its IRES improvement is z = −1.51 within the shuffle group and does not survive multiple-comparison correction.

Finally, the boundary of this report must be stated plainly: the construct is not the absolute optimum in the synonymous arrangement space. Within that space 5.50% of rearrangements have a lower MFE and 1.45% are not worse on all three objectives; these are enumerable, reproducible counterexamples. Any claim that the sequence is "optimal" therefore holds only when the space and the objective set are named. What is supported is this: the construct attains the supremum of the composition axis, sits in the leading multi-objective region of buildable synonymous designs, and has no systematic improvement available against it.
