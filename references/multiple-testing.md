# Multiple testing — which correction

First define the **family**: the set of tests whose results will be read together as one
claim (all pairwise comparisons after one ANOVA; all genes for one contrast; all GO terms).
Correct within the family. Separate contrasts in DE analyses are usually corrected separately.

## Control family-wise error (any false positive is costly; few tests)

| Method | Use when | R |
|---|---|---|
| **Holm** | Default FWER choice; always at least as powerful as Bonferroni, valid under any dependence | `p.adjust(p, "holm")` |
| Bonferroni | Simple thresholds (GWAS α / M_eff); otherwise prefer Holm | `p.adjust(p, "bonferroni")` |
| Tukey HSD | All pairwise means after ANOVA | `TukeyHSD`, `emmeans(adjust = "tukey")` |
| Dunnett | Each treatment vs one control | `emmeans(trt.vs.ctrl)` |
| Games-Howell | All pairwise, unequal variances | `rstatix::games_howell_test` |
| Hochberg / Hommel | Slightly more power than Holm; need independence or positive dependence | `p.adjust(p, "hochberg")` |

Typical: post hoc in a small experiment, a handful of planned contrasts, confirmatory tests.

## Control false discovery rate (many tests; discovery/screening)

| Method | Use when | R |
|---|---|---|
| **Benjamini-Hochberg (BH)** | Default for omics; valid for independent or positively dependent tests (PRDS), which covers most gene/taxon tests | `p.adjust(p, "BH")` |
| Benjamini-Yekutieli (BY) | Arbitrary/negative dependence must be guaranteed; very conservative | `p.adjust(p, "BY")` |
| Storey q-value | Large m with many true effects (estimates π₀ → more power than BH); needs a well-behaved p-value histogram | `qvalue::qvalue(p)$qvalues` |
| IHW | Have an informative covariate independent of p under H0 (e.g. mean count) | `IHW::ihw()`; DESeq2 `filterFun = ihw` |

FDR 0.05 is common; 0.1 is accepted for exploratory omics (DESeq2 default `alpha = 0.1`).
State the threshold used.

## Practical rules

- The terms of one designed ANOVA/factorial model (e.g. N, P, K main effects and their
  interactions) are conventionally **not** adjusted for multiplicity — each answers its own
  planned question. Adjust the post hoc comparisons *within* a term instead. If you choose to
  adjust anyway, report the unadjusted p too and say why.
- Check the p-value histogram for high-throughput results: uniform + spike near 0 is healthy;
  a hump near 1 or U-shape suggests a mis-specified model.
- Don't mix: report "BH-adjusted p (padj)" or "q-value" explicitly — they're not identical.
- Pairwise tests after Kruskal-Wallis: Dunn + Holm (few groups) or BH (many).
- Many traits × few contrasts: BH over all trait × contrast tests if they form one claim;
  otherwise per contrast — say which.
- Effect-size filter (|log2FC| > 1) after FDR is fine for display, but for formal tests use
  `lfcThreshold` (DESeq2) / `treat` (limma) instead of post-filtering.
