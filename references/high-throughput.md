# High-throughput data

Core rule: don't run Shapiro/Levene per feature and then pick t-test vs Wilcoxon feature by
feature. Use the framework built for the data type — it models the distribution (counts,
compositions, mean-variance trend) and borrows information across features, which matters a
lot at small n. Always include known batch/block factors in the design formula.

Contents: RNA-seq · single-cell · log-intensity omics · microbiome · GWAS/SNP · enrichment ·
many phenotypic traits

## Bulk RNA-seq (raw counts)

- **DESeq2** (`~ batch + condition`), **edgeR** quasi-likelihood (`glmQLFit` / `glmQLFTest`,
  good at small n), or **limma-voom** (fast, flexible designs, many samples).
- Input raw integer counts — never TPM/FPKM for DE.
- Correction: BH (`padj`). DESeq2's independent filtering and optional IHW are fine.
- Report log2 fold change; use shrinkage (`lfcShrink(type = "apeglm")`) for ranking/plots.
- ≥3 biological replicates per group; at 2 the results are weak — say so.

## Single-cell RNA-seq, differential expression between conditions

Cells are not replicates. Aggregate to **pseudobulk** per sample × cell type, then DESeq2/edgeR.
Per-cell Wilcoxon (Seurat default) is for marker finding, not condition comparisons.

## Microarray, proteomics, metabolomics (continuous intensities)

- log2 transform, normalize (quantile / median / VSN as appropriate), then **limma** moderated
  t / F (`lmFit` → `eBayes`), `topTable(adjust = "BH")`.
- Proteomics missing values are often missing-not-at-random (below detection): filter features
  with too many NAs, impute deliberately (e.g. left-censored for MNAR) and say how.
- Metabolomics: PCA/PLS-DA are exploratory; significance claims come from univariate limma/
  t-tests + BH (PLS-DA VIP is not a test).

## Microbiome (16S / shotgun)

- **Alpha diversity** (Shannon, observed, Faith PD): one value per sample → small-scale route
  (Welch/ANOVA on Shannon or Kruskal-Wallis + Dunn; lm/lmer if covariates/blocks).
- **Beta diversity**: PERMANOVA (`vegan::adonis2(dist ~ group, permutations = 999)`), and
  **always** check dispersion with `betadisper` + `permutest` — different spread can make
  PERMANOVA significant without a location shift, especially with unbalanced groups. Pairwise
  PERMANOVA → BH/Holm. Few samples (e.g. 2 per group): permutations are too few for a
  meaningful p-value — report it as descriptive.
- **Differential abundance**: data are compositional; methods disagree (Nearing et al. 2022).
  Recommend **ANCOM-BC2** (or ANCOM-BC in QIIME 2) plus one other — **ALDEx2**, **MaAsLin2/3**
  or **LinDA** — and highlight taxa found by both. Don't rarefy for DA; filter rare taxa
  (prevalence/abundance) before testing. BH correction (built into these tools).

## GWAS / SNP association

- Mixed model / multi-locus model controlling structure and kinship: GEMMA, GAPIT (MLM,
  FarmCPU, BLINK), EMMAX, or PLINK2 with PCs as covariates.
- Threshold: Bonferroni on the **effective** number of independent tests (LD-pruned SNP count,
  simpleM / GEC M_eff) is common in plants; human convention 5 × 10⁻⁸. FDR (BH) is sometimes
  used as a suggestive threshold — state which.
- Check QQ plot and genomic inflation λ (≈1 is good; >1.1 suggests uncontrolled structure).

## Enrichment (GO / KEGG)

- Over-representation: hypergeometric/Fisher with the **right background** (all genes that were
  testable/expressed, not the whole genome), BH across terms.
- Ranked list: GSEA (`fgsea`, `clusterProfiler::GSEA`), BH on gene sets.

## Many phenotypic traits / features with no dedicated tool

Run the small-scale model per feature (lm/Welch/Wilcoxon chosen once for all features by design
and data type, not per feature), then BH across features. Per-feature Shapiro is meaningless
here.
