# stat-test-advisor — Claude Code skill

A [Claude Code](https://claude.com/claude-code) skill that chooses, justifies, codes (R) and runs
the right hypothesis test for biological data — design first, not "Shapiro → pick a test".

- **Small-scale / phenotypic data:** Welch vs paired t, Mann-Whitney, Welch ANOVA / ANOVA /
  Kruskal-Wallis, factorial and blocked designs, mixed models, GLMs for counts and proportions,
  post hoc (Tukey, Games-Howell, Dunnett, Dunn, emmeans), small-n facts (e.g. Mann-Whitney 3 vs 3
  can never reach p < 0.05).
- **High-throughput data:** DESeq2 / edgeR / limma, pseudobulk scRNA-seq, microbiome
  (PERMANOVA + betadisper, ANCOM-BC2 / ALDEx2 / MaAsLin2), GWAS thresholds, enrichment.
- **Multiple testing:** when to use Holm, Tukey, Dunnett, BH, BY, Storey q-value, IHW.
- Shapiro / Levene are run and reported as descriptive diagnostics, not as test-selection gates.
- Output: design summary, recommended test with reasons, diagnostics, results with effect sizes
  and 95% CIs, and a ready-to-paste methods sentence.

## Install

```bash
git clone https://github.com/Tarobroccoli/stat-test-advisor.git ~/.claude/skills/stat-test-advisor
```

## Use

Ask Claude Code e.g. "which test should I use for this data?" with a CSV in the project, or type
`/stat-test-advisor`. `scripts/diagnose.R` (base R only) gives quick one-factor diagnostics:

```bash
Rscript scripts/diagnose.R data.csv <response_col> <group_col> [outdir]
```

## Evaluation (iteration 1)

7 small-scale cases (5 classic R datasets: sleep, PlantGrowth, InsectSprays, npk, CO2; plus a
generated qPCR n = 3 and a germination proportion dataset), with vs without the skill:
assertion pass rate 100% vs 89%. Main gains: dish-level inference for proportion data and a
methods sentence every time.
