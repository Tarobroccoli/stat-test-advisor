---
name: stat-test-advisor
description: Choose, justify, code (R) and run the right statistical hypothesis test for biological data — small phenotypic/experimental datasets (t-test vs Welch vs Mann-Whitney, paired tests, ANOVA/Welch ANOVA/Kruskal-Wallis, mixed models, GLMs for counts/proportions, post hoc Tukey/Games-Howell/Dunn/Dunnett/emmeans) and high-throughput data (RNA-seq, microbiome, proteomics/metabolomics, GWAS, enrichment) including which multiple-testing / FDR correction to use (BH, BY, q-value, Holm, Bonferroni). Use this whenever the user asks which test to use, mentions Shapiro/normality, Levene/equal variance, omnibus or post hoc tests, p-value adjustment, FDR, "is this significant", comparing groups/treatments, or wants to analyse a phenotype or omics table — even if they only paste data and ask "how should I analyse this".
---

# Statistical test advisor

Goal: pick a test that matches the **experimental design and the kind of data**, not one picked
by a chain of p-values. Then write R code, run it if data is available, and give a sentence the
user can paste into their methods.

## Why not "Shapiro → choose test"

The user may expect the classic flowchart (Shapiro p < 0.05 → non-parametric; Levene p < 0.05 →
…). Explain briefly, once, why this skill doesn't gate on those tests:
- Choosing a test from a pre-test changes the error rate of the final test.
- Shapiro has almost no power at n = 3–6 per group (passes skewed data) and rejects trivial
  deviations at large n (where t/ANOVA are robust anyway).
- Normality is an assumption about **residuals**, not raw values of each group.
- Welch t-test / Welch ANOVA handle unequal variances at almost no cost when variances are
  equal, so Levene doesn't need to decide anything.

So: the design + outcome type choose the family of test; diagnostics (residual QQ plot, group
SD ratio, outliers, knowledge of the variable — e.g. concentrations are usually log-normal)
decide between close alternatives. Shapiro (on residuals) and Brown-Forsythe/Levene are still
**run and reported**, because reviewers often ask for them — just labelled as descriptive.

## Workflow

1. **Intake.** Get these answers, from the data/script if possible, otherwise ask (all at once,
   short):
   - Outcome type: continuous / count / proportion or binary (x of n) / ordinal score /
     time-to-event — or high-throughput (many features: genes, ASVs, metabolites, SNPs).
   - Factors: how many groups, how many factors, is there a control to compare against?
   - Dependence: independent units, or paired/repeated measures, blocks, batches, nested
     sampling (several plants per pot, technical replicates, several leaves per plant)?
   - n per group (biological replicates — the unit that received the treatment).
   - Covariates? Hypothesis directional? Planned comparisons or all pairwise?
2. **Route.**
   - Small number of outcomes → read `references/small-scale.md`.
   - Many features → read `references/high-throughput.md`.
   - Any time more than one p-value is reported → read `references/multiple-testing.md`.
3. **Diagnose (if data given).** For a simple one-factor design run
   `Rscript <skill-dir>/scripts/diagnose.R <file.csv> <response_col> <group_col> [outdir]`
   (base R only; writes group summary, residual QQ/residual plots, Shapiro on residuals,
   Brown-Forsythe test, SD ratio). For other designs fit the planned model and check its
   residuals the same way. Look at the plots yourself before deciding.
4. **Check packages** with `Rscript -e 'requireNamespace("pkg")'` before using them. If one is
   missing, tell the user and give the install command (`install.packages()` or
   `BiocManager::install()`); offer a base-R fallback when one exists. Don't install silently.
5. **Write and run code.** Save the script in the project (e.g. `analysis/stats_<outcome>.R`)
   so it shows up later for the `pipeline-methods` skill. Report effect sizes with confidence
   intervals, not only p-values.
6. **Report** using the format below.

## Report format

```
Design: <outcome type>, <groups/factors>, <independent/paired/nested>, n = <per group>
Recommended: <omnibus test> → <post hoc / contrasts> → <correction>
Why: <2–4 bullets tied to the design and diagnostics>
Diagnostics: residual QQ <ok/heavy tails/skew>, SD ratio <x>, Shapiro (residuals) p = <>,
             Brown-Forsythe p = <> (descriptive only)
Alternative if <condition>: <test>
Results: <statistic, df, p, effect size [95% CI]; post hoc table>
Methods sentence: "<ready to paste, names test, software and version, correction>"
Caveats: <e.g. n = 3 per group: Mann-Whitney cannot reach p < 0.1>
```

Keep caveats honest: if n is too small for any test to say much, say so plainly rather than
hunting for a test that gives p < 0.05. Never switch tests after seeing p-values to get
significance; if the user asks for that, explain why it's p-hacking.
