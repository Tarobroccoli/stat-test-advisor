# Small-scale / phenotypic / experimental data

Contents: 1. Unit of replication · 2. Continuous outcome · 3. Counts · 4. Proportions/binary ·
5. Ordinal scores · 6. Correlation/regression · 7. Post hoc choice · 8. Small-n facts · 9. R snippets

## 1. Unit of replication (check first)

The replicate is the unit that independently received the treatment. Technical replicates,
several leaves from one plant, several plants in one pot/tray treated together are **not**
independent. Either average them to one value per true replicate, or model them as a random
effect (`lmer(y ~ trt + (1|pot))`). Treating them as n inflates significance (pseudoreplication).

## 2. Continuous outcome

| Design | Default | When diagnostics are poor |
|---|---|---|
| 2 independent groups | **Welch t-test** (`t.test(var.equal = FALSE)`) | Right-skewed positive data → log-transform then Welch (ratio of geometric means). Outliers / ordinal-like → Mann-Whitney U (`wilcox.test`) — tests P(X>Y), not medians unless shapes equal. Permutation test is a good option. |
| 2 paired (before/after, same plant, matched pairs) | **Paired t-test** (normality of the *differences*) | Wilcoxon signed-rank (`wilcox.test(paired = TRUE)`) |
| ≥3 independent groups, 1 factor | **Welch ANOVA** (`oneway.test`) → **Games-Howell**; if SDs similar, classic ANOVA (`aov`) → **Tukey HSD** | Kruskal-Wallis → **Dunn** test (Holm or BH) |
| ≥3 groups vs one control only | ANOVA → **Dunnett** (`emmeans(..., "trt.vs.ctrl")`) | Kruskal-Wallis → Dunn vs control only (Holm) |
| 2+ factors (e.g. genotype × treatment) | `lm(y ~ A * B)`; test interaction first; type II (no interaction) or type III with sum contrasts; post hoc via **emmeans** (simple effects within levels if interaction) | log-transform; aligned rank transform (`ARTool`). Don't run separate Kruskal tests per factor — that ignores the interaction. |
| Blocks / batches / fields | Include block: `lm(y ~ trt + block)` or `lmer(y ~ trt + (1|block))` | — |
| Repeated measures over time / nested sampling | **Linear mixed model** `lmer(y ~ trt * time + (1|plant))` with lmerTest (Satterthwaite / Kenward-Roger df) → emmeans | Friedman test only for simple complete RM with one factor |
| Covariate (initial size, age) | ANCOVA `lm(y ~ covariate + trt)`; check slopes similar | — |

Heteroscedasticity in models beyond one factor: `nlme::gls(..., weights = varIdent(form = ~1|group))`
or log-transform.

## 3. Count outcome (lesions, colonies, insects, seeds per plant)

Don't ANOVA raw counts when means are small. Use a GLM: Poisson `glm(family = poisson)`;
check overdispersion (residual deviance / df ≫ 1) → negative binomial (`MASS::glm.nb`,
`glmmTMB`); many zeros → zero-inflated model. Post hoc with emmeans (`type = "response"` for
ratios).

## 4. Proportion / binary outcome

- 2×2 or r×c table: chi-square; if any expected count < 5 → **Fisher's exact test**.
- Paired binary: McNemar.
- x successes out of n (germination 18/25, survival): **binomial GLM**
  `glm(cbind(x, n - x) ~ trt, family = binomial)`, not ANOVA on percentages. Overdispersed →
  quasibinomial or `glmer` with observation-level random effect.
- Few dishes/pots per group (≤ ~5): a `glmer(... + (1|dish))` is often singular (dish variance
  = 0) and then silently treats every seed as independent — too liberal. Prefer quasibinomial
  and test at the dish level: `emmeans(fit, ..., df = df.residual(fit))` (emmeans uses z / df =
  Inf for glm by default). Don't pick the seed-level model because its p-values are smaller.

## 5. Ordinal score (disease severity 0–5, visual rating)

Cumulative link model `ordinal::clm(score ~ trt)` (or `clmm` with random effects). Simple
alternative: Mann-Whitney / Kruskal-Wallis + Dunn. Don't treat a 0–5 scale as continuous with
t-tests unless averaged over many sub-samples.

## 6. Correlation / association

Pearson (linear, roughly bivariate normal) vs Spearman (monotonic, ranks, robust to outliers).
Many traits pairwise → adjust with BH. Correlation ≠ agreement (use Bland-Altman for methods).

## 7. Post hoc choice

- An omnibus test isn't required before **planned** contrasts; test those directly with Holm.
- All pairwise after ANOVA with equal-ish SDs → Tukey HSD.
- Unequal SDs / sizes → Games-Howell (`rstatix::games_howell_test` or `PMCMRplus::gamesHowellTest`).
- vs control only → Dunnett.
- After Kruskal-Wallis → Dunn (`FSA::dunnTest` or `rstatix::dunn_test`) with Holm (few
  comparisons) or BH (many).
- After models (lm/glm/lmer) → `emmeans` + `pairs()` (adjust = "tukey") or `contrast()` with
  "holm" / "mvt".
- Compact letter display (`multcomp::cld`) is fine for figures; state it hides effect sizes.

## 8. Small-n facts worth telling the user

- n = 3 vs 3: the smallest possible two-sided Mann-Whitney p is 0.1 → it can never be
  significant at 0.05. n = 4 vs 4: minimum ≈ 0.029. Rank tests are not "safer" at tiny n.
- With n ≤ 5 per group normality can't be checked from the data; rely on knowledge of the
  variable (log-transform concentrations, expression ratios, sizes that grow multiplicatively).
- Report the effect size with CI (mean difference, ratio, Hedges' g, rank-biserial r); a
  non-significant result at n = 3 is "inconclusive", not "no effect".

## 9. R snippets (base + emmeans/lme4 where possible)

```r
# Welch t-test + mean difference CI
t.test(y ~ group, data = d)              # var.equal = FALSE by default

# Welch ANOVA, then Games-Howell (rstatix) or pairwise Welch t-tests with Holm (base fallback)
oneway.test(y ~ group, data = d)
pairwise.t.test(d$y, d$group, pool.sd = FALSE, p.adjust.method = "holm")

# Classic ANOVA + Tukey
fit <- aov(y ~ group, data = d); summary(fit); TukeyHSD(fit)

# Kruskal + Dunn (base fallback: pairwise Wilcoxon with Holm)
kruskal.test(y ~ group, data = d)
pairwise.wilcox.test(d$y, d$group, p.adjust.method = "holm")

# Factorial + emmeans
library(emmeans)
fit <- lm(y ~ genotype * treatment, data = d)
anova(fit)                                   # sequential (type I); use car::Anova(type = 2/3) if unbalanced
emmeans(fit, pairwise ~ treatment | genotype, adjust = "tukey")

# Vs control (Dunnett-type)
emmeans(fit, trt.vs.ctrl ~ treatment, ref = 1)

# Mixed model
library(lme4)                                # lmerTest for p-values with Satterthwaite df
m <- lmer(y ~ treatment * time + (1 | plant), data = d)
emmeans(m, pairwise ~ treatment | time)

# Counts / proportions
glm(count ~ trt, family = poisson, data = d)              # check deviance/df
glm(cbind(germ, total - germ) ~ trt, family = binomial, data = d)
```
