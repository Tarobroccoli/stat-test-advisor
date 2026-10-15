#!/usr/bin/env Rscript
# Descriptive diagnostics for a one-factor comparison (base R only).
# Usage: Rscript diagnose.R <file.csv|tsv> <response_col> <group_col> [outdir]
# Writes: <outdir>/diagnostics.txt, <outdir>/diagnostics.png
# Shapiro and Brown-Forsythe are reported as descriptive information, not as test-selection gates.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript diagnose.R <file> <response_col> <group_col> [outdir]")
file <- args[1]; resp <- args[2]; grp <- args[3]
outdir <- if (length(args) >= 4) args[4] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sep <- if (grepl("\\.tsv$|\\.txt$", file)) "\t" else ","
d <- read.table(file, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
for (col in c(resp, grp)) if (!col %in% names(d)) stop("Column not found: ", col)
d <- d[!is.na(d[[resp]]) & !is.na(d[[grp]]), ]
y <- as.numeric(d[[resp]]); g <- factor(d[[grp]])

out <- character()
say <- function(...) out <<- c(out, sprintf(...))

# Group summary
say("Response: %s | Group: %s | N = %d | groups = %d", resp, grp, length(y), nlevels(g))
say("%-20s %4s %10s %10s %10s %10s %10s", "group", "n", "mean", "sd", "median", "min", "max")
for (lv in levels(g)) {
  v <- y[g == lv]
  say("%-20s %4d %10.4g %10.4g %10.4g %10.4g %10.4g", lv, length(v), mean(v), sd(v), median(v), min(v), max(v))
}
sds <- tapply(y, g, sd)
say("SD ratio (max/min): %.2f", max(sds, na.rm = TRUE) / min(sds, na.rm = TRUE))
ns <- table(g)
say("Smallest group n: %d%s", min(ns), if (min(ns) <= 5) "  (too small to judge normality from data)" else "")
if (all(y > 0)) say("All values > 0: log-transform is an option if right-skewed (skewness = %.2f)",
                    mean((y - mean(y))^3) / sd(y)^3)

# Residuals from one-way model
fit <- lm(y ~ g)
r <- residuals(fit)
if (length(r) >= 3 && length(r) <= 5000) {
  sw <- shapiro.test(r)
  say("Shapiro-Wilk on residuals: W = %.3f, p = %.3g (descriptive)", sw$statistic, sw$p.value)
}
# Brown-Forsythe (Levene with median): ANOVA on |y - group median|
z <- abs(y - ave(y, g, FUN = median))
bf <- anova(lm(z ~ g))
say("Brown-Forsythe (Levene, median-centred): F = %.3f, p = %.3g (descriptive)", bf$`F value`[1], bf$`Pr(>F)`[1])

# Outliers: |standardized residual| > 3
sr <- rstandard(fit)
idx <- which(abs(sr) > 3)
say("Outliers (|std residual| > 3): %s", if (length(idx)) paste(sprintf("row %d (%s, %.4g)", idx, g[idx], y[idx]), collapse = "; ") else "none")

writeLines(out, file.path(outdir, "diagnostics.txt"))
cat(out, sep = "\n")

png(file.path(outdir, "diagnostics.png"), width = 1500, height = 500, res = 130)
par(mfrow = c(1, 3))
stripchart(y ~ g, vertical = TRUE, method = "jitter", pch = 19, col = "grey30",
           main = "Data by group", ylab = resp)
boxplot(y ~ g, add = TRUE, border = "steelblue", col = NA, outline = FALSE)
qqnorm(r, main = "Normal QQ of residuals", pch = 19); qqline(r, col = "red")
plot(fitted(fit), r, pch = 19, xlab = "Fitted (group mean)", ylab = "Residual",
     main = "Residuals vs fitted"); abline(h = 0, lty = 2)
invisible(dev.off())
cat("Plots:", file.path(outdir, "diagnostics.png"), "\n")
