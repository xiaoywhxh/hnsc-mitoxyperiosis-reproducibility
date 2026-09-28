## ============================================================
## P1-1: External validation random-effects meta-analysis
## GPT scholar 5.2 MC22 / MC21
## Descriptive/exploratory random-effects meta of the 67-gene
## program HR across 3 external OS cohorts (GSE41613, GSE42743,
## GSE65858). Reports pooled HR, I2, tau2, and prediction
## interval (PI more important than pooled P).
## ============================================================
suppressMessages({ library(stats) })

OUT_DIR <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

## ---------- 1. Data: per-cohort 67-gene program HR (log HR + SE) ----------
## From P0A_program_external / program_external_67gene.csv
cohorts <- data.frame(
  Cohort = c("GSE41613", "GSE42743", "GSE65858"),
  N      = c(97, 74, 270),
  Events = c(51, 42, 94),
  HR     = c(1.203950, 0.960601, 0.917111),
  LCL    = c(0.905881, 0.689538, 0.751790),
  UCL    = c(1.600095, 1.338220, 1.118787),
  P      = c(0.200939, 0.812169, 0.393553),
  stringsAsFactors = FALSE
)
cohorts$logHR <- log(cohorts$HR)
cohorts$SE <- (log(cohorts$UCL) - log(cohorts$LCL)) / (2 * qnorm(0.975))
cat("Per-cohort logHR / SE:\n")
print(cohorts[, c("Cohort", "N", "Events", "logHR", "SE")], row.names = FALSE)

## ---------- 2. DerSimonian-Laird random-effects meta ----------
k <- nrow(cohorts)
w_fixed <- 1 / cohorts$SE^2
## Fixed-effect pooled
b_fixed <- sum(w_fixed * cohorts$logHR) / sum(w_fixed)
Q <- sum(w_fixed * (cohorts$logHR - b_fixed)^2)
df_Q <- k - 1
tau2 <- max(0, (Q - df_Q) / (sum(w_fixed) - sum(w_fixed^2) / sum(w_fixed)))
## Random-effects weights
w_re <- 1 / (cohorts$SE^2 + tau2)
b_re <- sum(w_re * cohorts$logHR) / sum(w_re)
se_re <- sqrt(1 / sum(w_re))
## 95% CI for pooled
z <- qnorm(0.975)
ci_low <- b_re - z * se_re; ci_high <- b_re + z * se_re
## Prediction interval (random effects, k-2 df t-distribution)
se_pred <- sqrt(se_re^2 + tau2)
t_crit <- qt(0.975, df = k - 2)
pi_low <- b_re - t_crit * se_pred; pi_high <- b_re + t_crit * se_pred
## Heterogeneity
I2 <- max(0, (Q - df_Q) / Q) * 100
H2 <- max(1, Q / df_Q)
## Test of pooled effect (z-test)
z_pool <- b_re / se_re
p_pool <- 2 * pnorm(-abs(z_pool))

cat("\n===== Random-effects meta (DerSimonian-Laird) =====\n")
cat(sprintf("tau2 = %.4f (tau = %.3f)\n", tau2, sqrt(tau2)))
cat(sprintf("Q = %.2f (df=%d), I2 = %.1f%%, H2 = %.2f\n", Q, df_Q, I2, H2))
cat(sprintf("Pooled HR = %.3f (95%% CI %.3f-%.3f), P = %.4g\n",
            exp(b_re), exp(ci_low), exp(ci_high), p_pool))
cat(sprintf("95%% Prediction interval: %.3f-%.3f\n", exp(pi_low), exp(pi_high)))
cat(sprintf("PI crosses 1: %s\n", ifelse(pi_low < 0 & pi_high > 0, "YES", ifelse(exp(pi_low) < 1 & exp(pi_high) > 1, "YES", "NO"))))

## ---------- 3. Also include TCGA as sensitivity (4 cohorts) ----------
tcga_row <- data.frame(Cohort = "TCGA (discovery)", N = 501, Events = 218,
                       HR = 1.192, LCL = 1.034, UCL = 1.375, P = 0.015,
                       stringsAsFactors = FALSE)
tcga_row$logHR <- log(tcga_row$HR)
tcga_row$SE <- (log(tcga_row$UCL) - log(tcga_row$LCL)) / (2 * qnorm(0.975))
cohorts4 <- rbind(tcga_row, cohorts)

w_fixed4 <- 1 / cohorts4$SE^2
b_fixed4 <- sum(w_fixed4 * cohorts4$logHR) / sum(w_fixed4)
Q4 <- sum(w_fixed4 * (cohorts4$logHR - b_fixed4)^2)
df_Q4 <- nrow(cohorts4) - 1
tau2_4 <- max(0, (Q4 - df_Q4) / (sum(w_fixed4) - sum(w_fixed4^2) / sum(w_fixed4)))
w_re4 <- 1 / (cohorts4$SE^2 + tau2_4)
b_re4 <- sum(w_re4 * cohorts4$logHR) / sum(w_re4)
se_re4 <- sqrt(1 / sum(w_re4))
ci4_low <- b_re4 - z * se_re4; ci4_high <- b_re4 + z * se_re4
se_pred4 <- sqrt(se_re4^2 + tau2_4)
t_crit4 <- qt(0.975, df = nrow(cohorts4) - 2)
pi4_low <- b_re4 - t_crit4 * se_pred4; pi4_high <- b_re4 + t_crit4 * se_pred4
I2_4 <- max(0, (Q4 - df_Q4) / Q4) * 100
z_pool4 <- b_re4 / se_re4; p_pool4 <- 2 * pnorm(-abs(z_pool4))

cat("\n===== Sensitivity: 4 cohorts incl TCGA =====\n")
cat(sprintf("tau2 = %.4f, I2 = %.1f%%\n", tau2_4, I2_4))
cat(sprintf("Pooled HR = %.3f (95%% CI %.3f-%.3f), P = %.4g\n",
            exp(b_re4), exp(ci4_low), exp(ci4_high), p_pool4))
cat(sprintf("95%% PI: %.3f-%.3f\n", exp(pi4_low), exp(pi4_high)))

## ---------- 4. Save ----------
meta_res <- data.frame(
  Analysis = c("External only (3 cohorts)", "External + TCGA (4 cohorts)"),
  k = c(k, nrow(cohorts4)),
  Pooled_HR = c(exp(b_re), exp(b_re4)),
  LCL = c(exp(ci_low), exp(ci4_low)), UCL = c(exp(ci_high), exp(ci4_high)),
  P = c(p_pool, p_pool4),
  I2_pct = c(I2, I2_4),
  tau2 = c(tau2, tau2_4),
  PI_low = c(exp(pi_low), exp(pi4_low)), PI_high = c(exp(pi_high), exp(pi4_high)),
  Q = c(Q, Q4), df = c(df_Q, df_Q4)
)
write.csv(meta_res, file.path(OUT_DIR, "P1-1_external_meta_67g.csv"), row.names = FALSE)
cat("\nSaved P1-1_external_meta_67g.csv\n")
print(meta_res, row.names = FALSE)

## ---------- 5. Forest plot ----------
png(file.path(OUT_DIR, "Fig_P1-1_external_meta_forest.png"), width = 2200, height = 1600, res = 300)
par(mar = c(5, 6, 3, 2))
plot_hr <- function(hr, lcl, ucl, labels, main) {
  xlim <- c(min(lcl, exp(pi_low), 0.6) * 0.9, max(ucl, exp(pi_high), 1.8) * 1.1)
  plot(NA, xlim = xlim, ylim = c(0.5, length(labels) + 1.5), xlab = "HR per 1 SD (95% CI)",
       ylab = "", main = main, yaxt = "n", log = "x")
  abline(v = 1, lty = 2, col = "gray50")
  for (i in seq_along(labels)) {
    y <- length(labels) - i + 1
    points(hr[i], y, pch = 15, cex = 1.3)
    segments(lcl[i], y, ucl[i], y, lwd = 2)
    text(xlim[1], y, labels = labels[i], adj = 1, xpd = TRUE, cex = 0.9)
    text(hr[i], y + 0.28, labels = sprintf("%.2f (%.2f-%.2f)", hr[i], lcl[i], ucl[i]),
         cex = 0.8, pos = 3)
  }
  ## Pooled diamond
  yp <- 0.8
  points(exp(b_re), yp, pch = 18, cex = 2, col = "red")
  segments(exp(ci_low), yp, exp(ci_high), yp, lwd = 4, col = "red")
  ## PI line
  segments(exp(pi_low), yp, exp(pi_high), yp, lwd = 1.5, col = "red", lty = 3)
  text(xlim[1], yp, labels = "RE pooled (diamond) / PI (dashed)", adj = 1, xpd = TRUE,
       cex = 0.8, col = "red")
  legend("topright", legend = c("Cohort", "Pooled", "Prediction interval"),
         pch = c(15, 18, NA), lty = c(NA, NA, 3), col = c("black", "red", "red"),
         cex = 0.8, bty = "n")
}
plot_hr(cohorts$HR, cohorts$LCL, cohorts$UCL, cohorts$Cohort,
        "External cohorts (67-gene program) — random-effects meta")
dev.off()
cat("Saved Fig_P1-1_external_meta_forest.png\n")

cat("\n=== P1-1 external meta DONE ===")
