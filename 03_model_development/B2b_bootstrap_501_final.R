## B2b: Bootstrap internal validation — FINAL 501-case 11-gene model (审稿意见 2 P0-1)
## 依据: Hm 审稿意见2 GPT.docx — "0.505 来自旧 226 例流程，必须在最终 501 例 canonical 73-gene 流程上重跑"
## 每次 bootstrap 完整重复: resampling → gene standardization → CV-LASSO → lambda selection
##                          → feature selection → coefficient estimation → boot apparent C → original test C
## 输出: apparent C / optimism / corrected C / distribution / 95% CI / selection frequency / failed iterations

suppressMessages({
  library(glmnet)
  library(survival)
  library(timeROC)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/05_model_development"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(20260720)
B <- 500

# ============================================================
# 1. Load FINAL 501-case full cohort + clinical
# ============================================================
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))  # 60660 x 502
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))  # 528 patients

clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
patient_ids <- substr(colnames(tumor_mat), 1, 12)
cm <- match(patient_ids, clin$bcr_patient_barcode)
surv_full <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm],
                        OS_status = clin$OS_status[cm], stringsAsFactors = FALSE)
surv501 <- surv_full[!is.na(surv_full$OS_time), ]
cat("FINAL cohort:", nrow(surv501), "samples,", sum(surv501$OS_status), "events\n")

# Canonical 73 genes (manuscript Methods 2.2)
core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1",
                "BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR",
                "TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB",
                "PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5",
                "PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2",
                 "CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5",
                 "STING1","CGAS","MAOB","ENDOG")
cand_genes <- c(core_genes, mtor_genes, mito_genes, metab_genes)
stopifnot(length(unique(cand_genes)) == 73)
cat("Candidate genes:", length(cand_genes), "\n")

# Expression: log2(TPM+1), gene-level Z (t(scale(t())))
expr <- log2(tumor_mat[intersect(cand_genes, rownames(tumor_mat)), surv501$sample] + 1)
stopifnot(nrow(expr) == 73)
X_raw <- t(expr)  # samples x genes
y <- surv501
X <- scale(X_raw)
attr(X, "scaled:center") <- NULL; attr(X, "scaled:scale") <- NULL
cat("Standardized X:", dim(X), "\n")

# ============================================================
# 2. Full-data LASSO (must reproduce locked_model_TCGA501.csv: 11 genes)
# ============================================================
fit_full <- cv.glmnet(X, Surv(y$OS_time, y$OS_status), family = "cox",
                      alpha = 1, nfolds = 10, foldid = rep(1:10, length.out = nrow(X)))
coef_full <- coef(fit_full, s = "lambda.min")[, 1]
sel_full <- names(coef_full)[coef_full != 0]
cat("Full-data selected genes:", length(sel_full), ":", paste(sel_full, collapse = ", "), "\n")

risk_full <- as.numeric(X[, sel_full, drop = FALSE] %*% coef_full[sel_full])
cfit_full <- coxph(Surv(y$OS_time, y$OS_status) ~ risk_full)
c_app <- summary(cfit_full)$concordance[1]
cat(sprintf("Apparent C-index: %.4f\n", c_app))

auc_app <- tryCatch({
  roc_full <- timeROC(T = y$OS_time, delta = y$OS_status, marker = risk_full,
                      cause = 1, times = c(365, 1095, 1825), iid = FALSE)
  c(roc_full$AUC[1], roc_full$AUC[2], roc_full$AUC[3])
}, error = function(e) c(NA, NA, NA))
cat(sprintf("Apparent AUC 1/3/5yr: %.3f / %.3f / %.3f\n", auc_app[1], auc_app[2], auc_app[3]))

# ============================================================
# 3. Bootstrap B=500 (full pipeline repetition)
# ============================================================
cat("\nStarting bootstrap:", B, "iterations...\n")
t0 <- Sys.time()

boot_cindex <- numeric(B)       # apparent C within bootstrap sample
boot_cindex_test <- numeric(B)  # boot-fit applied to original data (test C)
sel_freq <- setNames(numeric(ncol(X)), colnames(X))
coef_mat <- matrix(NA, nrow = B, ncol = ncol(X), dimnames = list(NULL, colnames(X)))
boot_nsel <- integer(B)         # number of selected genes per iteration
failed <- 0L

for (b in 1:B) {
  if (b %% 50 == 0) cat(sprintf("  boot %d/%d (%.1f min)\n", b, B, as.numeric(difftime(Sys.time(), t0, units="mins"))))
  idx <- sample(nrow(X), nrow(X), replace = TRUE)
  Xb <- X[idx, ]; yb <- y[idx, ]

  fit_b <- tryCatch(
    cv.glmnet(Xb, Surv(yb$OS_time, yb$OS_status), family = "cox",
              alpha = 1, nfolds = 10, foldid = rep(1:10, length.out = nrow(Xb))),
    error = function(e) NULL)
  if (is.null(fit_b)) { failed <- failed + 1L; next }
  coef_b <- coef(fit_b, s = "lambda.min")[, 1]
  sel_b <- names(coef_b)[coef_b != 0]
  boot_nsel[b] <- length(sel_b)
  coef_mat[b, ] <- coef_b
  if (length(sel_b) > 0) sel_freq[sel_b] <- sel_freq[sel_b] + 1

  if (length(sel_b) > 0) {
    rb <- as.numeric(Xb[, sel_b, drop = FALSE] %*% coef_b[sel_b])
    boot_cindex[b] <- tryCatch(summary(coxph(Surv(yb$OS_time, yb$OS_status) ~ rb))$concordance[1],
                               error = function(e) NA)
    rt <- as.numeric(X[, sel_b, drop = FALSE] %*% coef_b[sel_b])
    boot_cindex_test[b] <- tryCatch(summary(coxph(Surv(y$OS_time, y$OS_status) ~ rt))$concordance[1],
                                    error = function(e) NA)
  } else {
    failed <- failed + 1L  # no gene selected → uninformative iteration
  }
}

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\nBootstrap complete: %.1f min, failed/uninformative: %d\n", elapsed, failed))

# ============================================================
# 4. Optimism correction + full distribution reporting
# ============================================================
valid <- !is.na(boot_cindex) & !is.na(boot_cindex_test)
n_valid <- sum(valid)
cat(sprintf("Valid bootstrap iterations: %d / %d\n", n_valid, B))

optimism_vals <- boot_cindex[valid] - boot_cindex_test[valid]
optimism <- mean(optimism_vals, na.rm = TRUE)
optimism_sd <- sd(optimism_vals, na.rm = TRUE)
c_corrected <- c_app - optimism

ci_test <- quantile(boot_cindex_test[valid], c(0.025, 0.5, 0.975), na.rm = TRUE)
ci_opt <- quantile(optimism_vals, c(0.025, 0.5, 0.975), na.rm = TRUE)
ci_boot_app <- quantile(boot_cindex[valid], c(0.025, 0.5, 0.975), na.rm = TRUE)

cat(sprintf("\n=== Internal validation (FINAL 501-case 11-gene model) ===\n"))
cat(sprintf("Apparent C-index:            %.4f\n", c_app))
cat(sprintf("Optimism mean (SD):          %.4f (%.4f)\n", optimism, optimism_sd))
cat(sprintf("Optimism-corrected C-index:  %.4f\n", c_corrected))
cat(sprintf("Boot apparent C 95%% CI:      %.4f - %.4f (median %.4f)\n", ci_boot_app[1], ci_boot_app[3], ci_boot_app[2]))
cat(sprintf("Boot test C 95%% CI:          %.4f - %.4f (median %.4f)\n", ci_test[1], ci_test[3], ci_test[2]))
cat(sprintf("Optimism 95%% CI:             %.4f - %.4f (median %.4f)\n", ci_opt[1], ci_opt[3], ci_opt[2]))

# ============================================================
# 5. Outputs
# ============================================================
freq_df <- data.frame(Gene = names(sel_freq), Selection_Frequency = sel_freq / B)
freq_df <- freq_df[order(-freq_df$Selection_Frequency), ]
write.csv(freq_df, file.path(OUT_DIR, "bootstrap_selection_frequency_501.csv"), row.names = FALSE)
cat("\nTop selection frequency:\n")
print(head(freq_df[freq_df$Selection_Frequency > 0, ], 15))

coef_summary <- data.frame(
  Gene = colnames(X),
  Mean = colMeans(coef_mat, na.rm = TRUE),
  SD = apply(coef_mat, 2, sd, na.rm = TRUE),
  Pct_NonZero = apply(coef_mat, 2, function(x) mean(x != 0, na.rm = TRUE))
)
coef_summary <- coef_summary[order(-coef_summary$Pct_NonZero), ]
write.csv(coef_summary, file.path(OUT_DIR, "bootstrap_coefficients_501.csv"), row.names = FALSE)

perf_df <- data.frame(
  Metric = c("Apparent C-index", "Optimism_mean", "Optimism_SD",
             "Optimism-corrected C-index",
             "Boot_apparent_C_2.5pct", "Boot_apparent_C_median", "Boot_apparent_C_97.5pct",
             "Boot_test_C_2.5pct", "Boot_test_C_median", "Boot_test_C_97.5pct",
             "Optimism_2.5pct", "Optimism_median", "Optimism_97.5pct",
             "Apparent_AUC_1yr", "Apparent_AUC_3yr", "Apparent_AUC_5yr",
             "Valid_bootstraps", "Failed_uninformative"),
  Value = c(c_app, optimism, optimism_sd, c_corrected,
            ci_boot_app[1], ci_boot_app[2], ci_boot_app[3],
            ci_test[1], ci_test[2], ci_test[3],
            ci_opt[1], ci_opt[2], ci_opt[3],
            auc_app[1], auc_app[2], auc_app[3],
            n_valid, failed)
)
write.csv(perf_df, file.path(OUT_DIR, "internal_cindex_501_final.csv"), row.names = FALSE)
print(perf_df)

write.csv(data.frame(Gene = sel_full, Coefficient = unname(coef_full[sel_full])),
          file.path(OUT_DIR, "locked_model_coefficients_501_final.csv"), row.names = FALSE)

# 保存分布供图
saveRDS(list(app = boot_cindex[valid], test = boot_cindex_test[valid],
             optimism = optimism_vals, sel_freq = sel_freq,
             boot_nsel = boot_nsel[valid], n_valid = n_valid, B = B),
        file.path(OUT_DIR, "bootstrap_dist_501_final.Rds"))

cat("\n=== B2b DONE ===")
