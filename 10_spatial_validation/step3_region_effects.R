# Region (author region_4class) effect sizes, computed PER SAMPLE (no pooled spot-level p).
# Primary: core vs nc ; Secondary: edge vs nc, transitory vs nc.
# Per sample: median difference, Cliff's delta (rank-biserial), and a per-sample
# Mann-Whitney U as a *within-sample descriptive* statistic (NOT used for cross-sample
# inference). Cross-sample summary = sign test on direction + random-effects meta on delta.
suppressPackageStartupMessages({ library(data.table); library(metafor); library(jsonlite) })
IN  <- "D:/GSE208253_spatial/spatial_analysis/input"
SC  <- "D:/GSE208253_spatial/spatial_analysis/score"
OUT <- "D:/GSE208253_spatial/spatial_analysis/effects"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

spots  <- fread(file.path(IN, "spot_table_all.tsv"), sep = "\t", header = TRUE)
scores <- fread(file.path(SC, "spot_scores.tsv"), sep = "\t", header = TRUE)
d <- merge(spots, scores, by = c("sample", "barcode"), all.x = TRUE)
d <- d[!is.na(region_4class) & region_4class != ""]

cliff_delta <- function(a, b) {
  # P(a>b) - P(a<b) via rank-based computation
  n1 <- length(a); n2 <- length(b)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(a, b))
  r1 <- sum(r[seq_len(n1)])
  U <- r1 - n1 * (n1 + 1) / 2          # U for a
  (2 * U) / (n1 * n2) - 1               # Cliff's delta
}

comparisons <- list(
  core_vs_nc       = c("core", "nc"),
  edge_vs_nc       = c("edge", "nc"),
  transitory_vs_nc = c("transitory", "nc")
)

rows <- list()
for (cmp in names(comparisons)) {
  A <- comparisons[[cmp]][1]; B <- comparisons[[cmp]][2]
  for (s in paste0("s", 1:12)) {
    da <- d[sample == s & region_4class == A]
    db <- d[sample == s & region_4class == B]
    if (nrow(da) < 3 || nrow(db) < 3) {
      rows[[paste(cmp, s)]] <- data.table(comparison = cmp, sample = s,
        n_A = nrow(da), n_B = nrow(db),
        median_A = NA_real_, median_B = NA_real_, median_diff = NA_real_,
        cliffs_delta = NA_real_, MWU_p = NA_real_)
      next
    }
    md <- median(da$ucell) - median(db$ucell)
    cd <- cliff_delta(da$ucell, db$ucell)
    wt <- suppressWarnings(wilcox.test(da$ucell, db$ucell, exact = FALSE))
    rows[[paste(cmp, s)]] <- data.table(
      comparison = cmp, sample = s, n_A = nrow(da), n_B = nrow(db),
      median_A = median(da$ucell), median_B = median(db$ucell),
      median_diff = md, cliffs_delta = cd, MWU_p = wt$p.value)
  }
}
eff <- rbindlist(rows)
fwrite(eff, file.path(OUT, "region_effect_sizes.tsv"), sep = "\t")

# --- cross-sample summary per comparison ---
summ <- list()
for (cmp in names(comparisons)) {
  sub <- eff[comparison == cmp & !is.na(cliffs_delta)]
  n <- nrow(sub)
  n_pos <- sum(sub$cliffs_delta > 0)
  n_neg <- sum(sub$cliffs_delta < 0)
  st <- binom.test(n_pos, n, 0.5)      # sign test vs 50/50
  # random-effects meta-analysis on Cliff's delta (with sampling SE approximation)
  # SE(delta) approx = sqrt((n1+n2+1)/(3*n1*n2))  (Hanley-McNeil style for delta)
  sub[, se := sqrt((n_A + n_B + 1) / (3 * n_A * n_B))]
  m <- tryCatch(rma(yi = cliffs_delta, sei = se, data = as.data.frame(sub),
                    method = "REML"), error = function(e) NULL)
  summ[[cmp]] <- data.table(
    comparison = cmp, n_samples = n,
    n_delta_positive = n_pos, n_delta_negative = n_neg,
    sign_test_p = st$p.value,
    median_of_median_diff = median(sub$median_diff),
    mean_delta = mean(sub$cliffs_delta),
    meta_delta = if (!is.null(m)) as.numeric(m$beta) else NA_real_,
    meta_se = if (!is.null(m)) m$se else NA_real_,
    meta_ci_lo = if (!is.null(m)) m$ci.lb else NA_real_,
    meta_ci_hi = if (!is.null(m)) m$ci.ub else NA_real_,
    meta_p = if (!is.null(m)) m$pval else NA_real_,
    meta_tau2 = if (!is.null(m)) m$tau2 else NA_real_,
    meta_I2 = if (!is.null(m)) m$I2 else NA_real_)
}
ss <- rbindlist(summ)
fwrite(ss, file.path(OUT, "region_cross_sample_summary.tsv"), sep = "\t")

# --- leave-sample_3-out sensitivity (drop s3) ---
summ_l3 <- list()
for (cmp in names(comparisons)) {
  sub <- eff[comparison == cmp & !is.na(cliffs_delta) & sample != "s3"]
  n <- nrow(sub); n_pos <- sum(sub$cliffs_delta > 0)
  st <- binom.test(n_pos, n, 0.5)
  sub[, se := sqrt((n_A + n_B + 1) / (3 * n_A * n_B))]
  m <- tryCatch(rma(yi = cliffs_delta, sei = se, data = as.data.frame(sub),
                    method = "REML"), error = function(e) NULL)
  summ_l3[[cmp]] <- data.table(
    comparison = cmp, n_samples = n, n_delta_positive = n_pos,
    sign_test_p = st$p.value, meta_delta = if (!is.null(m)) as.numeric(m$beta) else NA_real_,
    meta_p = if (!is.null(m)) m$pval else NA_real_)
}
fwrite(rbindlist(summ_l3), file.path(OUT, "region_leave_s3_out_summary.tsv"), sep = "\t")

cat("=== REGION EFFECT SIZES (UCell) ===\n")
print(eff[order(comparison, sample), .(comparison, sample, n_A, n_B,
        median_A = round(median_A,3), median_B = round(median_B,3),
        median_diff = round(median_diff,3), cliffs_delta = round(cliffs_delta,3))])
cat("\n=== CROSS-SAMPLE SUMMARY ===\n")
print(ss)
cat("\n=== LEAVE-s3-OUT SUMMARY ===\n")
print(rbindlist(summ_l3))
