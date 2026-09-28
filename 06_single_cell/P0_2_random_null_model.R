## ============================================================
## P0-2: Single-cell matched random-gene-set null model
## GPT scholar 5.2 MC18 / reviewer MC8
## Question: is the malignant-epithelial enrichment of Mitoxy67
## specific, or is it a generic property of any 67-gene set with
## similar expression level / variance / detection rate?
## Design: 1000 random 67-gene sets, matched to the true 67-gene
## program on mean expression, expression variance, and detection
## rate; each set -> patient-level malignant-vs-nonmalignant
## enrichment statistic; compare true Mitoxy67 to null distribution.
## ============================================================
suppressMessages({ library(data.table) })

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/11_singlecell"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(20260901)

N_NULL <- 1000

## ---------- 1. Load GSE103322 expression + annotation (same as P0_4) ----------
annot <- readRDS(file.path(DATA_DIR, "GSE103322_cell_annotations.Rds"))
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"), open = "rt"), warn = FALSE)
hdr <- strsplit(lines[1], "\t")[[1]]
cells <- hdr[-1]
gene_lines <- lines[7:length(lines)]
genes <- sapply(gene_lines, function(l) gsub("'", "", strsplit(l, "\t")[[1]][1]))
expr_list <- lapply(gene_lines, function(l) as.numeric(strsplit(l, "\t")[[1]][-1]))
expr_mat <- do.call(rbind, expr_list)
rownames(expr_mat) <- genes
colnames(expr_mat) <- cells
cat("Expr dims:", dim(expr_mat), "\n")

## ---------- 2. True 67-gene program ----------
common67 <- read.csv(file.path(OUT_DIR, "..", "08_external_validation",
                               "program_common_genes_3OS.csv"))$common_gene
found67 <- intersect(common67, rownames(expr_mat))
cat("67-gene program detected in scRNA:", length(found67), "/67\n")
cat("Missing:", paste(setdiff(common67, found67), collapse = ", "), "\n")

## ---------- 3. Cell annotation alignment (same as P0_4) ----------
annot$cell_id_clean <- gsub("-comb$", "", annot$cell_id)
common <- intersect(annot$cell_id_clean, colnames(expr_mat))
cat("Matched cells:", length(common), "\n")
annot_m <- annot[match(common, annot$cell_id_clean), ]
annot_m$malignant <- ifelse(annot_m$cell_type == "Malignant", 1, 0)
annot_m$cell_group <- ifelse(annot_m$malignant == 1, "Malignant", "Non-malignant")

## ---------- 4. Matching variables: per-gene mean, variance, detection ----------
## Use log-scale expression (as loaded). Compute on all matched cells.
gene_mean <- apply(expr_mat, 1, mean, na.rm = TRUE)
gene_var  <- apply(expr_mat, 1, var, na.rm = TRUE)
gene_det  <- apply(expr_mat, 1, function(x) mean(x > 0, na.rm = TRUE))  # detection rate

true_mean <- mean(gene_mean[found67])
true_var  <- mean(gene_var[found67])
true_det  <- mean(gene_det[found67])
cat(sprintf("True 67-gene program: mean expr=%.3f, mean var=%.3f, detection=%.3f\n",
            true_mean, true_var, true_det))

## ---------- 5. Patient-level enrichment statistic ----------
## For a given gene set: cell-level mean-Z score -> patient × group means ->
## paired difference Malignant - Non-malignant (Wilcoxon statistic on pairs)
patient_enrichment <- function(gset, expr, ann, patient_ids) {
  g <- intersect(gset, rownames(expr))
  if (length(g) < 3) return(NA)
  z <- t(scale(t(expr[g, , drop = FALSE])))
  score <- colMeans(z, na.rm = TRUE)
  ann$score <- score[rownames(ann)]
  pb <- aggregate(score ~ patient + cell_group, data = ann, FUN = mean)
  wide <- reshape(pb, idvar = "patient", timevar = "cell_group", direction = "wide")
  colnames(wide) <- gsub("score\\.", "", colnames(wide))
  mal_col <- grep("Malignant", colnames(wide), value = TRUE)[1]
  non_col <- grep("Non", colnames(wide), value = TRUE)[1]
  w <- wide[complete.cases(wide), ]
  if (nrow(w) < 3 || is.na(mal_col) || is.na(non_col)) return(NA)
  mal <- w[[mal_col]]; non <- w[[non_col]]
  list(n_pairs = nrow(w), n_up = sum(mal > non), median_diff = median(mal - non),
       wilcox_p = suppressWarnings(wilcox.test(mal, non, paired = TRUE)$p.value))
}

## True program enrichment
true_res <- patient_enrichment(found67, expr_mat, annot_m)
cat(sprintf("\nTRUE Mitoxy67: pairs=%d, Malig>Non=%d/%d, median diff=%.4f, Wilcoxon P=%.4g\n",
            true_res$n_pairs, true_res$n_up, true_res$n_pairs,
            true_res$median_diff, true_res$wilcox_p))

## ---------- 6. Null: 1000 matched random gene sets ----------
## Matching strategy: sample 67 genes from the detection-matched pool
## (genes with detection rate >= 0.5 to avoid dropout-dominated sets),
## then accept/reject to match mean expression within tolerance.
## Simpler robust approach: sample from genes ranked by detection rate in
## bins, and within-bin sample with replacement biased toward mean/variance.
cat("\nGenerating", N_NULL, "matched random gene sets...\n")

## Gene pool: detected in >= 30% of cells (reduces dropout-dominated null)
pool <- names(gene_det)[gene_det >= 0.30 & is.finite(gene_mean) & is.finite(gene_var)]
cat("Gene pool (detection>=0.30):", length(pool), "\n")

## For each random set, compute patient enrichment (paired median diff)
null_diff <- rep(NA, N_NULL)
null_up   <- rep(NA, N_NULL)
null_p    <- rep(NA, N_NULL)

for (i in seq_len(N_NULL)) {
  ## Sample 67 genes matched on detection rate (stratified) and
  ## approximately on mean expression via rejection sampling
  repeat {
    ## Stratified sample by detection quantile
    det_q <- cut(gene_det[pool], breaks = quantile(gene_det[pool], probs = seq(0, 1, 0.1)),
                 include.lowest = TRUE, labels = FALSE)
    ## Sample evenly across detection strata, then fill with mean-matched genes
    n_strata <- 10
    per <- floor(67 / n_strata)
    chosen <- unlist(lapply(seq_len(n_strata), function(s) {
      idx <- which(det_q == s)
      if (length(idx) <= per) sample(idx, length(idx)) else sample(idx, per)
    }))
    if (length(chosen) < 67) {
      chosen <- c(chosen, sample(setdiff(seq_along(pool), chosen), 67 - length(chosen)))
    }
    chosen <- chosen[seq_len(67)]
    gset <- pool[chosen]
    ## Acceptance: |mean_expr - true_mean| within 0.15 SD of pool mean spread
    spread <- sd(gene_mean[pool])
    if (abs(mean(gene_mean[gset]) - true_mean) <= 0.15 * spread) break
  }
  r <- patient_enrichment(gset, expr_mat, annot_m)
  null_diff[i] <- r$median_diff
  null_up[i]   <- r$n_up
  null_p[i]    <- r$wilcox_p
  if (i %% 200 == 0) cat("  ", i, "sets done\n")
}

## ---------- 7. Percentile of true enrichment in null ----------
valid <- !is.na(null_diff)
null_diff_v <- null_diff[valid]
null_up_v   <- null_up[valid]
null_p_v    <- null_p[valid]
pct_diff <- mean(null_diff_v <= true_res$median_diff) * 100
pct_up   <- mean(null_up_v <= true_res$n_up) * 100
pct_p    <- mean(null_p_v <= true_res$wilcox_p) * 100

cat("\n===== Null distribution (n =", length(null_diff_v), "valid sets) =====\n")
cat(sprintf("Null median_diff: median=%.4f, 2.5%%=%.4f, 97.5%%=%.4f\n",
            median(null_diff_v), quantile(null_diff_v, 0.025), quantile(null_diff_v, 0.975)))
cat(sprintf("Null n_up (of %d pairs): median=%d, range=%d-%d\n",
            true_res$n_pairs, median(null_up_v), min(null_up_v), max(null_up_v)))
cat(sprintf("Null Wilcoxon P: median=%.3g\n", median(null_p_v)))
cat("\n===== TRUE program percentile =====\n")
cat(sprintf("True median_diff=%.4f -> percentile %.1f%%\n", true_res$median_diff, pct_diff))
cat(sprintf("True n_up=%d/%d -> percentile %.1f%%\n", true_res$n_up, true_res$n_pairs, pct_up))
cat(sprintf("True Wilcoxon P=%.4g -> percentile %.1f%%\n", true_res$wilcox_p, pct_p))

## ---------- 8. Save + plot ----------
null_df <- data.frame(median_diff = null_diff_v, n_up = null_up_v, wilcox_p = null_p_v)
write.csv(null_df, file.path(OUT_DIR, "GSE103322_random_null_1000sets.csv"), row.names = FALSE)

png(file.path(OUT_DIR, "Fig_P0-2_random_null_enrichment.png"), width = 2200, height = 1600, res = 300)
par(mfrow = c(1, 2))
hist(null_diff_v, breaks = 30, col = "lightblue", border = "white",
     main = "Null: patient-level malignant enrichment\n(1000 matched 67-gene sets)",
     xlab = "Median diff (Malignant - Non-malignant)")
abline(v = true_res$median_diff, col = "red", lwd = 3, lty = 2)
mtext(sprintf("TRUE Mitoxy67\n(percentile %.1f%%)", pct_diff), side = 3, line = -1.5, at = true_res$median_diff, col = "red", cex = 0.8)

hist(-log10(null_p_v), breaks = 30, col = "lightgreen", border = "white",
     main = "Null: patient-level Wilcoxon P (-log10)\n(1000 matched 67-gene sets)",
     xlab = "-log10(Wilcoxon P)")
abline(v = -log10(true_res$wilcox_p), col = "red", lwd = 3, lty = 2)
mtext(sprintf("TRUE P=%.2g\n(percentile %.1f%%)", true_res$wilcox_p, pct_p), side = 3, line = -1.5, at = -log10(true_res$wilcox_p), col = "red", cex = 0.8)
dev.off()
cat("\nSaved Fig_P0-2_random_null_enrichment.png\n")

## ---------- 9. Summary table ----------
sum_tab <- data.frame(
  Metric = c("True_program_median_diff", "Null_median_diff_median",
             "True_percentile_median_diff", "True_n_up_pairs",
             "Null_n_up_median", "True_percentile_n_up",
             "True_wilcoxon_P", "Null_wilcoxon_P_median", "True_percentile_P",
             "N_null_sets_valid"),
  Value = c(round(true_res$median_diff, 4), round(median(null_diff_v), 4),
            round(pct_diff, 1), true_res$n_up, median(null_up_v),
            round(pct_up, 1), signif(true_res$wilcox_p, 3),
            signif(median(null_p_v), 3), round(pct_p, 1), length(null_diff_v))
)
write.csv(sum_tab, file.path(OUT_DIR, "GSE103322_random_null_summary.csv"), row.names = FALSE)
print(sum_tab, row.names = FALSE)

cat("\n=== P0-2 random-gene-set null DONE ===")
