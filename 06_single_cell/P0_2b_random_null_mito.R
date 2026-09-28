## ============================================================
## P0-2 Null-2 (efficient version): matched random-gene-set null
## with mitochondrial + ribosomal annotation proportion matching.
## Stratified sampling (no rejection loop): genes stratified by
## {mito, ribo} x detection-decile; sample to match target
## proportions and expression statistics.
## ============================================================
suppressMessages({
  library(data.table); library(org.Hs.eg.db)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/11_singlecell"
set.seed(20260902)
N_NULL <- 1000

## ---------- 1. Load GSE103322 ----------
annot <- readRDS(file.path(DATA_DIR, "GSE103322_cell_annotations.Rds"))
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"), open = "rt"), warn = FALSE)
hdr <- strsplit(lines[1], "\t")[[1]]
cells <- hdr[-1]
gene_lines <- lines[7:length(lines)]
genes <- sapply(gene_lines, function(l) gsub("'", "", strsplit(l, "\t")[[1]][1]))
expr_list <- lapply(gene_lines, function(l) as.numeric(strsplit(l, "\t")[[1]][-1]))
expr_mat <- do.call(rbind, expr_list)
rownames(expr_mat) <- genes; colnames(expr_mat) <- cells

common67 <- read.csv(file.path(OUT_DIR, "..", "08_external_validation",
                               "program_common_genes_3OS.csv"))$common_gene
found67 <- intersect(common67, rownames(expr_mat))

annot$cell_id_clean <- gsub("-comb$", "", annot$cell_id)
common <- intersect(annot$cell_id_clean, colnames(expr_mat))
annot_m <- annot[match(common, annot$cell_id_clean), ]
annot_m$malignant <- ifelse(annot_m$cell_type == "Malignant", 1, 0)
annot_m$cell_group <- ifelse(annot_m$malignant == 1, "Malignant", "Non-malignant")

## ---------- 2. Annotation flags ----------
mito_go <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = "GO:0005739",
                                 column = "SYMBOL", keytype = "GO", multiVals = "list")
mito_genes <- unique(unlist(mito_go))
ribo_go <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = "GO:0005840",
                                 column = "SYMBOL", keytype = "GO", multiVals = "list")
ribo_genes <- unique(unlist(ribo_go))
cat("Mito:", length(mito_genes), "Ribo:", length(ribo_genes), "\n")

is_mito <- function(g) g %in% mito_genes
is_ribo <- function(g) g %in% ribo_genes
true_mito <- mean(is_mito(found67)); true_ribo <- mean(is_ribo(found67))
cat(sprintf("TRUE: mito prop=%.3f, ribo prop=%.3f\n", true_mito, true_ribo))

## ---------- 3. Gene stats + pool ----------
gene_mean <- apply(expr_mat, 1, mean, na.rm = TRUE)
gene_var  <- apply(expr_mat, 1, var, na.rm = TRUE)
gene_det  <- apply(expr_mat, 1, function(x) mean(x > 0, na.rm = TRUE))
pool <- names(gene_det)[gene_det >= 0.30 & is.finite(gene_mean) & is.finite(gene_var)]
cat("Pool:", length(pool), "\n")

true_mean <- mean(gene_mean[found67]); true_var <- mean(gene_var[found67]); true_det <- mean(gene_det[found67])

## ---------- 4. Stratified matched sampling (no rejection) ----------
## Targets: mito count, ribo count (round), detection quantile-balanced, mean-matched
target_n_mito <- round(true_mito * 67)
target_n_ribo <- round(true_ribo * 67)
cat(sprintf("Target counts (67): mito=%d, ribo=%d\n", target_n_mito, target_n_ribo))

## Mitophagy genes overlap with ribo? assume disjoint categories: mito-only, ribo-only, both, neither
in_mito_pool <- pool %in% mito_genes
in_ribo_pool <- pool %in% ribo_genes
cat("Pool composition: mito+ribo+ =", sum(in_mito_pool & in_ribo_pool),
    "mito+ribo- =", sum(in_mito_pool & !in_ribo_pool),
    "mito-ribo+ =", sum(!in_mito_pool & in_ribo_pool),
    "neither =", sum(!in_mito_pool & !in_ribo_pool), "\n")

sample_stratified <- function() {
  ## Strict two-layer stratified sampling:
  ## Layer 1: mito genes (target n_mito) vs non-mito genes (67 - n_mito)
  ## Layer 2: within each layer, balance detection deciles
  n_mito <- target_n_mito
  n_non <- 67 - n_mito
  mito_idx <- which(in_mito_pool)
  non_idx  <- which(!in_mito_pool)

  ## Layer-2 helper: sample n genes from idx balanced across detection deciles
  sample_balanced <- function(idx, n) {
    if (length(idx) <= n) return(idx)
    det_vals <- gene_det[pool[idx]]
    brk <- quantile(det_vals, probs = seq(0, 1, 0.2), na.rm = TRUE)  # 5 deciles
    q <- cut(det_vals, breaks = unique(brk), include.lowest = TRUE, labels = FALSE)
    per <- floor(n / length(unique(q)))
    out <- unlist(lapply(sort(unique(q)), function(d) {
      ii <- idx[q == d]
      if (length(ii) <= per) ii else sample(ii, per)
    }))
    if (length(out) < n) out <- c(out, sample(setdiff(idx, out), n - length(out)))
    out[seq_len(n)]
  }

  sel_mito <- sample_balanced(mito_idx, n_mito)
  sel_non  <- sample_balanced(non_idx, n_non)
  sel <- c(sel_mito, sel_non)
  sel <- sel[seq_len(67)]
  gset <- pool[sel]
  list(gset = gset, mito = mean(in_mito_pool[sel]), ribo = mean(in_ribo_pool[sel]))
}

## ---------- 5. Patient enrichment ----------
patient_enrichment <- function(gset, expr, ann) {
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

true_res <- patient_enrichment(found67, expr_mat, annot_m)
cat(sprintf("TRUE: pairs=%d, up=%d/%d, diff=%.4f, P=%.4g\n",
            true_res$n_pairs, true_res$n_up, true_res$n_pairs,
            true_res$median_diff, true_res$wilcox_p))

## ---------- 6. Run 1000 null-2 ----------
cat("Generating Null-2...\n")
null2_diff <- numeric(N_NULL); null2_up <- integer(N_NULL); null2_p <- numeric(N_NULL)
null2_mito <- numeric(N_NULL); null2_ribo <- numeric(N_NULL)
null2_mean <- numeric(N_NULL); null2_var <- numeric(N_NULL); null2_det <- numeric(N_NULL)

for (i in seq_len(N_NULL)) {
  sm <- sample_stratified()
  gset <- sm$gset
  r <- patient_enrichment(gset, expr_mat, annot_m)
  null2_diff[i] <- r$median_diff; null2_up[i] <- r$n_up; null2_p[i] <- r$wilcox_p
  null2_mito[i] <- sm$mito; null2_ribo[i] <- sm$ribo
  null2_mean[i] <- mean(gene_mean[gset]); null2_var[i] <- mean(gene_var[gset]); null2_det[i] <- mean(gene_det[gset])
  if (i %% 200 == 0) cat("  ", i, "done\n")
}

valid <- !is.na(null2_diff)
pct_diff <- mean(null2_diff[valid] <= true_res$median_diff) * 100
pct_up   <- mean(null2_up[valid] <= true_res$n_up) * 100

cat("\n===== Null-2 (mito+ribo stratified) =====\n")
cat(sprintf("Matched: mito %.3f (tgt %.3f), ribo %.3f (tgt %.3f)\n",
            mean(null2_mito), true_mito, mean(null2_ribo), true_ribo))
cat(sprintf("Stats: mean %.3f (tgt %.3f), det %.3f (tgt %.3f)\n",
            mean(null2_mean), true_mean, mean(null2_det), true_det))
cat(sprintf("Null-2 median_diff: median=%.4f, 97.5%%=%.4f\n",
            median(null2_diff[valid]), quantile(null2_diff[valid], 0.975)))
cat(sprintf("TRUE diff=%.4f -> percentile %.1f%%\n", true_res$median_diff, pct_diff))
cat(sprintf("TRUE up=%d/%d -> percentile %.1f%%\n", true_res$n_up, true_res$n_pairs, pct_up))

## ---------- 7. Save + plot ----------
null2_df <- data.frame(median_diff = null2_diff[valid], n_up = null2_up[valid],
                       wilcox_p = null2_p[valid], mito_prop = null2_mito[valid],
                       ribo_prop = null2_ribo[valid], mean_expr = null2_mean[valid],
                       var_expr = null2_var[valid], detection = null2_det[valid])
write.csv(null2_df, file.path(OUT_DIR, "GSE103322_random_null2_1000sets.csv"), row.names = FALSE)

png(file.path(OUT_DIR, "Fig_P0-2b_random_null2_enrichment.png"), width = 2200, height = 1600, res = 300)
par(mfrow = c(1, 2))
hist(null2_df$median_diff, breaks = 30, col = "lightblue", border = "white",
     main = "Null-2 (mito+ribo matched): enrichment\n(1000 matched 67-gene sets)",
     xlab = "Median diff (Malignant - Non-malignant)")
abline(v = true_res$median_diff, col = "red", lwd = 3, lty = 2)
mtext(sprintf("TRUE Mitoxy67\n(percentile %.1f%%)", pct_diff), side = 3, line = -1.5,
      at = true_res$median_diff, col = "red", cex = 0.8)
hist(-log10(null2_df$wilcox_p), breaks = 30, col = "lightgreen", border = "white",
     main = "Null-2: Wilcoxon P (-log10)", xlab = "-log10(P)")
abline(v = -log10(true_res$wilcox_p), col = "red", lwd = 3, lty = 2)
dev.off()
cat("Saved Fig_P0-2b_random_null2_enrichment.png\n")

sum2 <- data.frame(
  Metric = c("True_median_diff", "Null2_median", "Null2_97.5", "Percentile_diff",
             "True_n_up", "Null2_n_up_median", "Percentile_up",
             "Matched_mito", "Target_mito", "Matched_ribo", "Target_ribo",
             "Matched_mean", "Target_mean", "Matched_det", "Target_det", "N_valid"),
  Value = c(round(true_res$median_diff, 4), round(median(null2_df$median_diff), 4),
            round(quantile(null2_df$median_diff, 0.975), 4), round(pct_diff, 1),
            true_res$n_up, median(null2_df$n_up), round(pct_up, 1),
            round(mean(null2_df$mito_prop), 3), round(true_mito, 3),
            round(mean(null2_df$ribo_prop), 3), round(true_ribo, 3),
            round(mean(null2_df$mean_expr), 3), round(true_mean, 3),
            round(mean(null2_df$detection), 3), round(true_det, 3), length(valid))
)
write.csv(sum2, file.path(OUT_DIR, "GSE103322_random_null2_summary.csv"), row.names = FALSE)
print(sum2, row.names = FALSE)

cat("\n=== P0-2 Null-2 DONE ===")
