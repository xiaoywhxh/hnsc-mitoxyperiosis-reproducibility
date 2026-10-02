# Pathology-based independent validation of malignant localization.
# pathologist_anno_raw: SCC (malignant) vs pre-defined clearly non-SCC tissue.
#
# PRE-DEFINED non-SCC comparator (declared BEFORE looking at scores):
#   INCLUDED (host/normal tissue):
#     Lymphocyte Negative Stroma, Lymphocyte Positive Stroma, Muscle,
#     Glandular Stroma, Non-cancerous Mucosa, Lymphocyte Positive Muscles,
#     Artery/Vein
#   EXCLUDED (not tissue / technical / ambiguous, to avoid comparator contamination):
#     Artifact, Cautery, Fold, Edge Effects  (technical artefacts)
#     Keratin                                (keratinisation can be tumour-associated)
#
# Per sample: median difference, Cliff's delta. Cross-sample: sign test + RE meta.
# sample_3: uses only its 476 coordinate-matched spots (as is).
suppressPackageStartupMessages({ library(data.table); library(metafor) })
IN  <- "D:/GSE208253_spatial/spatial_analysis/input"
SC  <- "D:/GSE208253_spatial/spatial_analysis/score"
OUT <- "D:/GSE208253_spatial/spatial_analysis/effects"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

spots  <- fread(file.path(IN, "spot_table_all.tsv"), sep = "\t", header = TRUE)
scores <- fread(file.path(SC, "spot_scores.tsv"), sep = "\t", header = TRUE)
d <- merge(spots, scores, by = c("sample", "barcode"), all.x = TRUE)

NON_SCC <- c("Lymphocyte Negative Stroma", "Lymphocyte Positive Stroma",
             "Muscle", "Glandular Stroma", "Non-cancerous Mucosa",
             "Lymphocyte Positive Muscles", "Artery/Vein")
EXCLUDED <- c("Artifact", "Cautery", "Fold", "Edge Effects", "Keratin")

d[, pclass := fifelse(pathologist_anno_raw == "SCC", "SCC",
              fifelse(pathologist_anno_raw %in% NON_SCC, "nonSCC",
              fifelse(pathologist_anno_raw %in% EXCLUDED, "excluded",
              fifelse(pathologist_anno_raw == "", "unannotated", "other"))))]

cat("=== spot class counts ===\n"); print(d[, .N, by = pclass][order(-N)])

cliff_delta <- function(a, b) {
  n1 <- length(a); n2 <- length(b); if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(a, b)); U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  (2 * U) / (n1 * n2) - 1
}

rows <- list()
for (s in paste0("s", 1:12)) {
  da <- d[sample == s & pclass == "SCC"]
  db <- d[sample == s & pclass == "nonSCC"]
  if (nrow(da) < 3 || nrow(db) < 3) next
  wt <- suppressWarnings(wilcox.test(da$ucell, db$ucell, exact = FALSE))
  rows[[s]] <- data.table(sample = s, n_SCC = nrow(da), n_nonSCC = nrow(db),
    median_SCC = median(da$ucell), median_nonSCC = median(db$ucell),
    median_diff = median(da$ucell) - median(db$ucell),
    cliffs_delta = cliff_delta(da$ucell, db$ucell), MWU_p = wt$p.value)
}
pe <- rbindlist(rows)
fwrite(pe, file.path(OUT, "pathology_effect_sizes.tsv"), sep = "\t")

n <- nrow(pe); npos <- sum(pe$cliffs_delta > 0)
st <- binom.test(npos, n, 0.5)
pe[, se := sqrt((n_SCC + n_nonSCC + 1) / (3 * n_SCC * n_nonSCC))]
m <- rma(yi = cliffs_delta, sei = se, data = as.data.frame(pe), method = "REML")
sump <- data.table(
  n_samples = n, n_delta_positive = npos, n_delta_negative = sum(pe$cliffs_delta < 0),
  sign_test_p = st$p.value,
  median_of_median_diff = median(pe$median_diff), mean_delta = mean(pe$cliffs_delta),
  meta_delta = as.numeric(m$beta), meta_se = m$se,
  meta_ci_lo = m$ci.lb, meta_ci_hi = m$ci.ub, meta_p = m$pval,
  meta_tau2 = m$tau2, meta_I2 = m$I2)
fwrite(sump, file.path(OUT, "pathology_cross_sample_summary.tsv"), sep = "\t")

# leave-s3-out
pe3 <- pe[sample != "s3"]; n3 <- nrow(pe3); p3 <- sum(pe3$cliffs_delta > 0)
m3 <- rma(yi = cliffs_delta, sei = se, data = as.data.frame(pe3), method = "REML")
sum3 <- data.table(n_samples = n3, n_delta_positive = p3,
  sign_test_p = binom.test(p3, n3, 0.5)$p.value,
  meta_delta = as.numeric(m3$beta), meta_p = m3$pval)
fwrite(sum3, file.path(OUT, "pathology_leave_s3_out_summary.tsv"), sep = "\t")

cat("\n=== PATHOLOGY EFFECT SIZES (SCC vs non-SCC) ===\n")
print(pe[, .(sample, n_SCC, n_nonSCC, median_SCC = round(median_SCC,3),
             median_nonSCC = round(median_nonSCC,3),
             median_diff = round(median_diff,3), cliffs_delta = round(cliffs_delta,3))])
cat("\n=== CROSS-SAMPLE ===\n"); print(sump)
cat("\n=== LEAVE-s3-OUT ===\n"); print(sum3)
