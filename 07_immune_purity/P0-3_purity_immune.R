## ============================================================
## P0-3: Tumor-purity-adjusted immune analysis
## GPT scholar 5.2 MC24 / MC6
## Are the immune associations (CD4 memory resting rho=0.26,
## PD-L1 rho=0.183) driven by tumor purity / composition?
## 1) ESTIMATE tumor purity (immune + stromal scores)
## 2) Partial correlation: cor(Program, CD4Memory | purity)
## 3) Regression: CD4Memory ~ Program + Purity + Site
##               CD274 ~ Program + IFNG + Purity + Site
## ============================================================
suppressMessages({
  library(survival); library(estimate)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------- 1. Load data ----------
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat),
                 OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm],
                 age = clin$age[cm], gender = clin$gender[cm],
                 stage = clin$ajcc_stage[cm], stringsAsFactors = FALSE)
df <- df[!is.na(df$OS_time), ]
cat("Cohort N =", nrow(df), "\n")

## Site mapping
site_map <- read.delim(file.path(OUT_DIR, "..", "08_external_validation", "tcga_hnsc_primary_site.tsv"),
                       check.names = FALSE, quote = "")
site_map$patient <- substr(site_map$submitter_id, 1, 12)
df$patient <- substr(df$sample, 1, 12)
sm <- match(df$patient, site_map$patient)
ps <- tolower(site_map$primary_site[sm])
df$site <- ifelse(grepl("tongue|mouth|gum|palate|lip|floor", ps), "Oral cavity",
           ifelse(grepl("larynx|hypopharyn", ps), "Larynx/Hypopharynx",
           ifelse(grepl("tonsil|oropharyn|base of tongue", ps), "Oropharynx",
           "Other")))

## ---------- 2. Mitoxy67 score ----------
common67 <- read.csv(file.path(OUT_DIR, "..", "08_external_validation",
                               "program_common_genes_3OS.csv"))$common_gene
calc_score <- function(genes, expr_mat, samples) {
  g <- intersect(genes, rownames(expr_mat))
  ex <- log2(expr_mat[g, samples, drop = FALSE] + 1)
  z <- t(scale(t(ex)))
  colMeans(z, na.rm = TRUE)
}
df$mitoxy67 <- calc_score(common67, tumor_mat, df$sample)

## ---------- 3. CIBERSORT fractions (502 x 25) ----------
cib <- readRDS(file.path(OUT_DIR, "..", "20_figures_main", "Step8_cibersort_502.Rds"))
cib_df <- as.data.frame(cib)
cib_df$sample <- rownames(cib)
df <- merge(df, cib_df, by = "sample", all.x = TRUE)
cat("Merged CIBERSORT samples:", sum(!is.na(df[["T cells CD4 memory resting"]])), "\n")
## CD4 memory resting column name (contains spaces)
cd4_col <- "T cells CD4 memory resting"
cat("CD4 memory resting column:", cd4_col, "\n")

## ---------- 4. ESTIMATE purity ----------
## Write input matrix (gene symbol x sample, linear scale), de-duplicate genes
## by median expression (1233 duplicated symbols in TCGA matrix)
est_mat <- 2^tumor_mat - 1  # back to linear scale
est_mat[est_mat < 0] <- 0
est_mat <- est_mat[!duplicated(rownames(est_mat)), , drop = FALSE]
cat("ESTIMATE input genes after dedup:", nrow(est_mat), "\n")
write.table(est_mat, file.path(OUT_DIR, "P0-3_estimate_input.txt"),
            sep = "\t", quote = FALSE)
filterCommonGenes(input.f = file.path(OUT_DIR, "P0-3_estimate_input.txt"),
                  output.f = file.path(OUT_DIR, "P0-3_estimate_genes.gct"),
                  id = "GeneSymbol")
estimateScore(file.path(OUT_DIR, "P0-3_estimate_genes.gct"),
              output.ds = file.path(OUT_DIR, "P0-3_estimate_scores.gct"),
              platform = "illumina")
## Read ESTIMATE scores (skip=2: #1.2 + dimension line; line 3 = header)
est_dat <- read.delim(file.path(OUT_DIR, "P0-3_estimate_scores.gct"),
                      skip = 2, row.names = 1, check.names = FALSE)
## est_dat rows = StromalScore/ImmuneScore/ESTIMATEScore; cols = samples
est_scores <- as.data.frame(t(est_dat[, -1]))  # drop Description col
est_scores$sample_raw <- rownames(est_scores)
est_scores$sample <- gsub("\\.", "-", est_scores$sample_raw)
cat("ESTIMATE scores samples:", nrow(est_scores), "cols:", paste(colnames(est_scores), collapse=","), "\n")
df <- merge(df, est_scores, by = "sample", all.x = TRUE)
df$est_purity <- cos(0.6049872018 + 0.0001467884 * df$ESTIMATEScore)
cat("ESTIMATE purity: median =", round(median(df$est_purity, na.rm = TRUE), 3),
    "range =", round(range(df$est_purity, na.rm = TRUE), 3), "\n")

## ---------- 5. Partial correlation: cor(Program, CD4Mem | purity) ----------
## Residual-based partial Spearman (no ppcor dependency):
## rank-transform both variables, regress out purity ranks, correlate residuals
partial_spearman <- function(x, y, z) {
  rx <- rank(x); ry <- rank(y); rz <- rank(z)
  ## residuals of rx on rz, ry on rz
  rx_res <- residuals(lm(rx ~ rz))
  ry_res <- residuals(lm(ry ~ rz))
  r <- cor(rx_res, ry_res)
  n <- length(rx)
  df_adj <- n - 3
  t_stat <- r * sqrt(df_adj / (1 - r^2))
  p_val <- 2 * pt(-abs(t_stat), df = df_adj)
  list(rho = r, p.value = p_val, n = n)
}

{
  cd4 <- df[[cd4_col]]
  keep <- complete.cases(df$mitoxy67, cd4, df$est_purity)
  cat("Complete cases for partial cor:", sum(keep), "\n")
  if (sum(keep) > 50) {
    pcor_res <- partial_spearman(df$mitoxy67[keep], cd4[keep], df$est_purity[keep])
    cat(sprintf("Partial cor (Spearman) Program vs CD4Mem | purity: rho=%.3f, P=%.4g (n=%d)\n",
                pcor_res$rho, pcor_res$p.value, pcor_res$n))
  }
}

## ---------- 6. Regression models ----------
cat("\n===== Regression models =====\n")

## CD4Memory ~ Program + Purity + Site (column name has spaces -> backticks)
cd4_keep <- complete.cases(df$mitoxy67, df[[cd4_col]], df$est_purity, df$site)
cat("Complete cases for CD4 regression:", sum(cd4_keep), "\n")
d_cd4 <- df[cd4_keep, ]
d_cd4$cd4_mem <- d_cd4[[cd4_col]]
fit_cd4_uni <- lm(cd4_mem ~ mitoxy67, data = d_cd4)
fit_cd4_adj <- lm(cd4_mem ~ mitoxy67 + est_purity + site, data = d_cd4)
s1 <- summary(fit_cd4_uni); s2 <- summary(fit_cd4_adj)
cat(sprintf("CD4Mem ~ Mitoxy (unadjusted):      beta=%.4f, P=%.4g\n",
            s1$coefficients["mitoxy67", 1], s1$coefficients["mitoxy67", 4]))
cat(sprintf("CD4Mem ~ Mitoxy + Purity + Site:   beta=%.4f, P=%.4g\n",
            s2$coefficients["mitoxy67", 1], s2$coefficients["mitoxy67", 4]))

## CD274 ~ Program + IFNG + Purity + Site
if ("CD274" %in% rownames(tumor_mat) && "IFNG" %in% rownames(tumor_mat)) {
  df$CD274 <- as.numeric(log2(tumor_mat["CD274", df$sample] + 1))
  df$IFNG  <- as.numeric(log2(tumor_mat["IFNG", df$sample] + 1))
  pd_keep <- complete.cases(df$mitoxy67, df$CD274, df$IFNG, df$est_purity, df$site)
  cat("Complete cases for CD274 regression:", sum(pd_keep), "\n")
  d_pd <- df[pd_keep, ]
  fit_pd_uni <- lm(CD274 ~ mitoxy67, data = d_pd)
  fit_pd_adj <- lm(CD274 ~ mitoxy67 + IFNG + est_purity + site, data = d_pd)
  p1 <- summary(fit_pd_uni); p2 <- summary(fit_pd_adj)
  cat(sprintf("CD274 ~ Mitoxy (unadjusted):          beta=%.4f, P=%.4g\n",
              p1$coefficients["mitoxy67", 1], p1$coefficients["mitoxy67", 4]))
  cat(sprintf("CD274 ~ Mitoxy + IFNG + Purity + Site: beta=%.4f, P=%.4g\n",
              p2$coefficients["mitoxy67", 1], p2$coefficients["mitoxy67", 4]))
}

## ---------- 7. Spearman raw vs partial summary ----------
cat("\n===== Summary: raw vs purity-adjusted =====\n")
raw_cd4 <- cor(df$mitoxy67, df[[cd4_col]], method = "spearman", use = "complete.obs")
cat(sprintf("Raw Spearman: Program vs CD4Mem rho=%.3f\n", raw_cd4))
if (exists("pcor_res")) {
  cat(sprintf("Partial (|purity): rho=%.3f P=%.4g\n", pcor_res$estimate, pcor_res$p.value))
}
if (exists("fit_pd_uni")) {
  raw_pd <- cor(df$mitoxy67, df$CD274, method = "spearman", use = "complete.obs")
  cat(sprintf("Raw Spearman: Program vs CD274 rho=%.3f\n", raw_pd))
}

## ---------- 8. Save results ----------
res_list <- list(
  CD4_raw_spearman = if (exists("raw_cd4")) round(raw_cd4, 4) else NA,
  CD4_partial_rho = if (exists("pcor_res")) round(pcor_res$rho, 4) else NA,
  CD4_partial_P = if (exists("pcor_res")) signif(pcor_res$p.value, 3) else NA,
  CD4_reg_beta_uni = if (exists("s1")) round(s1$coefficients["mitoxy67", 1], 5) else NA,
  CD4_reg_P_uni = if (exists("s1")) signif(s1$coefficients["mitoxy67", 4], 3) else NA,
  CD4_reg_beta_adj = if (exists("s2")) round(s2$coefficients["mitoxy67", 1], 5) else NA,
  CD4_reg_P_adj = if (exists("s2")) signif(s2$coefficients["mitoxy67", 4], 3) else NA,
  CD274_raw_spearman = if (exists("raw_pd")) round(raw_pd, 4) else NA,
  CD274_reg_beta_uni = if (exists("p1")) round(p1$coefficients["mitoxy67", 1], 5) else NA,
  CD274_reg_P_uni = if (exists("p1")) signif(p1$coefficients["mitoxy67", 4], 3) else NA,
  CD274_reg_beta_adj = if (exists("p2")) round(p2$coefficients["mitoxy67", 1], 5) else NA,
  CD274_reg_P_adj = if (exists("p2")) signif(p2$coefficients["mitoxy67", 4], 3) else NA,
  ESTIMATE_purity_median = round(median(df$est_purity, na.rm = TRUE), 3)
)
res_df <- data.frame(Metric = names(res_list), Value = unlist(res_list))
write.csv(res_df, file.path(OUT_DIR, "P0-3_purity_adjusted_immune.csv"), row.names = FALSE)
print(res_df, row.names = FALSE)

cat("\n=== P0-3 purity-adjusted immune DONE ===")
