## ============================================================
## P0-1: Specificity / null-model benchmark (GPT scholar 5.2 MC17)
## MitoxyScore vs 9 established biological signatures:
##   Hallmark OXPHOS / Glycolysis / Hypoxia / PI3K-AKT-mTOR /
##   Apoptosis / ROS / cell-cycle(proliferation) / mitophagy /
##   mitochondrial dynamics
## 1) Score correlation matrix (MitoxyScore vs each signature)
## 2) Multivariable Cox:
##    OS ~ MitoxyScore + OXPHOS + Apoptosis + Glycolysis +
##        Proliferation + Age + Sex + Stage + Site
## 3) PD-L1: CD274 ~ MitoxyScore + IFNG + purity + site (purity in P0-3)
## Output: CSV + correlation heatmap PNG
## ============================================================
suppressMessages({
  library(survival); library(org.Hs.eg.db)
  library(pheatmap); library(timeROC)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"
GMT_FILE <- file.path(OUT_DIR, "msigdb", "h.all.v2024.1.Hs.symbols.gmt")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ---------- 0. Read MSigDB Hallmark GMT (local file) ----------
read_gmt <- function(path) {
  lines <- readLines(path)
  res <- lapply(lines, function(l) {
    parts <- strsplit(l, "\t")[[1]]
    list(name = parts[1], desc = parts[2], genes = parts[-(1:2)])
  })
  sets <- lapply(res, function(x) x$genes)
  names(sets) <- sapply(res, function(x) x$name)
  sets
}
hall_list <- read_gmt(GMT_FILE)
cat("Hallmark sets loaded:", length(hall_list), "\n")

## ---------- 1. Load data (same pipeline as V5_P0_1_tcga67.R) ----------
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
cat("Cohort N =", nrow(df), "events =", sum(df$OS_status), "\n")

## ---------- 2. Site mapping (same as P1_supp_program_adjusted.R) ----------
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
cat("Site distribution:\n"); print(table(df$site, useNA = "ifany"))

## ---------- 3. Score computation helper (program-mean Z, same rule) ----------
calc_score <- function(genes, expr_mat, samples) {
  g <- intersect(genes, rownames(expr_mat))
  if (length(g) < 3) return(rep(NA, length(samples)))
  ex <- log2(expr_mat[g, samples, drop = FALSE] + 1)
  z <- t(scale(t(ex)))
  colMeans(z, na.rm = TRUE)
}

## ---------- 4. Mitoxy program (67-gene common) ----------
common67 <- read.csv(file.path(OUT_DIR, "..", "08_external_validation",
                               "program_common_genes_3OS.csv"))$common_gene
df$mitoxy67 <- calc_score(common67, tumor_mat, df$sample)
cat("Mitoxy 67-gene score: N valid =", sum(!is.na(df$mitoxy67)), "\n")

## ---------- 5. 9 comparison signatures via local Hallmark GMT ----------
## NOTE (2026-09-01 review): E2F_TARGETS used as transcriptional proxy for
## proliferative state (cell-cycle / DNA replication); not a complete
## proliferation signature (alternatives: MKI67 metagene, CIN70, PAM50 proliferation).
sig_names <- c(
  OXPHOS        = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  Glycolysis    = "HALLMARK_GLYCOLYSIS",
  Hypoxia       = "HALLMARK_HYPOXIA",
  PI3K_mTOR     = "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  Apoptosis     = "HALLMARK_APOPTOSIS",
  ROS           = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  Proliferation = "HALLMARK_E2F_TARGETS"   # proxy for proliferative state
)
## mitophagy + mitochondrial dynamics: literature-based (not in Hallmark)
## Sources: PINK1/Parkin canonical machinery (Pickrell & Youle 2015), BNIP3/BNIP3L
## hypoxia-induced mitophagy (Zhang & Ney 2009), FUNDC1 (Liu et al 2012),
## receptor OPTN/NBR1/TAX1BP1 (Lazarou et al 2015), core ATG machinery.
mitophagy_genes <- c("PINK1","PRKN","PARK7","BNIP3","BNIP3L","FUNDC1","SQSTM1",
                     "OPTN","NBR1","TAX1BP1","ATG5","ATG7","LC3B","MAP1LC3A",
                     "MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2","ULK1","ULK2")
## CYCS (cytochrome c) participates in both mitochondrial respiration and
## apoptosis execution; retained but flagged as a dynamics/apoptosis interface gene.
mitodyn_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1",
                   "MARCHF5","MARCH5","SLC25A3","SLC25A4","SLC25A5","SLC25A6",
                   "PHB","PHB2","OMA1","YME1L1")   # CYCS excluded from core dynamics set

## Evidence-source table for literature-based signatures (Supplementary)
mito_src <- data.frame(
  Gene = mitophagy_genes,
  Signature = "Mitophagy",
  Evidence_source = rep("", length(mitophagy_genes)),
  stringsAsFactors = FALSE
)
mito_src$Evidence_source[mito_src$Gene %in% c("PINK1","PRKN")] <- "Canonical PINK1/Parkin machinery (Pickrell & Youle 2015)"
mito_src$Evidence_source[mito_src$Gene == "PARK7"] <- "Parkinsonism-associated deglycase (Parkin pathway)"
mito_src$Evidence_source[mito_src$Gene %in% c("BNIP3","BNIP3L")] <- "Hypoxia-induced mitophagy receptor (Zhang & Ney 2009)"
mito_src$Evidence_source[mito_src$Gene == "FUNDC1"] <- "Hypoxia receptor (Liu et al 2012)"
mito_src$Evidence_source[mito_src$Gene %in% c("SQSTM1","OPTN","NBR1","TAX1BP1")] <- "Ubiquitin/LC3 receptor (Lazarou et al 2015)"
mito_src$Evidence_source[mito_src$Gene %in% c("ATG5","ATG7","LC3B")] <- "Core ATG conjugation machinery"
mito_src$Evidence_source[mito_src$Gene %in% c("MAP1LC3A","MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2")] <- "LC3/GABARAP family"
mito_src$Evidence_source[mito_src$Gene %in% c("ULK1","ULK2")] <- "ULK complex"
write.csv(mito_src, file.path(OUT_DIR, "P0-1_mitophagy_gene_sources.csv"), row.names = FALSE)

mitodyn_src <- data.frame(
  Gene = mitodyn_genes,
  Signature = "MitoDynamics",
  Evidence_source = c(rep("Fission machinery (DNM1L/FIS1/MFF/MIEF1/2)", 5),
                      rep("Fusion machinery (MFN1/2/OPA1)", 3),
                      rep("E3 ligase/adaptor (MARCHF5/MARCH5)", 2),
                      rep("Mitochondrial carriers (SLC25A3-6)", 4),
                      rep("Scaffold/proteostasis (PHB/PHB2/OMA1/YME1L1)", 4)),
  stringsAsFactors = FALSE
)
write.csv(mitodyn_src, file.path(OUT_DIR, "P0-1_mitodyn_gene_sources.csv"), row.names = FALSE)

## Build score matrix
sig_genes <- lapply(sig_names, function(nm) hall_list[[nm]])
sig_genes$Mitophagy <- mitophagy_genes
sig_genes$MitoDynamics <- mitodyn_genes
names(sig_genes) <- c(names(sig_names), "Mitophagy", "MitoDynamics")
cat("Signature gene counts:\n"); print(sapply(sig_genes, length))

score_mat <- sapply(names(sig_genes), function(nm) {
  calc_score(sig_genes[[nm]], tumor_mat, df$sample)
})
df <- cbind(df, score_mat)
cat("\nScore matrix dim:", dim(score_mat), "\n")

## ---------- 6. Correlation matrix: Mitoxy67 vs all signatures ----------
cor_vars <- c("mitoxy67", names(sig_genes))
cor_mat <- cor(df[, cor_vars], use = "pairwise.complete.obs")
cat("\n===== Correlation with Mitoxy67 =====\n")
print(round(cor_mat["mitoxy67", ], 3))

## Spearman (more robust) + P values
sp_rho <- sapply(names(sig_genes), function(nm) {
  cor(df$mitoxy67, df[[nm]], method = "spearman", use = "complete.obs")
})
sp_p <- sapply(names(sig_genes), function(nm) {
  ct <- cor.test(df$mitoxy67, df[[nm]], method = "spearman")
  ct$p.value
})
cor_res <- data.frame(Signature = names(sig_genes),
                      Pearson_r = round(cor_mat["mitoxy67", names(sig_genes)], 4),
                      Spearman_rho = round(sp_rho, 4),
                      Spearman_P = formatC(sp_p, format = "e", digits = 2))
cat("\n===== Spearman with P =====\n"); print(cor_res, row.names = FALSE)

## Correlation heatmap
png(file.path(OUT_DIR, "P0-1_specificity_corheatmap.png"), width = 2600, height = 2200, res = 300)
pheatmap(cor_mat,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = 10, fontsize = 11,
         color = colorRampPalette(c("#2e6fde", "white", "#d64545"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "A  Score correlation matrix (TCGA-HNSC, n=501)")
dev.off()
cat("\nSaved heatmap: P0-1_specificity_corheatmap.png\n")

## ---------- 7. Multivariable Cox: OS ~ Mitoxy + OXPHOS + Apop + Glycol + Prolif + clin ----------
df$gender2 <- factor(df$gender)
df$stage2 <- ifelse(grepl("III|IV", df$stage), "III/IV",
             ifelse(grepl("I|II", df$stage), "I/II", NA))

run_cox <- function(fit, label) {
  s <- summary(fit)
  co <- coef(fit); ci <- exp(confint(fit))
  cat(sprintf("%-42s HR=%.3f (%.3f-%.3f) P=%.4g\n", label,
              exp(co["scale(mitoxy67)"]), ci["scale(mitoxy67)", 1],
              ci["scale(mitoxy67)", 2], s$coefficients["scale(mitoxy67)", 5]))
  invisible(data.frame(Model = label,
    HR = exp(co["scale(mitoxy67)"]), LCL = ci["scale(mitoxy67)", 1],
    UCL = ci["scale(mitoxy67)", 2], P = s$coefficients["scale(mitoxy67)", 5]))
}

cat("\n===== Mitoxy67 HR across models (per 1 SD) =====\n")
res <- data.frame()
res <- rbind(res, run_cox(coxph(Surv(OS_time, OS_status) ~ scale(mitoxy67), data = df),
                          "1. Unadjusted"))
res <- rbind(res, run_cox(coxph(Surv(OS_time, OS_status) ~ scale(mitoxy67) + OXPHOS + Glycolysis + Apoptosis + Proliferation,
                          data = df), "2. + OXPHOS+Glyco+Apop+Prolif"))
res <- rbind(res, run_cox(coxph(Surv(OS_time, OS_status) ~ scale(mitoxy67) + OXPHOS + Glycolysis + Apoptosis + Proliferation + Hypoxia + ROS,
                          data = df), "3. + Hypoxia+ROS"))
res <- rbind(res, run_cox(coxph(Surv(OS_time, OS_status) ~ scale(mitoxy67) + OXPHOS + Glycolysis + Apoptosis + Proliferation + Hypoxia + ROS + Mitophagy + MitoDynamics,
                          data = df), "4. + Mitophagy+MitoDyn"))
df_c <- df[!is.na(df$stage2), ]
res <- rbind(res, run_cox(coxph(Surv(OS_time, OS_status) ~ scale(mitoxy67) + OXPHOS + Glycolysis + Apoptosis + Proliferation + Hypoxia + ROS + Mitophagy + MitoDynamics + age + gender2 + stage2 + site,
                          data = df_c), "5. + clinical (age,sex,stage,site)"))

write.csv(res, file.path(OUT_DIR, "P0-1_specificity_cox.csv"), row.names = FALSE)
cat("\n===== Final specificity Cox table =====\n"); print(res, row.names = FALSE)

## ---------- 7b. VIF / collinearity diagnostic (Model 2 structure) ----------
## Reviewer 5.2 requirement: demonstrate that HR attenuation reflects
## collinearity (variance inflation) rather than disappearance of signal.
vif_diag <- function(df, vars) {
  ## linear regression of each predictor on the others -> R^2 -> VIF
  v <- sapply(vars, function(v) {
    others <- setdiff(vars, v)
    f <- as.formula(paste(v, "~", paste(others, collapse = " + ")))
    r2 <- summary(lm(f, data = df))$r.squared
    c(VIF = 1 / (1 - r2), Tolerance = 1 - r2)
  })
  t(v)
}
cat("\n===== VIF for Model 2 predictors (mitoxy67 + 4 hallmark programs) =====\n")
m2_vars <- c("mitoxy67", "OXPHOS", "Glycolysis", "Apoptosis", "Proliferation")
vif2 <- vif_diag(df, m2_vars)
print(round(vif2, 2))
write.csv(data.frame(Variable = rownames(vif2), round(vif2, 3)),
          file.path(OUT_DIR, "P0-1_specificity_vif.csv"), row.names = FALSE)

cat("\n===== VIF for Model 4 predictors (all 10) =====\n")
m4_vars <- c("mitoxy67", "OXPHOS", "Glycolysis", "Apoptosis", "Proliferation",
             "Hypoxia", "ROS", "Mitophagy", "MitoDynamics")
vif4 <- vif_diag(df, m4_vars)
print(round(vif4, 2))
write.csv(data.frame(Variable = rownames(vif4), round(vif4, 3)),
          file.path(OUT_DIR, "P0-1_specificity_vif_model4.csv"), row.names = FALSE)

## ---------- 8. Individual module HR after full adjustment (diagnostic) ----------
cat("\n===== Module scores after full adjustment (model 5 structure) =====\n")
mods <- list(Core = intersect(c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1"), common67),
             mTOR = intersect(c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD"), common67),
             MitoDynamics = intersect(c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS"), common67),
             MetaboImmune = intersect(c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG"), common67))
for (m in names(mods)) {
  sc <- calc_score(mods[[m]], tumor_mat, df$sample)
  df_c$tmp <- sc[df_c$sample]
  fit <- coxph(Surv(OS_time, OS_status) ~ scale(tmp) + OXPHOS + Glycolysis + Apoptosis + Proliferation +
                 Hypoxia + ROS + Mitophagy + MitoDynamics + age + gender2 + stage2 + site, data = df_c)
  s <- summary(fit); ci <- exp(confint(fit))
  cat(sprintf("  %-14s HR=%.3f (%.3f-%.3f) P=%.4g\n", m,
              exp(coef(fit)["scale(tmp)"]), ci["scale(tmp)", 1], ci["scale(tmp)", 2],
              s$coefficients["scale(tmp)", 5]))
}

## ---------- 9. CD274/PD-L1 partial (purity adjustment placeholder) ----------
## Full purity adjustment in P0-3; here report unadjusted correlation for context
if ("CD274" %in% rownames(tumor_mat)) {
  cd274 <- as.numeric(log2(tumor_mat["CD274", df$sample] + 1))
  cat("\nCD274 (PD-L1) vs mitoxy67: r =",
      round(cor(df$mitoxy67, cd274, method = "spearman", use = "complete.obs"), 3), "\n")
}

cat("\n=== P0-1 specificity DONE ===")
