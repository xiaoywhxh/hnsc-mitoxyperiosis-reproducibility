## P0-3: mean-Z program score sensitivity analysis (审稿意见 2 P0-3 / R1-M5 / R3-M10)
## 1) module overlap matrix (82 assignments / 73 unique) 2) 去共享基因后模块相关性
## 3) direction-aware (signed) 敏感性 4) score 方法对比 (mean-Z vs ssGSEA 已有)

suppressMessages({ library(survival) })

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/02_gene_program"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 加载 501 全队列
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
sv <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
sv <- sv[!is.na(sv$OS_time), ]
cat("501 cohort:", nrow(sv), "\n")

core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
modules <- list(Core = core_genes, mTOR = mtor_genes, MitoDynamics = mito_genes,
                MetaboImmune = metab_genes)
modnames <- c("Core","mTOR","MitoDynamics","MetaboImmune")

# ============================================================
# 1. Module overlap matrix
# ============================================================
cat("\n===== 1. Module overlap =====\n")
overlap <- matrix(0, 4, 4, dimnames = list(modnames, modnames))
shared_mat <- matrix("", 4, 4, dimnames = list(modnames, modnames))
for (i in 1:4) for (j in 1:4) {
  si <- intersect(modules[[i]], modules[[j]])
  overlap[i, j] <- length(si)
  shared_mat[i, j] <- if (length(si) > 0) paste(si, collapse = ", ") else ""
}
print(overlap)
cat("\nShared gene details:\n")
print(shared_mat)
cat("\n82 assignments breakdown:\n")
print(sapply(modules, length))
cat("Unique genes:", length(unique(unlist(modules))), "\n")
write.csv(data.frame(Module = modnames, N_genes = sapply(modules, length)), 
          file.path(OUT_DIR, "module_sizes.csv"), row.names = FALSE)

# ============================================================
# 2. Expression + mean-Z scores (full vs no-shared)
# ============================================================
expr <- log2(tumor_mat[unique(unlist(modules)), sv$sample] + 1)
z <- t(scale(t(expr)))
score_full <- sapply(modnames, function(m) {
  gi <- intersect(modules[[m]], rownames(z))
  colMeans(z[gi, , drop = FALSE], na.rm = TRUE)
})

# 去共享基因：每个模块的"专属基因"（不在其他模块的并集中）
union_other <- function(m) {
  others <- setdiff(modnames, m)
  setdiff(modules[[m]], unlist(modules[others]))
}
exclusive <- lapply(modnames, union_other)
names(exclusive) <- modnames
cat("\nExclusive gene counts:\n")
print(sapply(exclusive, length))

score_excl <- sapply(modnames, function(m) {
  gi <- intersect(exclusive[[m]], rownames(z))
  if (length(gi) == 0) return(rep(NA, ncol(z)))
  colMeans(z[gi, , drop = FALSE], na.rm = TRUE)
})

# ============================================================
# 3. MetaboImmune vs mTOR correlation: full vs exclusive
# ============================================================
cat("\n===== 2. MetaboImmune vs mTOR correlation sensitivity =====\n")
r_full <- cor(score_full[, "MetaboImmune"], score_full[, "mTOR"], method = "spearman")
r_excl <- cor(score_excl[, "MetaboImmune"], score_excl[, "mTOR"], method = "spearman")
cat(sprintf("Full modules:     MetaboImmune vs mTOR Spearman r = %.3f\n", r_full))
cat(sprintf("Exclusive genes:  MetaboImmune vs mTOR Spearman r = %.3f (Metabo excl %d genes, mTOR excl %d genes)\n",
            r_excl, length(exclusive[["MetaboImmune"]]), length(exclusive[["mTOR"]])))

# ============================================================
# 4. Direction-aware (signed) score 敏感性
# ============================================================
cat("\n===== 3. Direction-aware (signed) score =====\n")
# Tier A/B/C expected direction (从 B3 证据表简化: 促死亡/促程序 = +1, 抗死亡/抑制 = -1)
# 仅对已知方向基因做 signed mean-Z
neg_genes <- c("BCL2","BCL2L1","MCL1","PTEN","TSC1","TSC2","DEPTOR")  # 抑制程序 → 取负
dir_w <- setNames(rep(1, length(unique(unlist(modules)))), unique(unlist(modules)))
dir_w[neg_genes] <- -1

score_signed <- sapply(modnames, function(m) {
  gi <- intersect(modules[[m]], rownames(z))
  colMeans(dir_w[gi] * z[gi, , drop = FALSE], na.rm = TRUE)
})
r_signed <- cor(score_signed[, "MetaboImmune"], score_signed[, "mTOR"], method = "spearman")
cat(sprintf("Signed score:     MetaboImmune vs mTOR Spearman r = %.3f\n", r_signed))

# ============================================================
# 5. Score 方法对比汇总（含已有 ssGSEA 值）
# ============================================================
sens_df <- data.frame(
  Method = c("mean-Z (full modules)", "mean-Z (exclusive genes only)",
             "mean-Z (direction-signed)", "ssGSEA (226-case, GSVA v1.50.5)"),
  MetaboImmune_vs_mTOR_r = c(r_full, r_excl, r_signed, -0.31)
)
write.csv(sens_df, file.path(OUT_DIR, "module_correlation_sensitivity.csv"), row.names = FALSE)
print(sens_df)

cat("\n=== P0-3 DONE ===")
