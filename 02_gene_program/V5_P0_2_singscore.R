## ============================================================
## P0-2: 标准 singscore rank-based sensitivity（审稿意见 4 编辑必做 2）
## 全转录组 rankGenes + simpleScore（官方实现），重算 5 个模块
## 诊断 previous manual "overall constant" 问题
## ============================================================
suppressMessages({library(singscore); library(survival)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/02_gene_program"

core <- c('PRKN','VDAC1','VDAC2','VDAC3','BAX','BAK1','BID','BBC3','PMAIP1','BCL2','BCL2L1','MCL1','TSPO','AIFM1')
mtor <- c('MTOR','RICTOR','RPTOR','MLST8','MAPKAP1','PRR5','PRR5L','DEPTOR','TSC1','TSC2','RHEB','AKT1','AKT2','AKT3','PTEN','PIK3CA','PIK3CB','PIK3CD','RRAGA','RRAGB','RRAGC','RRAGD')
mito <- c('DNM1L','FIS1','MFF','MIEF1','MIEF2','MFN1','MFN2','OPA1','MARCHF5','PINK1','SLC25A3','SLC25A5','SLC25A6','CYCS')
metab <- c('HK1','HK2','PFKFB3','PKM','LDHA','LDHB','IDH1','IDH2','MDH1','MDH2','CS','GLS','GLS2','GOT1','GOT2','NLRP3','CASP1','CASP4','CASP5','STING1','CGAS','MAOB','ENDOG')
all73 <- unique(c(core, mtor, mito, metab))
mods <- list(Core = core, mTOR = mtor, MitoDynamics = mito, MetaboImmune = metab, All = all73)

## 数据
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]

## 全转录组 log2(TPM+1)，去除零方差基因（singscore 要求）
expr <- log2(tumor_mat[, df$sample] + 1)
sdv <- apply(expr, 1, sd, na.rm = TRUE)
keep <- rownames(expr)[sdv > 0 & is.finite(sdv)]
expr <- as.matrix(expr[keep, ])
cat("全转录组用于 rankGenes:", nrow(expr), "genes x", ncol(expr), "samples\n")

## rankGenes（每样本内全转录组 rank；singscore 官方：rank 相对整个可测转录组）
ranked <- rankGenes(expr)
cat("ranked dims:", dim(ranked), " range:", range(ranked), "\n")

## simpleScore for each module（仅 upSet，singscore 官方签名）
res <- data.frame()
for (m in names(mods)) {
  gs <- mods[[m]]
  gs_present <- intersect(gs, rownames(expr))
  if (length(gs_present) < 3) next
  sc <- simpleScore(ranked, upSet = gs_present)
  score <- sc$TotalScore
  names(score) <- rownames(sc)
  # 匹配样本
  score <- score[df$sample]
  cat(sprintf("%s: %d genes, score range %.4f-%.4f, SD %.4f, constant=%s\n",
      m, length(gs_present), min(score), max(score), sd(score), sd(score) < 1e-10))
  fit <- coxph(Surv(df$OS_time, df$OS_status) ~ scale(score))
  s <- summary(fit); ci <- exp(confint(fit))
  res <- rbind(res, data.frame(Module = m, Genes = length(gs_present),
    HR = exp(coef(fit)), LCL = ci[1], UCL = ci[2], P = s$coefficients[5],
    Score_SD = sd(score)))
  cat(sprintf("  HR=%.3f (%.3f-%.3f) P=%.4f\n", exp(coef(fit)), ci[1], ci[2], s$coefficients[5]))
}
res$P_q <- p.adjust(res$P, method = "BH")
write.csv(res, file.path(OUT_DIR, "singscore_module_survival_501.csv"), row.names = FALSE)

## mean-Z 对照
z_all <- t(scale(t(log2(tumor_mat[intersect(all73, rownames(tumor_mat)), df$sample] + 1))))
meanZ <- function(gs) colMeans(z_all[intersect(gs, rownames(z_all)), , drop = FALSE], na.rm = TRUE)
cat("\n=== mean-Z 对照 ===\n")
for (m in names(mods)) {
  sc <- meanZ(mods[[m]])
  fit <- coxph(Surv(df$OS_time, df$OS_status) ~ scale(sc))
  cat(sprintf("%s: HR=%.3f P=%.4f\n", m, exp(coef(fit)), summary(fit)$coefficients[5]))
}

## 关键诊断：overall 73g rank score 是否恒定（审稿人核心疑虑）
sc_all <- simpleScore(ranked, upSet = intersect(all73, rownames(expr)))
ov <- sc_all$TotalScore
cat(sprintf("\n[诊断] overall 73g rank score: SD=%.6f, range %.4f-%.4f, constant=%s\n",
    sd(ov), min(ov), max(ov), sd(ov) < 1e-10))
cat("==> 官方 singscore 中 overall 73g 不再恒定（rank 相对全转录组）\n")

cat("\nDONE\n")
