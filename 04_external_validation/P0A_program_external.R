## P0-A: Program-level external survival analysis (审稿意见 3 核心 P0)
## 目标: 真正回答标题问题 "program transportability" — 用 73-gene overall program（非 11-gene LASSO）
## 在外部队列验证: common-gene program → 同 program-mean Z scoring → per-1-SD Cox → OS
## 1) 73 基因在 GSE41613/GSE42743/GSE65858 的覆盖 2) common intersection
## 3) 每队列 overall program HR 4) 与 TCGA HR=1.18 对比（program vs program, 同 estimand）

suppressMessages({ library(survival) })

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
all73 <- unique(c(core_genes, mtor_genes, mito_genes, metab_genes))

# ============================================================
# 1. 各队列 73 基因覆盖
# ============================================================
# GSE41613 (complete series matrix rebuild: 73/73; PRKN direct, no PARK2 mapping needed)
expr41613 <- readRDS(file.path(DATA_DIR, "GSE41613_expr_full.Rds"))
surv41613 <- readRDS(file.path(DATA_DIR, "GSE41613_surv.Rds"))
rn41613 <- rownames(expr41613)
cov41613 <- setNames(rep(FALSE, 73), all73)
for (g in all73) {
  cov41613[g] <- (g %in% rn41613) || (g == "PRKN" && "PARK2" %in% rn41613)
}

# GSE42743
expr42743 <- readRDS(file.path(DATA_DIR, "GSE42743_processed/GSE42743_expression_log2RMA_gene.Rds"))
clin42743 <- read.csv(file.path(DATA_DIR, "GSE42743_processed/GSE42743_clinical.csv"), check.names = FALSE)
rn42743 <- rownames(expr42743)
cov42743 <- setNames(all73 %in% rn42743, all73)

# GSE65858
expr65858 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))
surv65858 <- readRDS(file.path(DATA_DIR, "GSE65858_surv.Rds"))
rn65858 <- rownames(expr65858)
cov65858 <- setNames(all73 %in% rn65858, all73)

cat("73-gene coverage:\n")
cat(sprintf("  GSE41613: %d/73\n", sum(cov41613)))
cat(sprintf("  GSE42743: %d/73\n", sum(cov42743)))
cat(sprintf("  GSE65858: %d/73\n", sum(cov65858)))
common3 <- names(which(cov41613 & cov42743 & cov65858))
cat("3 OS 队列共同可测 (common):", length(common3), "genes\n")
cat("  缺失 per 队列: GSE41613", paste(all73[!cov41613], collapse = ","),
    "| GSE42743", paste(all73[!cov42743], collapse = ","),
    "| GSE65858", paste(all73[!cov65858], collapse = ","), "\n")

# GSE27020 (DFS secondary) — probe-level 提取后检查
lines27020 <- readLines(gzfile(file.path(DATA_DIR, "GSE27020_series_matrix.txt.gz"), open = "rt"), warn = FALSE)
start27020 <- grep("!series_matrix_table_begin", lines27020)

# ============================================================
# 2. Program score per cohort (common-gene program-mean Z)
# ============================================================
prog_cox <- function(expr_sub, time, status, common) {
  gi <- intersect(common, rownames(expr_sub))
  z <- t(scale(t(expr_sub[gi, , drop = FALSE])))
  prog <- colMeans(z, na.rm = TRUE)
  fit <- coxph(Surv(time, status) ~ scale(prog))
  ss <- summary(fit)
  c(HR = exp(coef(fit)), LCL = exp(confint(fit))[1], UCL = exp(confint(fit))[2],
    P = ss$coefficients[5], N = length(time), Events = sum(status), Genes = length(gi))
}

# GSE41613 (PRKN→PARK2)
expr41613_use <- expr41613
rn_map <- rn41613
# 若 PRKN 缺失用 PARK2 行重命名
if (!"PRKN" %in% rn_map && "PARK2" %in% rn_map) {
  rownames(expr41613_use)[rownames(expr41613_use) == "PARK2"] <- "PRKN"
}
s416 <- surv41613[!is.na(surv41613$OS_time_months), ]
e416 <- expr41613_use[, s416$sample]
r416 <- prog_cox(e416, s416$OS_time_months, s416$OS_status, common3)
cat(sprintf("\nGSE41613 program: HR=%.3f (%.3f-%.3f) P=%.4f N=%d E=%d (common %d genes)\n",
            r416[1], r416[2], r416[3], r416[4], r416[5], r416[6], r416[7]))

# GSE42743
c427 <- clin42743[!is.na(clin42743$os_time), ]
e427 <- expr42743[, as.character(c427$sample_id)]
r427 <- prog_cox(e427, c427$os_time, c427$os_status, common3)
cat(sprintf("GSE42743 program: HR=%.3f (%.3f-%.3f) P=%.4f N=%d E=%d (common %d genes)\n",
            r427[1], r427[2], r427[3], r427[4], r427[5], r427[6], r427[7]))

# GSE65858
s658 <- surv65858[!is.na(surv65858$OS_time_days), ]
e658 <- expr65858[, s658$sample]
r658 <- prog_cox(e658, s658$OS_time_days, s658$OS_status, common3)
cat(sprintf("GSE65858 program: HR=%.3f (%.3f-%.3f) P=%.4f N=%d E=%d (common %d genes)\n",
            r658[1], r658[2], r658[3], r658[4], r658[5], r658[6], r658[7]))

# GSE27020 (DFS, secondary) — probe 提取 73 基因工作量太大，用 locked 队列已有的 9 基因子集近似？
# 暂不做 GSE27020 program（DFS secondary；正文注明）

# ============================================================
# 3. 输出
# ============================================================
res <- data.frame(
  Cohort = c("TCGA (discovery)", "GSE41613", "GSE42743", "GSE65858"),
  Endpoint = c("OS", "OS", "OS", "OS"),
  N = c(501, r416[5], r427[5], r658[5]),
  Events = c(218, r416[6], r427[6], r658[6]),
  Genes_used = c(73, r416[7], r427[7], r658[7]),
  HR_per_SD = c(1.176, r416[1], r427[1], r658[1]),
  LCL = c(1.022, r416[2], r427[2], r658[2]),
  UCL = c(1.354, r416[3], r427[3], r658[3]),
  P = c(0.0238, r416[4], r427[4], r658[4])
)
write.csv(res, file.path(OUT_DIR, "program_external_common73.csv"), row.names = FALSE)
print(res, digits = 4)

# 保存 common genes
write.csv(data.frame(common_gene = common3), file.path(OUT_DIR, "program_common_genes_3OS.csv"), row.names = FALSE)

cat("\n=== P0-A DONE ===")
