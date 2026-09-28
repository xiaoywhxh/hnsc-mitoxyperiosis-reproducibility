## ============================================================
## P0-1: TCGA 67-gene common program 重算（审稿意见 4 编辑必做 1）
## TCGA 73g vs 67g coverage sensitivity；Figure 2/Table 1 统一 67g estimand
## ============================================================
suppressMessages({library(survival); library(timeROC)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

common67 <- read.csv(file.path(OUT_DIR, "program_common_genes_3OS.csv"))$common_gene
cat("common 67 genes:", length(common67), "\n")

tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]

## 73-gene list（全模块）
core <- c('PRKN','VDAC1','VDAC2','VDAC3','BAX','BAK1','BID','BBC3','PMAIP1','BCL2','BCL2L1','MCL1','TSPO','AIFM1')
mtor <- c('MTOR','RICTOR','RPTOR','MLST8','MAPKAP1','PRR5','PRR5L','DEPTOR','TSC1','TSC2','RHEB','AKT1','AKT2','AKT3','PTEN','PIK3CA','PIK3CB','PIK3CD','RRAGA','RRAGB','RRAGC','RRAGD')
mito <- c('DNM1L','FIS1','MFF','MIEF1','MIEF2','MFN1','MFN2','OPA1','MARCHF5','PINK1','SLC25A3','SLC25A5','SLC25A6','CYCS')
metab <- c('HK1','HK2','PFKFB3','PKM','LDHA','LDHB','IDH1','IDH2','MDH1','MDH2','CS','GLS','GLS2','GOT1','GOT2','NLRP3','CASP1','CASP4','CASP5','STING1','CGAS','MAOB','ENDOG')
all73 <- unique(c(core, mtor, mito, metab))

## 计算 73g 和 67g program score（program-mean Z，同一评分规则）
calc_prog <- function(genes) {
  ex <- log2(tumor_mat[genes, df$sample] + 1)
  z <- t(scale(t(ex)))
  colMeans(z[intersect(genes, rownames(z)), , drop = FALSE], na.rm = TRUE)
}
prog73 <- calc_prog(all73)
prog67 <- calc_prog(common67)

fit73 <- coxph(Surv(OS_time, OS_status) ~ scale(prog73), data = df)
fit67 <- coxph(Surv(OS_time, OS_status) ~ scale(prog67), data = df)
s73 <- summary(fit73); s67 <- summary(fit67)
ci73 <- exp(confint(fit73)); ci67 <- exp(confint(fit67))

cat("\n=== TCGA 73g vs 67g program（coverage sensitivity）===\n")
cat(sprintf("73g: HR=%.3f (%.3f-%.3f) P=%.4f\n", exp(coef(fit73)), ci73[1], ci73[2], s73$coefficients[5]))
cat(sprintf("67g: HR=%.3f (%.3f-%.3f) P=%.4f\n", exp(coef(fit67)), ci67[1], ci67[2], s67$coefficients[5]))

## 模块级 67g
mods <- list(Core = intersect(core, common67), mTOR = intersect(mtor, common67),
             MitoDynamics = intersect(mito, common67), MetaboImmune = intersect(metab, common67))
cat("\n67g 模块基因数:", sapply(mods, length), "\n")
for (m in names(mods)) {
  sc <- calc_prog(mods[[m]])
  fit <- coxph(Surv(OS_time, OS_status) ~ scale(sc), data = df)
  cat(sprintf("  67g %s: HR=%.3f P=%.4f\n", m, exp(coef(fit)), summary(fit)$coefficients[5]))
}

## 时间依赖 AUC（67g）
roc67 <- timeROC(T = df$OS_time, delta = df$OS_status, marker = prog67, cause = 1,
                 times = c(365, 1095, 1825), iid = FALSE)
cat("\n67g AUC 1/3/5yr:", round(roc67$AUC, 3), "\n")

## 保存
res <- data.frame(
  Program = c("TCGA 73-gene", "TCGA 67-gene"),
  Genes = c(length(all73), length(common67)),
  HR = c(exp(coef(fit73)), exp(coef(fit67))),
  LCL = c(ci73[1], ci67[1]), UCL = c(ci73[2], ci67[2]),
  P = c(s73$coefficients[5], s67$coefficients[5]),
  N = nrow(df), Events = sum(df$OS_status)
)
write.csv(res, file.path(OUT_DIR, "TCGA_73g_vs_67g_program.csv"), row.names = FALSE)
cat("\n已保存 TCGA_73g_vs_67g_program.csv\n")
print(res, digits = 4, row.names = FALSE)
cat("\nDONE\n")
