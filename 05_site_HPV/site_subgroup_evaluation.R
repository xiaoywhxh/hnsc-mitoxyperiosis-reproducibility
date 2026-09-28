## ============================================================
## 亚组分析可行性评估: site-stratified program Cox
## TCGA (501) + GSE65858 (266) 内部亚组
## ============================================================
suppressMessages({library(survival)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

core <- c('PRKN','VDAC1','VDAC2','VDAC3','BAX','BAK1','BID','BBC3','PMAIP1','BCL2','BCL2L1','MCL1','TSPO','AIFM1')
mtor <- c('MTOR','RICTOR','RPTOR','MLST8','MAPKAP1','PRR5','PRR5L','DEPTOR','TSC1','TSC2','RHEB','AKT1','AKT2','AKT3','PTEN','PIK3CA','PIK3CB','PIK3CD','RRAGA','RRAGB','RRAGC','RRAGD')
mito <- c('DNM1L','FIS1','MFF','MIEF1','MIEF2','MFN1','MFN2','OPA1','MARCHF5','PINK1','SLC25A3','SLC25A5','SLC25A6','CYCS')
metab <- c('HK1','HK2','PFKFB3','PKM','LDHA','LDHB','IDH1','IDH2','MDH1','MDH2','CS','GLS','GLS2','GOT1','GOT2','NLRP3','CASP1','CASP4','CASP5','STING1','CGAS','MAOB','ENDOG')
all73 <- unique(c(core, mtor, mito, metab))

## ---------- 1. TCGA site-stratified ----------
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]

site_os <- read.csv(file.path(OUT_DIR, "site_TCGA_501_os.csv"))
df$site <- site_os$site_group[match(substr(df$sample, 1, 12), site_os$submitter_id)]

expr <- log2(tumor_mat[all73, df$sample] + 1)
z <- t(scale(t(expr)))
score_fn <- function(gs) colMeans(z[intersect(gs, rownames(z)), , drop = FALSE], na.rm = TRUE)
df$prog_all <- score_fn(all73)
df$prog_core <- score_fn(core)
df$prog_mito <- score_fn(mito)
df$prog_metab <- score_fn(metab)

cat("===== TCGA site-stratified program Cox (per-1-SD) =====\n")
res_tcga <- data.frame()
for (st in c("Oral cavity", "Oropharynx", "Larynx")) {
  sub <- df[df$site == st, ]
  if (nrow(sub) < 30) next
  fit <- coxph(Surv(OS_time, OS_status) ~ scale(prog_all), data = sub)
  s <- summary(fit)
  res_tcga <- rbind(res_tcga, data.frame(
    Cohort = "TCGA", Site = st, N = nrow(sub), Events = sum(sub$OS_status),
    HR = exp(coef(fit)), LCL = exp(confint(fit))[1], UCL = exp(confint(fit))[2],
    P = s$coefficients[5]))
  cat(sprintf("  %-12s N=%-4d E=%-4d HR=%.3f (%.3f-%.3f) P=%.4f\n",
      st, nrow(sub), sum(sub$OS_status), exp(coef(fit)),
      exp(confint(fit))[1], exp(confint(fit))[2], s$coefficients[5]))
}

# 模块级（Oral cavity 探索）
cat("\n-- TCGA Oral cavity 模块级 --\n")
for (m in c("prog_core", "prog_mito", "prog_metab")) {
  sub <- df[df$site == "Oral cavity", ]
  fit <- coxph(Surv(OS_time, OS_status) ~ scale(sub[[m]]), data = sub)
  s <- summary(fit)
  cat(sprintf("  %-10s HR=%.3f (%.3f-%.3f) P=%.4f\n", m, exp(coef(fit)),
      exp(confint(fit))[1], exp(confint(fit))[2], s$coefficients[5]))
}

## ---------- 2. GSE65858 site-stratified ----------
pd <- readRDS(file.path(DATA_DIR, "GSE65858_pdata.Rds"))
map <- function(x) ifelse(is.na(x), "NA",
  ifelse(grepl("Cavum Oris|Oral", x), "Oral cavity",
  ifelse(grepl("Oropharynx", x), "Oropharynx",
  ifelse(grepl("Larynx", x), "Larynx",
  ifelse(grepl("Hypopharynx", x), "Hypopharynx", "Other")))))
pd$site_group <- map(pd[["tumor_site:ch1"]])
pd$os_time <- as.numeric(pd[["os:ch1"]])
pd$os_event <- as.integer(pd[["os_event:ch1"]] == "TRUE")  # "TRUE"/"FALSE" 字符串

expr658 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))
common <- intersect(colnames(expr658), pd$geo_accession)
cat("\n===== GSE65858 site-stratified program Cox (per-1-SD) =====\n")
cat("匹配样本:", length(common), "\n")
if (length(common) > 0) {
  ex <- expr658[intersect(all73, rownames(expr658)), common]
  cat("可测基因:", nrow(ex), "/", length(all73), "\n")
  z658 <- t(scale(t(ex)))
  score658 <- function(gs) colMeans(z658[intersect(gs, rownames(z658)), , drop = FALSE], na.rm = TRUE)
  prog_all <- score658(all73)
  pd2 <- pd[match(common, pd$geo_accession), ]
  pd2$prog_all <- prog_all
  res658 <- data.frame()
  for (st in c("Oral cavity", "Oropharynx", "Larynx", "Hypopharynx")) {
    sub <- pd2[pd2$site_group == st & !is.na(pd2$os_event), ]
    if (nrow(sub) < 20) next
    fit <- coxph(Surv(os_time, os_event) ~ scale(prog_all), data = sub)
    s <- summary(fit)
    res658 <- rbind(res658, data.frame(
      Cohort = "GSE65858", Site = st, N = nrow(sub), Events = sum(sub$os_event),
      HR = exp(coef(fit)), LCL = exp(confint(fit))[1], UCL = exp(confint(fit))[2],
      P = s$coefficients[5]))
    cat(sprintf("  %-12s N=%-4d E=%-4d HR=%.3f (%.3f-%.3f) P=%.4f\n",
        st, nrow(sub), sum(sub$os_event), exp(coef(fit)),
        exp(confint(fit))[1], exp(confint(fit))[2], s$coefficients[5]))
  }
  write.csv(rbind(res_tcga, res658), file.path(OUT_DIR, "site_stratified_program_cox.csv"), row.names = FALSE)
}

write.csv(res_tcga, file.path(OUT_DIR, "site_stratified_TCGA.csv"), row.names = FALSE)
cat("\nDONE\n")
