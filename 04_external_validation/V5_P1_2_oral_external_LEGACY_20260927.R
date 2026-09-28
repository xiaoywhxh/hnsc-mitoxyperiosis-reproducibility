## LEGACY — DO NOT USE FOR MANUSCRIPT RESULTS. SUPERSEDED BY V5_P1_2_oral_external_CANONICAL_20260927.R (provenance only).
## ⚠️ LEGACY — 已被 CANONICAL 版本取代（2026-09-27）
## 取代原因：本实现对 TCGA 原始计数矩阵漏做预规 log2(x+1) 变换，且 subgroup 基因 z 以子集而非父队列为参照，
##           导致 S20/Figure S13 的 TCGA HR 系统性偏低（TCGA oral 1.312 vs canonical 1.3407）。
## supersedes 见 V5_P1_2_oral_external_CANONICAL_20260927.R；本文件仅作 provenance 存档，禁止再运行。
## Canonical scoring specification v1.0：见 Canonical_scoring_specification_v1.0.md

## ============================================================
## P1-2: 口腔外部比较（审稿意见 4 P1-2 + Reviewer2 P1-1）
## TCGA oral / GSE41613 (100% oral) / GSE42743 (96% oral) / GSE65858 oral
## 统一 67-gene common program, program-mean Z, per-1-SD Cox
## ============================================================
suppressMessages({library(survival)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

common67 <- read.csv(file.path(OUT_DIR, "program_common_genes_3OS.csv"))$common_gene

prog_cox <- function(expr, time, status, genes = common67) {
  sub <- expr[intersect(genes, rownames(expr)), , drop = FALSE]
  z <- t(scale(t(sub)))
  sc <- colMeans(z, na.rm = TRUE)
  fit <- coxph(Surv(time, status) ~ scale(sc))
  s <- summary(fit); ci <- exp(confint(fit))
  list(HR = exp(coef(fit)), LCL = ci[1], UCL = ci[2], P = s$coefficients[5], N = length(time), E = sum(status))
}

res <- data.frame()

## 1. TCGA oral subset（67g）
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
oral <- df[df$site == "Oral cavity", ]
r <- prog_cox(tumor_mat[, oral$sample], oral$OS_time, oral$OS_status)
res <- rbind(res, data.frame(Cohort = "TCGA oral subset", N = r$N, Events = r$E, HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))
cat(sprintf("TCGA oral (67g): N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4f\n", r$N, r$E, r$HR, r$LCL, r$UCL, r$P))

## 2. GSE41613 overall（100% oral）
expr416 <- readRDS(file.path(DATA_DIR, "GSE41613_expr_full.Rds"))
surv416 <- readRDS(file.path(DATA_DIR, "GSE41613_surv.Rds"))
r <- prog_cox(expr416, surv416$OS_time, surv416$OS_status)
res <- rbind(res, data.frame(Cohort = "GSE41613 (100% oral)", N = r$N, Events = r$E, HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))
cat(sprintf("GSE41613 (67g): N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4f\n", r$N, r$E, r$HR, r$LCL, r$UCL, r$P))

## 3. GSE42743（96% oral: overall + oral subset）
expr427 <- readRDS(file.path(DATA_DIR, "GSE42743_processed/GSE42743_expression_log2RMA_gene.Rds"))
clin427 <- read.csv(file.path(DATA_DIR, "GSE42743_processed/GSE42743_clinical.csv"), check.names = FALSE)
common <- intersect(colnames(expr427), clin427$sample_id)
clin427 <- clin427[match(common, clin427$sample_id), ]
r <- prog_cox(expr427[, common], clin427$OS_time, clin427$OS_status)
res <- rbind(res, data.frame(Cohort = "GSE42743 overall (96% oral)", N = r$N, Events = r$E, HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))
cat(sprintf("GSE42743 overall (67g): N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4f\n", r$N, r$E, r$HR, r$LCL, r$UCL, r$P))
oral427 <- clin427$tumor_site == "Oral cavity"
r <- prog_cox(expr427[, common[oral427]], clin427$OS_time[oral427], clin427$OS_status[oral427])
res <- rbind(res, data.frame(Cohort = "GSE42743 oral subset", N = r$N, Events = r$E, HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))
cat(sprintf("GSE42743 oral subset (67g): N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4f\n", r$N, r$E, r$HR, r$LCL, r$UCL, r$P))

## 4. GSE65858 oral subset
pd <- readRDS(file.path(DATA_DIR, "GSE65858_pdata.Rds"))
expr658 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))
pd$site <- ifelse(grepl("Cavum Oris|Oral", pd[["tumor_site:ch1"]]), "Oral cavity",
           ifelse(grepl("Oropharynx", pd[["tumor_site:ch1"]]), "Oropharynx",
           ifelse(grepl("Larynx", pd[["tumor_site:ch1"]]), "Larynx",
           ifelse(grepl("Hypopharynx", pd[["tumor_site:ch1"]]), "Hypopharynx", NA))))
pd$os_time <- as.numeric(pd[["os:ch1"]])
pd$os_event <- as.integer(pd[["os_event:ch1"]] == "TRUE")
oral658 <- which(pd$site == "Oral cavity" & !is.na(pd$os_event))
r <- prog_cox(expr658[, pd$geo_accession[oral658]], pd$os_time[oral658], pd$os_event[oral658])
res <- rbind(res, data.frame(Cohort = "GSE65858 oral subset", N = r$N, Events = r$E, HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))
cat(sprintf("GSE65858 oral (67g): N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4f\n", r$N, r$E, r$HR, r$LCL, r$UCL, r$P))

## 全队列对比（TCGA overall 67g 作参考）
df_all <- df
r0 <- prog_cox(tumor_mat[, df_all$sample], df_all$OS_time, df_all$OS_status)
res <- rbind(data.frame(Cohort = "TCGA overall (reference)", N = r0$N, Events = r0$E, HR = r0$HR, LCL = r0$LCL, UCL = r0$UCL, P = r0$P), res)
cat(sprintf("TCGA overall (67g ref): N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4f\n", r0$N, r0$E, r0$HR, r0$LCL, r0$UCL, r0$P))

write.csv(res, file.path(OUT_DIR, "oral_external_comparison_67g.csv"), row.names = FALSE)
cat("\n已保存 oral_external_comparison_67g.csv\n")
print(res, digits = 4, row.names = FALSE)
cat("\nDONE\n")
