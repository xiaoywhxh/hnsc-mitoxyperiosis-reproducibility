## ============================================================
## P1-2 CANONICAL (2026-09-27) — supersedes V5_P1_2_oral_external.R (LEGACY)
## 口径 = Canonical scoring specification v1.0：
##   固定 67-gene + TCGA log2(x+1) + GEO 直接 + 基因 z 在「父队列」内 + 子集内 per-1-SD Cox
## 与旧实现的差异：
##   (1) 旧版对 TCGA 原始计数矩阵漏做 log2(x+1)  -> 主因（TCGA oral 1.312 -> 1.3407）
##   (2) 旧版 subgroup 基因 z 用子集而非父队列参照 -> GSE65858 oral 0.6693 -> 0.6800
##   (3) 旧版未对 site 未映射(NA)样本做 ghost-row 防护 -> 统一 which()
## 输出：oral_external_comparison_67g_CANONICAL.csv + FigureS13 源数据（与 corrected_S20 一致）
## ============================================================
suppressMessages({library(survival)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

common67 <- read.csv(file.path(OUT_DIR, "program_common_genes_3OS.csv"))$common_gene

## 规范评分：gene z 在「父队列」内，均值得分；Cox 在子集内 per-1-SD
prog_cox_parentz <- function(expr, samples_parent, samples_sub, time, status, do_log = FALSE) {
  ex <- expr[intersect(common67, rownames(expr)), samples_parent, drop = FALSE]
  if (do_log) ex <- log2(ex + 1)
  z <- t(scale(t(ex)))
  sc_parent <- colMeans(z, na.rm = TRUE)
  sc_sub <- sc_parent[samples_sub]
  fit <- coxph(Surv(time, status) ~ scale(sc_sub))
  s <- summary(fit); ci <- exp(confint(fit))
  list(HR = exp(coef(fit)), LCL = ci[1], UCL = ci[2], P = s$coefficients[5],
       N = length(time), E = sum(status))
}

res <- data.frame()

## 1. TCGA overall（父队列 = 501 全队）
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

r0 <- prog_cox_parentz(tumor_mat, df$sample, df$sample, df$OS_time, df$OS_status, do_log = TRUE)
res <- rbind(res, data.frame(Cohort = "TCGA overall (reference)", N = r0$N, Events = r0$E,
                             HR = r0$HR, LCL = r0$LCL, UCL = r0$UCL, P = r0$P))

## 2. TCGA oral subset（父队列 z = 501，子集 Cox = oral 309）
oral_idx <- which(df$site == "Oral cavity")
r <- prog_cox_parentz(tumor_mat, df$sample, df$sample[oral_idx],
                      df$OS_time[oral_idx], df$OS_status[oral_idx], do_log = TRUE)
res <- rbind(res, data.frame(Cohort = "TCGA oral subset", N = r$N, Events = r$E,
                             HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))

## 3. GSE41613（100% oral；GEO 已 log，父=子）
expr416 <- readRDS(file.path(DATA_DIR, "GSE41613_expr_full.Rds"))
surv416 <- readRDS(file.path(DATA_DIR, "GSE41613_surv.Rds"))
s416 <- colnames(expr416)
r <- prog_cox_parentz(expr416, s416, s416, surv416$OS_time, surv416$OS_status, do_log = FALSE)
res <- rbind(res, data.frame(Cohort = "GSE41613 (100% oral)", N = r$N, Events = r$E,
                             HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))

## 4. GSE42743 overall + oral subset（父队列 z = 74）
expr427 <- readRDS(file.path(DATA_DIR, "GSE42743_processed/GSE42743_expression_log2RMA_gene.Rds"))
clin427 <- read.csv(file.path(DATA_DIR, "GSE42743_processed/GSE42743_clinical.csv"), check.names = FALSE)
common427 <- intersect(colnames(expr427), clin427$sample_id)
clin427 <- clin427[match(common427, clin427$sample_id), ]
r <- prog_cox_parentz(expr427, common427, common427, clin427$OS_time, clin427$OS_status, do_log = FALSE)
res <- rbind(res, data.frame(Cohort = "GSE42743 overall (96% oral)", N = r$N, Events = r$E,
                             HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))
o427 <- which(clin427$tumor_site == "Oral cavity")
r <- prog_cox_parentz(expr427, common427, common427[o427],
                      clin427$OS_time[o427], clin427$OS_status[o427], do_log = FALSE)
res <- rbind(res, data.frame(Cohort = "GSE42743 oral subset", N = r$N, Events = r$E,
                             HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))

## 5. GSE65858 oral subset（父队列 z = 全部表达样本；ghost-row 用 which()）
pd <- readRDS(file.path(DATA_DIR, "GSE65858_pdata.Rds"))
expr658 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))
pd$site_group <- ifelse(grepl("Cavum Oris|Oral", pd[["tumor_site:ch1"]]), "Oral cavity",
                 ifelse(grepl("Oropharynx", pd[["tumor_site:ch1"]]), "Oropharynx",
                 ifelse(grepl("Larynx", pd[["tumor_site:ch1"]]), "Larynx",
                 ifelse(grepl("Hypopharynx", pd[["tumor_site:ch1"]]), "Hypopharynx", NA))))
pd$os_time <- as.numeric(pd[["os:ch1"]])
pd$os_event <- as.integer(pd[["os_event:ch1"]] == "TRUE")
common658 <- intersect(colnames(expr658), pd$geo_accession)
oral658 <- which(pd$geo_accession %in% common658 & pd$site_group == "Oral cavity" & !is.na(pd$os_event))
r <- prog_cox_parentz(expr658, common658, pd$geo_accession[oral658],
                      pd$os_time[oral658], pd$os_event[oral658], do_log = FALSE)
res <- rbind(res, data.frame(Cohort = "GSE65858 oral subset", N = r$N, Events = r$E,
                             HR = r$HR, LCL = r$LCL, UCL = r$UCL, P = r$P))

write.csv(res, file.path(OUT_DIR, "oral_external_comparison_67g_CANONICAL.csv"), row.names = FALSE)
write.csv(res, file.path(OUT_DIR, "FigureS13_oral_external_data_CANONICAL.csv"), row.names = FALSE)
cat("\n=== CANONICAL 结果 ===\n")
print(res, digits = 6, row.names = FALSE)
cat("\n已保存 oral_external_comparison_67g_CANONICAL.csv + FigureS13_oral_external_data_CANONICAL.csv\nDONE\n")
