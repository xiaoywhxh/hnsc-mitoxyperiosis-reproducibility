## 跨队列 HPV+ 分组比较（67-gene common program）
## 2026-08-31 | V5 配套
suppressMessages({library(survival)})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

# ============================================================
# 0. 67-gene common program（从 P0A 输出读取）
# ============================================================
common_genes <- read.csv(file.path(OUT_DIR, "program_common_genes_3OS.csv"))
genes67 <- common_genes$common_gene
cat("67-gene common program:", length(genes67), "genes\n")

# program-mean Z score 函数
prog_score <- function(expr) {
  # expr: genes x samples (log 尺度), 返回每样本 mean-Z
  z <- t(scale(t(expr)))
  colMeans(z, na.rm = TRUE)
}

cox_hr <- function(time, status, score) {
  d <- data.frame(time = time, status = status, score = score)
  d <- d[!is.na(d$time) & !is.na(d$score), ]
  if (nrow(d) < 10 || sum(d$status) < 3) return(data.frame(N = nrow(d), Events = sum(d$status), HR = NA, LCL = NA, UCL = NA, P = NA))
  fit <- tryCatch(coxph(Surv(time, status) ~ scale(score), data = d), error = function(e) NULL)
  if (is.null(fit)) return(data.frame(N = nrow(d), Events = sum(d$status), HR = NA, LCL = NA, UCL = NA, P = NA))
  s <- summary(fit)
  data.frame(N = nrow(d), Events = sum(d$status),
             HR = exp(coef(fit)[1]), LCL = exp(confint(fit)[1, 1]), UCL = exp(confint(fit)[1, 2]),
             P = s$coefficients[5])
}

res_all <- data.frame()

# 统一 rbind 辅助（补 Interaction_P 列）
rbind_res <- function(res_all, r) {
  if (!"Interaction_P" %in% colnames(r)) r$Interaction_P <- NA
  rbind(res_all, r)
}

# ============================================================
# 1. TCGA（HPV status from cBioPortal hnsc_tcga_pub）
# ============================================================
cat("\n===== TCGA =====\n")
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), patient = pid,
                 OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]

hpv <- read.csv(file.path(OUT_DIR, "tcga_hpv_status.csv"), stringsAsFactors = FALSE)
hpv$patient <- substr(hpv$sample_id, 1, 12)
hpv_map <- setNames(hpv$hpv_status, hpv$patient)
df$hpv <- hpv_map[df$patient]
cat("TCGA 501 例 HPV 状态覆盖:", sum(!is.na(df$hpv)), "/", nrow(df), "\n")
print(table(df$hpv, useNA = "ifany"))

# 表达: log2(TPM+1) 中 67 基因
expr <- log2(tumor_mat[intersect(genes67, rownames(tumor_mat)), df$sample] + 1)
cat("TCGA 67-gene 覆盖:", nrow(expr), "/", length(genes67), "\n")
score_all <- prog_score(expr)

for (grp in c("Positive", "Negative")) {
  idx <- which(df$hpv == grp)
  r <- cox_hr(df$OS_time[idx], df$OS_status[idx], score_all[idx])
  r <- cbind(Cohort = "TCGA", HPV = grp, r)
  res_all <- rbind_res(res_all, r)
  cat(sprintf("  HPV%s: HR=%.3f (%.3f-%.3f) P=%.4f (N=%d, E=%d)\n",
              grp, r$HR, r$LCL, r$UCL, r$P, r$N, r$Events))
}

# HPV x program interaction（TCGA 501 内, 可及 HPV 样本）
d2 <- df[!is.na(df$hpv), ]
d2$score <- score_all[d2$sample]
d2$hpv_bin <- ifelse(d2$hpv == "Positive", 1, 0)
fit_int <- tryCatch(coxph(Surv(OS_time, OS_status) ~ scale(score) * hpv_bin, data = d2), error = function(e) NULL)
if (!is.null(fit_int)) {
  p_int <- summary(fit_int)$coefficients["scale(score):hpv_bin", "Pr(>|z|)"]
  cat(sprintf("  TCGA HPV x program interaction P = %.4f\n", p_int))
  res_all$Interaction_P[res_all$Cohort == "TCGA" & res_all$HPV == "Positive"] <- p_int
}

# ============================================================
# 2. GSE65858（hpv_dna:ch1）
# ============================================================
cat("\n===== GSE65858 =====\n")
pd <- readRDS(file.path(DATA_DIR, "GSE65858_pdata.Rds"))
surv658 <- readRDS(file.path(DATA_DIR, "GSE65858_surv.Rds"))
expr658 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))

hpv658 <- as.character(pd[["hpv_dna:ch1"]])
hpv658_bin <- ifelse(hpv658 %in% c("HPV16", "Other HPV"), "Positive",
              ifelse(hpv658 == "Negative", "Negative", NA))
cat("GSE65858 HPV 分组:\n"); print(table(hpv658_bin, useNA = "ifany"))

# 表达
g <- intersect(genes67, rownames(expr658))
cat("GSE65858 67-gene 覆盖:", length(g), "/", length(genes67), "\n")
expr658g <- expr658[g, ]
score658 <- prog_score(expr658g)
common <- intersect(names(score658), surv658$sample)
s658 <- surv658[match(common, surv658$sample), ]
s658$score <- score658[common]
s658$hpv <- hpv658_bin[match(common, rownames(pd))]

for (grp in c("Positive", "Negative")) {
  idx <- which(s658$hpv == grp)
  r <- cox_hr(s658$OS_time_days[idx], s658$OS_status[idx], s658$score[idx])
  r <- cbind(Cohort = "GSE65858", HPV = grp, r)
  res_all <- rbind_res(res_all, r)
  cat(sprintf("  HPV%s: HR=%.3f (%.3f-%.3f) P=%.4f (N=%d, E=%d)\n",
              grp, r$HR, r$LCL, r$UCL, r$P, r$N, r$Events))
}
# interaction
d658 <- s658[!is.na(s658$hpv), ]
d658$hpv_bin <- ifelse(d658$hpv == "Positive", 1, 0)
fit_int658 <- tryCatch(coxph(Surv(OS_time_days, OS_status) ~ scale(score) * hpv_bin, data = d658), error = function(e) NULL)
if (!is.null(fit_int658)) {
  p_int <- summary(fit_int658)$coefficients["scale(score):hpv_bin", "Pr(>|z|)"]
  cat(sprintf("  GSE65858 HPV x program interaction P = %.4f\n", p_int))
  res_all$Interaction_P[res_all$Cohort == "GSE65858" & res_all$HPV == "Positive"] <- p_int
}

# ============================================================
# 3. GSE41613（100% HPV-negative OSCC）与 GSE42743（无 HPV 数据）
# ============================================================
cat("\n===== GSE41613（100% HPV-） =====\n")
expr416 <- readRDS(file.path(DATA_DIR, "GSE41613_expr_full.Rds"))
surv416 <- readRDS(file.path(DATA_DIR, "GSE41613_surv.Rds"))
g416 <- intersect(genes67, rownames(expr416))
cat("GSE41613 67-gene 覆盖:", length(g416), "/", length(genes67), "\n")
score416 <- prog_score(expr416[g416, ])
common416 <- intersect(names(score416), surv416$sample)
s416 <- surv416[match(common416, surv416$sample), ]
s416$score <- score416[common416]
# GSE41613 OS_time 单位: 月（P0A 用 OS_time_months 或 OS_time? 检查列名）
time_col <- ifelse("OS_time_months" %in% colnames(s416), "OS_time_months",
            ifelse("OS_time" %in% colnames(s416), "OS_time", "os_time"))
cat("GSE41613 时间列:", time_col, "\n")
r416 <- cox_hr(s416[[time_col]], s416$OS_status, s416$score)
r416 <- cbind(Cohort = "GSE41613", HPV = "Negative", r416)
res_all <- rbind_res(res_all, r416)
cat(sprintf("  HPV-(100%%): HR=%.3f (%.3f-%.3f) P=%.4f (N=%d, E=%d)\n",
            r416$HR, r416$LCL, r416$UCL, r416$P, r416$N, r416$Events))

cat("\n===== GSE42743（无 HPV 数据，口腔为主） =====\n")
expr427 <- readRDS(file.path(DATA_DIR, "GSE42743_processed", "GSE42743_expression_log2RMA_gene.Rds"))
clin427 <- read.csv(file.path(DATA_DIR, "GSE42743_processed", "GSE42743_clinical.csv"), check.names = FALSE)
g427 <- intersect(genes67, rownames(expr427))
cat("GSE42743 67-gene 覆盖:", length(g427), "/", length(genes67), "\n")
score427 <- prog_score(expr427[g427, ])
# clinical 样本匹配
s427 <- data.frame(sample = clin427$sample_id, OS_time = clin427$OS_time, OS_status = clin427$OS_status)
common427 <- intersect(names(score427), s427$sample)
s427 <- s427[match(common427, s427$sample), ]
s427$score <- score427[common427]
r427 <- cox_hr(s427$OS_time, s427$OS_status, s427$score)
r427 <- cbind(Cohort = "GSE42743", HPV = "Not assayed", r427)
res_all <- rbind_res(res_all, r427)
cat(sprintf("  HPV未测: HR=%.3f (%.3f-%.3f) P=%.4f (N=%d, E=%d)\n",
            r427$HR, r427$LCL, r427$UCL, r427$P, r427$N, r427$Events))

# ============================================================
# 4. 输出
# ============================================================
write.csv(res_all, file.path(OUT_DIR, "hpv_stratified_program_cox_67g.csv"), row.names = FALSE)
cat("\n===== 跨队列 HPV 分组汇总 =====\n")
print(res_all, digits = 4, row.names = FALSE)
cat("\nDONE\n")
