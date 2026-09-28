## P0: 主图 501 全队列重做（V2 RED 版）
## 产出: Fig1A/B, Fig2, Fig3A/B/C, Fig4A/B/C, Fig5B/C（Fig5A 由 P0_cibersort_502.R 后台产出）
## 与 V2 论文数值严格一致（引用 locked_model_TCGA501.csv / TCGA501_program_cox.csv / external_locked_TCGA501.csv）

suppressMessages({
  library(survival)
  library(ggplot2)
  library(pheatmap)
  library(timeROC)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/20_figures_main"
V2_DIR   <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

NPG_RED   <- "#E64B35"
NPG_BLUE  <- "#4DBBD5"
NPG_DARK  <- "#3C5488"
NPG_GREY  <- "#8491B4"

# ============================================================
# 0. 加载 TCGA 501 全队列
# ============================================================
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))  # 60660 x 502
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))  # 528

clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
patient_ids <- substr(colnames(tumor_mat), 1, 12)
cm <- match(patient_ids, clin$bcr_patient_barcode)
surv_full <- data.frame(sample = colnames(tumor_mat), patient = patient_ids,
                        OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm],
                        age = clin$age[cm], gender = clin$gender[cm], stage = clin$ajcc_stage[cm],
                        stringsAsFactors = FALSE)
surv501 <- surv_full[!is.na(surv_full$OS_time), ]
cat("TCGA 501 cohort:", nrow(surv501), "samples,", sum(surv501$OS_status), "events\n")

# 73 基因程序（canonical, 与 B1/TCGA501_program_cox.csv 一致）
core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
modules <- list(Core = core_genes, mTOR = mtor_genes, MitoDynamics = mito_genes,
                MetaboImmune = metab_genes, All = unique(c(core_genes, mtor_genes, mito_genes, metab_genes)))

# 均值 Z 程序分数（501）
expr501 <- log2(tumor_mat[intersect(unique(unlist(modules)), rownames(tumor_mat)), surv501$sample] + 1)
z501 <- t(scale(t(expr501)))
prog_scores <- sapply(modules, function(gs) {
  gi <- intersect(gs, rownames(z501))
  colMeans(z501[gi, , drop = FALSE], na.rm = TRUE)
})
# 验证与已报告数值一致
for (m in c("Core","mTOR","MitoDynamics","MetaboImmune","All")) {
  fit <- coxph(Surv(surv501$OS_time, surv501$OS_status) ~ scale(prog_scores[, m]))
  cat(sprintf("  %s: HR=%.3f P=%.4f\n", m, exp(coef(fit)), summary(fit)$coefficients[5]))
}

# 11-gene locked risk score（501）
locked <- read.csv(file.path(V2_DIR, "locked_model_TCGA501.csv"))
coefs <- setNames(locked$Coefficient, locked$Gene)
sel <- locked$Gene
gi <- intersect(sel, rownames(z501))
risk501 <- colSums(coefs[gi] * z501[gi, , drop = FALSE])
fit_risk <- coxph(Surv(surv501$OS_time, surv501$OS_status) ~ scale(risk501))
cat(sprintf("11-gene risk C-index=%.3f HR=%.3f P=%.4g\n",
            summary(fit_risk)$concordance[1], exp(coef(fit_risk)), summary(fit_risk)$coefficients[5]))

# ============================================================
# 1. Fig1A: Overall program KM（中位切分）
# ============================================================
cat(">>> Fig1A: Overall program KM (501)\n")
ov <- prog_scores[, "All"]
grp <- ifelse(ov >= median(ov), "High", "Low")
df1 <- data.frame(OS_time = surv501$OS_time, OS_status = surv501$OS_status, group = grp)
fit1 <- survfit(Surv(OS_time, OS_status) ~ group, data = df1)
lr <- survdiff(Surv(OS_time, OS_status) ~ group, data = df1)
p1 <- 1 - pchisq(lr$chisq, df = 1)
png(file.path(OUT_DIR, "Fig1A_program_KM_501.png"), width = 5.5, height = 5, res = 300, units = "in")
# survfit 颜色按 group 字母序 (High < Low) 分配; 蓝=High, 红=Low
plot(fit1, col = c(NPG_BLUE, NPG_RED), lwd = 2, xlab = "Time (days)", ylab = "Overall survival",
     main = paste0("Mitoxyperiosis program score (median split)\nlog-rank P = ", signif(p1, 3)),
     cex.main = 0.95)
legend("topright", legend = c("High program", "Low program"),
       col = c(NPG_BLUE, NPG_RED), lwd = 2, bty = "n")
mtext(paste0("n = ", nrow(df1)), side = 3, line = 0.2, adj = 0.02, cex = 0.8)
dev.off()
cat("  log-rank P =", signif(p1, 3), "\n")

# ============================================================
# 2. Fig1B: Module score heatmap（501, 均值 Z）
# ============================================================
cat(">>> Fig1B: Module heatmap (501)\n")
hm <- t(prog_scores[, c("Core","mTOR","MitoDynamics","MetaboImmune","All")])
png(file.path(OUT_DIR, "Fig1B_module_heatmap_501.png"), width = 8, height = 3.2, res = 300, units = "in")
pheatmap(hm, cluster_cols = TRUE, cluster_rows = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c(NPG_BLUE, "white", NPG_RED))(100),
         main = "Mitoxyperiosis module scores (501 tumors)",
         fontsize_row = 10)
dev.off()

# ============================================================
# 3. Fig2: Multivariable Cox forest（11-gene risk, 调整 age/gender/stage）
# ============================================================
cat(">>> Fig2: Multivariable Cox forest\n")
df2 <- surv501
df2$risk <- risk501[surv501$sample]
df2$gender2 <- factor(df2$gender)
df2$stage2 <- ifelse(grepl("III|IV", df2$stage), "III/IV", ifelse(grepl("I|II", df2$stage), "I/II", NA))
fit2 <- coxph(Surv(OS_time, OS_status) ~ scale(risk) + age + gender2 + stage2, data = df2)
s2 <- summary(fit2)
coefs2 <- coef(fit2); ci2 <- confint(fit2)
terms2 <- c("Risk score (per SD)", "Age (per year)", "Gender (Female vs Male)", "Stage (III/IV vs I/II)")
hr2 <- exp(coefs2); lcl2 <- exp(ci2[, 1]); ucl2 <- exp(ci2[, 2])
p2 <- s2$coefficients[, 5]
lab2 <- c("Risk score (per SD)", "Age (per year)", "Gender", "Stage (III/IV vs I/II)")
df_forest <- data.frame(term = factor(lab2, levels = rev(lab2)), HR = hr2, LCL = lcl2, UCL = ucl2,
                        P = sprintf("%.4g", p2), stringsAsFactors = FALSE)
png(file.path(OUT_DIR, "Fig2_multivariable_forest_501.png"), width = 7, height = 4, res = 300, units = "in")
par(mar = c(4, 9, 3, 3))
x <- rev(seq_along(lab2))
plot(NA, xlim = c(0, max(ucl2) * 1.15), ylim = c(0.5, length(lab2) + 0.6), yaxt = "n",
     xlab = "Hazard ratio (95% CI)", ylab = "", main = "Multivariable Cox (11-gene risk score, 501 tumors)")
abline(v = 1, lty = 2, col = "grey50")
segments(lcl2, x, ucl2, x, lwd = 2, col = NPG_DARK)
points(hr2, x, pch = 18, cex = 1.8, col = ifelse(p2 < 0.05, NPG_RED, NPG_DARK))
axis(2, at = x, labels = rev(lab2), las = 1, cex.axis = 0.85)
text(rep(max(ucl2) * 1.15, length(x)), x, labels = df_forest$P, pos = 2, cex = 0.8, col = "grey30")
mtext("P", side = 3, line = 0.2, adj = 0.98, cex = 0.8)
dev.off()
cat("  risk per-SD HR=", round(hr2[1], 3), "P=", signif(p2[1], 3), "\n")
write.csv(df_forest, file.path(OUT_DIR, "Fig2_multivariable_cox_501.csv"), row.names = FALSE)

# ============================================================
# 4. Fig3A: timeROC（11-gene risk, 501）
# ============================================================
cat(">>> Fig3A: timeROC (501)\n")
roc3 <- timeROC(T = surv501$OS_time, delta = surv501$OS_status, marker = risk501,
                cause = 1, times = c(365, 1095, 1825), ROC = TRUE, iid = FALSE)
cat("  AUC 1/3/5yr:", round(roc3$AUC, 3), "\n")
png(file.path(OUT_DIR, "Fig3A_timeROC_501.png"), width = 6, height = 5, res = 300, units = "in")
plot(roc3, time = 365, col = NPG_RED, lwd = 2, xlab = "1 - Specificity", ylab = "Sensitivity",
     main = "Time-dependent ROC (11-gene risk score, 501 tumors)")
plot(roc3, time = 1095, col = NPG_BLUE, lwd = 2, add = TRUE)
plot(roc3, time = 1825, col = NPG_DARK, lwd = 2, add = TRUE)
abline(0, 1, lty = 2, col = "grey60")
legend("bottomright",
       legend = c(sprintf("1-yr AUC=%.3f", roc3$AUC[1]),
                  sprintf("3-yr AUC=%.3f", roc3$AUC[2]),
                  sprintf("5-yr AUC=%.3f", roc3$AUC[3])),
       col = c(NPG_RED, NPG_BLUE, NPG_DARK), lwd = 2, bty = "n")
dev.off()

# ============================================================
# 5. Fig3B: KM by median risk（11-gene）
# ============================================================
cat(">>> Fig3B: KM by median risk (501)\n")
grp3 <- ifelse(risk501 >= median(risk501), "High", "Low")
df3 <- data.frame(OS_time = surv501$OS_time, OS_status = surv501$OS_status, group = grp3)
fit3 <- survfit(Surv(OS_time, OS_status) ~ group, data = df3)
lr3 <- survdiff(Surv(OS_time, OS_status) ~ group, data = df3)
p3 <- 1 - pchisq(lr3$chisq, df = 1)
png(file.path(OUT_DIR, "Fig3B_risk_KM_501.png"), width = 5.5, height = 5, res = 300, units = "in")
# 字母序: High < Low → 蓝=High, 红=Low
plot(fit3, col = c(NPG_BLUE, NPG_RED), lwd = 2, xlab = "Time (days)", ylab = "Overall survival",
     main = paste0("Risk score (median split, 501 tumors)\nlog-rank P = ", signif(p3, 3)),
     cex.main = 0.95)
legend("topright", legend = c("High risk", "Low risk"),
       col = c(NPG_BLUE, NPG_RED), lwd = 2, bty = "n")
dev.off()

# ============================================================
# 6. Fig3C: Risk score distribution by OS status
# ============================================================
cat(">>> Fig3C: Risk distribution by OS status\n")
df3c <- data.frame(risk = risk501, status = ifelse(surv501$OS_status == 1, "Dead", "Alive"))
p3c <- ggplot(df3c, aes(x = status, y = risk, fill = status)) +
  geom_boxplot(outlier.size = 0.6, width = 0.5) +
  scale_fill_manual(values = c("Alive" = NPG_BLUE, "Dead" = NPG_RED), guide = "none") +
  labs(x = "", y = "Risk score (11-gene)", title = "Risk score by survival status (501 tumors)") +
  theme_bw() + theme(plot.title = element_text(size = 11))
ggsave(file.path(OUT_DIR, "Fig3C_risk_distribution_501.png"), p3c, width = 4.5, height = 4.5, dpi = 300)
w <- wilcox.test(risk ~ status, data = df3c)
cat("  Wilcoxon P =", signif(w$p.value, 3), "\n")

# ============================================================
# 7. Fig4A/B/C: 外部验证 forest + meta + C-index（501 locked）
# ============================================================
cat(">>> Fig4: External evaluation (501 locked model)\n")

# --- 用 V2_full_cohort_analysis.R 一致的逻辑加载外部队列并补算 C-index ---
ext_cindex <- data.frame()
ext_forest <- read.csv(file.path(V2_DIR, "external_locked_TCGA501.csv"))
meta_os <- read.csv(file.path(V2_DIR, "external_meta_TCGA501.csv"))

# GSE41613
expr41613 <- readRDS(file.path(DATA_DIR, "GSE41613_expr.Rds"))
surv41613 <- readRDS(file.path(DATA_DIR, "GSE41613_surv.Rds"))
rn41613 <- rownames(expr41613)
gene_map41613 <- c(PRKN = "PARK2", VDAC1 = "VDAC1", SLC25A5 = "SLC25A5", HK1 = "HK1", MAOB = "MAOB")
avail41613 <- sel[sel %in% names(gene_map41613) & gene_map41613[sel] %in% rn41613]
sub41613 <- expr41613[gene_map41613[avail41613], , drop = FALSE]
z41613 <- t(scale(t(sub41613)))
risk41613 <- colSums(coefs[avail41613] * z41613)
s416 <- surv41613[match(colnames(sub41613), surv41613$sample), ]
s416$risk <- risk41613[colnames(sub41613)]
s416 <- s416[!is.na(s416$risk), ]
# 带方向 C-index = 1 - concordance()（concordance() 约定: x 小→生存差; 取反后: x 大→生存差, TCGA 方向）
c416 <- 1 - concordance(Surv(s416$OS_time_months, s416$OS_status) ~ s416$risk)$concordance
ext_cindex <- rbind(ext_cindex, data.frame(Cohort = "GSE41613", C_index = c416, N = nrow(s416)))

# GSE42743
expr42743 <- readRDS(file.path(DATA_DIR, "GSE42743_processed/GSE42743_expression_log2RMA_gene.Rds"))
clin42743 <- read.csv(file.path(DATA_DIR, "GSE42743_processed/GSE42743_clinical.csv"), check.names = FALSE)
rn42743 <- rownames(expr42743)
avail42743 <- sel[sel %in% rn42743]
sub42743 <- expr42743[avail42743, , drop = FALSE]
z42743 <- t(scale(t(sub42743)))
risk42743 <- colSums(coefs[avail42743] * z42743)
c427 <- clin42743
c427$risk <- risk42743[as.character(clin42743$sample_id)]
c427 <- c427[!is.na(c427$risk), ]
c427_c <- 1 - concordance(Surv(c427$os_time, c427$os_status) ~ c427$risk)$concordance
ext_cindex <- rbind(ext_cindex, data.frame(Cohort = "GSE42743", C_index = c427_c, N = nrow(c427)))

# GSE27020 (DFS)
lines27020 <- readLines(gzfile(file.path(DATA_DIR, "GSE27020_series_matrix.txt.gz"), open = "rt"), warn = FALSE)
start27020 <- grep("!series_matrix_table_begin", lines27020)
sample27020 <- gsub('"', "", strsplit(lines27020[start27020 + 1], "\t")[[1]][-1])
probes27020 <- c(PRKN = "207058_s_at", VDAC1 = "212038_s_at", SLC25A5 = "200657_at",
                 HK1 = "200697_at", MAOB = "204041_at")
extract_row <- function(probe) {
  hit <- grep(paste0('"', probe, '"'), lines27020[(start27020 + 2):length(lines27020)], fixed = TRUE)
  if (length(hit) == 0) return(NULL)
  as.numeric(strsplit(lines27020[start27020 + 1 + hit[1]], "\t")[[1]][-1])
}
avail27020 <- sel[sel %in% names(probes27020)]
expr27020_sub <- matrix(NA, nrow = length(avail27020), ncol = length(sample27020),
                        dimnames = list(avail27020, sample27020))
for (g in avail27020) {
  v <- extract_row(probes27020[[g]])
  if (!is.null(v)) expr27020_sub[g, ] <- v
}
keep27020 <- apply(expr27020_sub, 2, function(x) !any(is.na(x)))
if (sum(keep27020) > 30 && nrow(expr27020_sub) >= 4) {
  z27020 <- t(scale(t(expr27020_sub[, keep27020])))
  risk27020 <- colSums(coefs[avail27020] * z27020)
  dfs27020 <- read.csv(file.path(DATA_DIR, "output_v5/GSE27020_clinical_risk_data.csv"), check.names = FALSE)
  dfs27020$risk <- risk27020[as.character(dfs27020$sample)]
  dfs27020 <- dfs27020[!is.na(dfs27020$risk), ]
  c27020 <- 1 - concordance(Surv(dfs27020$DFS_time, dfs27020$DFS_event) ~ dfs27020$risk)$concordance
  ext_cindex <- rbind(ext_cindex, data.frame(Cohort = "GSE27020 (DFS)", C_index = c27020, N = nrow(dfs27020)))
}

# GSE65858
expr65858 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))
surv65858 <- readRDS(file.path(DATA_DIR, "GSE65858_surv.Rds"))
rn65858 <- rownames(expr65858)
avail65858 <- sel[sel %in% rn65858]
sub65858 <- expr65858[avail65858, , drop = FALSE]
z65858 <- t(scale(t(sub65858)))
risk65858 <- colSums(coefs[avail65858] * z65858)
s658 <- surv65858[match(colnames(sub65858), surv65858$sample), ]
s658$risk <- risk65858[colnames(sub65858)]
s658 <- s658[!is.na(s658$risk), ]
c658 <- 1 - concordance(Surv(s658$OS_time_days, s658$OS_status) ~ s658$risk)$concordance
ext_cindex <- rbind(ext_cindex, data.frame(Cohort = "GSE65858", C_index = c658, N = nrow(s658)))
# TCGA 参考（带方向, HR>1 → coxph C 即带方向 C）
ext_cindex <- rbind(ext_cindex, data.frame(Cohort = "TCGA (discovery)", C_index = summary(fit_risk)$concordance[1], N = nrow(surv501)))
write.csv(ext_cindex, file.path(OUT_DIR, "Fig4C_external_cindex_501.csv"), row.names = FALSE)
cat("  External C-indices:\n"); print(ext_cindex)

# --- Fig4A: forest（TCGA ref + 4 external）---
tcga_hr <- exp(coef(fit_risk)); tcga_ci <- exp(confint(fit_risk))
fa <- data.frame(
  Cohort = c("TCGA (discovery)", as.character(ext_forest$Cohort)),
  Model  = c("11-gene locked", as.character(ext_forest$Model)),
  HR = c(tcga_hr, ext_forest$HR_per_SD),
  LCL = c(tcga_ci[1], ext_forest$LCL),
  UCL = c(tcga_ci[2], ext_forest$UCL),
  P = c(summary(fit_risk)$coefficients[5], ext_forest$P),
  stringsAsFactors = FALSE
)
fa$is_tcga <- fa$Cohort == "TCGA (discovery)"
fa$Cohort2 <- ifelse(fa$Cohort == "TCGA (discovery)", "TCGA (discovery, 501 tumors)",
                     paste0(fa$Cohort, " (", fa$Model, ")"))
fa$Cohort2 <- factor(fa$Cohort2, levels = rev(fa$Cohort2))
p4a <- ggplot(fa, aes(x = HR, y = Cohort2)) +
  geom_point(size = 3.5, aes(color = is_tcga)) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2, linewidth = 0.9) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("TRUE" = NPG_RED, "FALSE" = NPG_DARK), guide = "none") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 3)) +
  labs(title = "Locked-model external evaluation (HR per 1 SD, 501-cohort model)",
       x = "Hazard ratio (log scale)", y = "") +
  theme_bw() + theme(axis.text.y = element_text(size = 9))
ggsave(file.path(OUT_DIR, "Fig4A_external_forest_501.png"), p4a, width = 7.5, height = 4.5, dpi = 300)

# --- Fig4B: meta forest（3 OS cohorts）---
os_rows <- ext_forest[ext_forest$Endpoint == "OS", ]
mb <- data.frame(
  Cohort = c(as.character(os_rows$Cohort), "Pooled (REML)"),
  HR = c(os_rows$HR_per_SD, meta_os$pooled_HR),
  LCL = c(os_rows$LCL, meta_os$LCL),
  UCL = c(os_rows$UCL, meta_os$UCL),
  stringsAsFactors = FALSE
)
mb$Cohort <- factor(mb$Cohort, levels = rev(mb$Cohort))
p4b <- ggplot(mb, aes(x = HR, y = Cohort)) +
  geom_point(size = 3.5, color = c(rep(NPG_DARK, nrow(os_rows)), NPG_RED)) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.2, linewidth = 0.9) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  scale_x_log10() +
  labs(title = sprintf("Random-effects meta-analysis (3 OS cohorts)\npooled HR=%.2f (95%% CI %.2f-%.2f), I2=%.1f%%",
                       meta_os$pooled_HR, meta_os$LCL, meta_os$UCL, meta_os$I2),
       x = "Hazard ratio (log scale)", y = "") +
  theme_bw() + theme(axis.text.y = element_text(size = 9))
ggsave(file.path(OUT_DIR, "Fig4B_meta_501.png"), p4b, width = 7, height = 4, dpi = 300)

# --- Fig4C: cross-cohort C-index ---
ext_cindex$Cohort <- factor(ext_cindex$Cohort, levels = rev(ext_cindex$Cohort))
p4c <- ggplot(ext_cindex, aes(x = C_index, y = Cohort)) +
  geom_point(size = 3.5, color = NPG_DARK) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  scale_x_continuous(limits = c(0.3, 0.8), breaks = seq(0.3, 0.8, 0.1)) +
  labs(title = "Cross-cohort C-index (locked 11-gene model)", x = "C-index", y = "") +
  theme_bw() + theme(axis.text.y = element_text(size = 9))
ggsave(file.path(OUT_DIR, "Fig4C_cindex_501.png"), p4c, width = 6, height = 4, dpi = 300)
cat("  Fig4A/B/C saved\n")

# ============================================================
# 8. Fig5B: MetaboImmune vs mTOR（501 重算）
# ============================================================
cat(">>> Fig5B: MetaboImmune vs mTOR (501)\n")
r_mm <- cor(prog_scores[, "MetaboImmune"], prog_scores[, "mTOR"], method = "spearman")
n_mm <- length(prog_scores[, "MetaboImmune"])
t_mm <- r_mm * sqrt((n_mm - 2) / (1 - r_mm^2))
p_mm <- 2 * pt(-abs(t_mm), df = n_mm - 2)
df5b <- data.frame(x = prog_scores[, "MetaboImmune"], y = prog_scores[, "mTOR"])
p5b <- ggplot(df5b, aes(x = x, y = y)) +
  geom_point(color = NPG_DARK, alpha = 0.5, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, color = NPG_RED, linewidth = 0.8) +
  labs(x = "MetaboImmune module score", y = "mTOR module score",
       title = sprintf("Spearman r = %.3f, P = %.2e", r_mm, p_mm)) +
  theme_bw() + theme(plot.title = element_text(size = 10))
ggsave(file.path(OUT_DIR, "Fig5B_metabo_mTOR_501.png"), p5b, width = 5.5, height = 5, dpi = 300)
cat("  r =", round(r_mm, 3), "P =", signif(p_mm, 3), "\n")
write.csv(data.frame(Variable1 = "MetaboImmune", Variable2 = "mTOR",
                     Spearman_r = r_mm, P = p_mm, N = n_mm),
          file.path(OUT_DIR, "Fig5B_metabo_mTOR_501.csv"), row.names = FALSE)

# ============================================================
# 9. Fig5C: CD274 vs overall program（501 重算；CD274 单独从 tumor_mat 取 log2 TPM）
# ============================================================
cat(">>> Fig5C: CD274 vs overall program (501)\n")
if ("CD274" %in% rownames(tumor_mat)) {
  cd274 <- log2(tumor_mat["CD274", surv501$sample] + 1)  # 原始 TPM → log2
  cd274 <- as.numeric(cd274)
  r_pdl <- cor(cd274, prog_scores[, "All"], method = "spearman")
  t_pdl <- r_pdl * sqrt((n_mm - 2) / (1 - r_pdl^2))
  p_pdl <- 2 * pt(-abs(t_pdl), df = n_mm - 2)
  df5c <- data.frame(x = prog_scores[, "All"], y = cd274)
  p5c <- ggplot(df5c, aes(x = x, y = y)) +
    geom_point(color = NPG_DARK, alpha = 0.5, size = 1.6) +
    geom_smooth(method = "lm", se = TRUE, color = NPG_RED, linewidth = 0.8) +
    labs(x = "Overall mitoxyperiosis program score", y = "CD274 (PD-L1) expression (log2 TPM)",
         title = sprintf("Spearman rho = %.3f, P = %.2e", r_pdl, p_pdl)) +
    theme_bw() + theme(plot.title = element_text(size = 10))
  ggsave(file.path(OUT_DIR, "Fig5C_CD274_program_501.png"), p5c, width = 5.5, height = 5, dpi = 300)
  cat("  rho =", round(r_pdl, 3), "P =", signif(p_pdl, 3), "\n")
  write.csv(data.frame(Var1 = "CD274", Var2 = "Overall program",
                       Spearman_rho = r_pdl, P = p_pdl, N = n_mm),
            file.path(OUT_DIR, "Fig5C_CD274_program_501.csv"), row.names = FALSE)
} else {
  cat("  CD274 不在表达矩阵（跳过 Fig5C）\n")
}

cat("\n=== P0 FIGURES 501 DONE ===\n")
