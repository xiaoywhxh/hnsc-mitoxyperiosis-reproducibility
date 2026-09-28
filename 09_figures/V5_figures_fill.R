## ============================================================
## V5 缺失主图补齐：Fig1A / Fig3A / Fig3B / Fig3D
## 数据源与 V5 论文 Results 3.2/3.6 一致（67-gene 主口径）
## ============================================================
suppressMessages({library(survival); library(ggplot2)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/20_figures_main"

NPG_RED  <- "#E64B35"; NPG_BLUE <- "#3C5488"; NPG_TEAL <- "#4DBBD5"; NPG_GREY <- "#8491B4"

tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]

common67 <- read.csv(file.path(DATA_DIR, "..", "V2_Reanalysis", "08_external_validation", "program_common_genes_3OS.csv"))$common_gene

core <- c('PRKN','VDAC1','VDAC2','VDAC3','BAX','BAK1','BID','BBC3','PMAIP1','BCL2','BCL2L1','MCL1','TSPO','AIFM1')
mtor <- c('MTOR','RICTOR','RPTOR','MLST8','MAPKAP1','PRR5','PRR5L','DEPTOR','TSC1','TSC2','RHEB','AKT1','AKT2','AKT3','PTEN','PIK3CA','PIK3CB','PIK3CD','RRAGA','RRAGB','RRAGC','RRAGD')
mito <- c('DNM1L','FIS1','MFF','MIEF1','MIEF2','MFN1','MFN2','OPA1','MARCHF5','PINK1','SLC25A3','SLC25A5','SLC25A6','CYCS')
metab <- c('HK1','HK2','PFKFB3','PKM','LDHA','LDHB','IDH1','IDH2','MDH1','MDH2','CS','GLS','GLS2','GOT1','GOT2','NLRP3','CASP1','CASP4','CASP5','STING1','CGAS','MAOB','ENDOG')
all73 <- unique(c(core, mtor, mito, metab))

calc_prog <- function(genes) {
  genes <- intersect(genes, rownames(tumor_mat))
  ex <- log2(tumor_mat[genes, df$sample] + 1)
  z <- t(scale(t(ex)))
  colMeans(z, na.rm = TRUE)
}

# ============================================================
# Fig1A: 模块级 + overall 森林图（67-gene 主口径，BH q-values）
# ============================================================
cat(">>> Fig1A: module forest (67-gene)\n")
mods <- list(Core = intersect(core, common67), mTOR = intersect(mtor, common67),
             MitoDynamics = intersect(mito, common67), MetaboImmune = intersect(metab, common67),
             All = common67)
mod_res <- lapply(names(mods), function(m) {
  sc <- calc_prog(mods[[m]])
  fit <- coxph(Surv(OS_time, OS_status) ~ scale(sc), data = df)
  ci <- exp(confint(fit))
  data.frame(Module = m, Genes = length(mods[[m]]), HR = exp(coef(fit)), LCL = ci[1], UCL = ci[2], P = summary(fit)$coefficients[5])
})
mod_df <- do.call(rbind, mod_res)
mod_df$q <- p.adjust(mod_df$P, method = "BH")
# 顺序: MetaboImmune, MitoDynamics, Core, All, mTOR（按 HR 降序，与正文一致）
mod_df$Module <- factor(mod_df$Module, levels = c("mTOR", "All", "Core", "MitoDynamics", "MetaboImmune"))
mod_df$sig <- mod_df$q < 0.05

p1a <- ggplot(mod_df, aes(x = HR, y = Module)) +
  geom_point(size = 3.6, aes(color = sig)) +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.18, linewidth = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  scale_color_manual(values = c("TRUE" = NPG_RED, "FALSE" = NPG_GREY), guide = "none") +
  scale_x_continuous(breaks = seq(0.8, 1.5, 0.1)) +
  labs(x = "Per-1-SD HR (95% CI), 67-gene common program", y = "",
       title = "Mitoxyperiosis program scores and OS in TCGA-HNSC (n=501)\nBH q-values: MetaboImmune 0.007, MitoDynamics 0.009, Core 0.040, All 0.026; mTOR NS") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(size = 10.5))
ggsave(file.path(OUT_DIR, "Fig1A_module_forest_67g.png"), p1a, width = 6.2, height = 3.4, dpi = 300)
write.csv(mod_df, file.path(OUT_DIR, "Fig1A_module_forest_67g.csv"), row.names = FALSE)
cat("  Fig1A saved; q-values:", paste(round(mod_df$q, 4), collapse = ", "), "\n")

# ============================================================
# Fig3A: CIBERSORT continuous 关联（22 细胞 ρ 点图）
# ============================================================
cat(">>> Fig3A: CIBERSORT continuous\n")
cont <- read.csv(file.path(OUT_DIR, "Fig5A_cibersort_continuous_501.csv"))
cont$sig <- cont$FDR < 0.05
cont$CellType <- reorder(cont$CellType, cont$rho)

p3a <- ggplot(cont, aes(x = rho, y = CellType)) +
  geom_point(size = 3.2, aes(color = sig)) +
  geom_vline(xintercept = 0, color = "grey45", linetype = "dashed") +
  scale_color_manual(values = c("TRUE" = NPG_RED, "FALSE" = NPG_GREY), name = "BH-FDR < 0.05") +
  labs(x = "Spearman rho vs overall program score", y = "",
       title = "CIBERSORT (LM22, relative mode): continuous associations with program score") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(size = 10.5))
ggsave(file.path(OUT_DIR, "Fig3A_cibersort_continuous.png"), p3a, width = 6.8, height = 5.2, dpi = 300)
cat("  Fig3A saved; FDR<0.05:", sum(cont$sig), "/22\n")

# ============================================================
# Fig3B: CD4 memory resting T 分数 by program subgroup
# ============================================================
cat(">>> Fig3B: CD4 memory resting T by program subgroup\n")
ciber <- readRDS(file.path(OUT_DIR, "Step8_cibersort_502.Rds"))
prog67 <- calc_prog(common67)
common_s <- intersect(names(prog67), rownames(ciber))
sub_df <- data.frame(sample = common_s, prog = prog67[common_s],
                     cd4mr = ciber[common_s, "T cells CD4 memory resting"])
sub_df$group <- ifelse(sub_df$prog > median(sub_df$prog), "High program", "Low program")
sub_df$group <- factor(sub_df$group, levels = c("Low program", "High program"))
p3b <- ggplot(sub_df, aes(x = group, y = cd4mr, fill = group)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85, width = 0.55) +
  scale_fill_manual(values = c("Low program" = NPG_BLUE, "High program" = NPG_RED), guide = "none") +
  labs(x = "", y = "CD4 memory resting T-cell fraction (CIBERSORT)",
       title = sprintf("CD4 memory resting T cells by program subgroup (Wilcoxon P=%.2g)", wilcox.test(cd4mr ~ group, data = sub_df)$p.value)) +
  theme_bw(base_size = 11) + theme(plot.title = element_text(size = 10.5))
ggsave(file.path(OUT_DIR, "Fig3B_CD4_memory_subgroup.png"), p3b, width = 4.6, height = 4.4, dpi = 300)
cat(sprintf("  Fig3B saved; High=%.3f vs Low=%.3f\n",
            mean(sub_df$cd4mr[sub_df$group == "High program"]), mean(sub_df$cd4mr[sub_df$group == "Low program"])))

# ============================================================
# Fig3D: 模块相关算法敏感性（mean-Z / singscore / ssGSEA）
# ============================================================
cat(">>> Fig3D: module correlation sensitivity\n")
sens <- data.frame(
  Method = c("Program-mean Z", "Standard singscore", "Earlier ssGSEA"),
  r = c(0.536, -0.08, -0.31)
)
sens$Method <- factor(sens$Method, levels = rev(c("Program-mean Z", "Standard singscore", "Earlier ssGSEA")))
p3d <- ggplot(sens, aes(x = r, y = Method)) +
  geom_col(aes(fill = r > 0), width = 0.55, alpha = 0.9) +
  geom_vline(xintercept = 0, color = "grey40") +
  scale_fill_manual(values = c("TRUE" = NPG_RED, "FALSE" = NPG_BLUE), guide = "none") +
  geom_text(aes(label = sprintf("%+.2f", r)), hjust = ifelse(sens$r > 0, -0.2, 1.2), size = 4) +
  scale_x_continuous(limits = c(-1, 0.9)) +
  labs(x = "Spearman r (MetaboImmune vs mTOR)", y = "",
       title = "Module-correlation sensitivity to scoring method") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(size = 10.5))
ggsave(file.path(OUT_DIR, "Fig3D_module_corr_sensitivity.png"), p3d, width = 5.4, height = 3.0, dpi = 300)
cat("  Fig3D saved\n")

cat("\nDONE: Fig1A/Fig3A/Fig3B/Fig3D all generated\n")
