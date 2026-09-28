## v10 Supplementary 图件：S22 B/C 面板 + S24 purity 图
suppressMessages({library(ggplot2)})
OUT <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"

## ---------- S22 B: Cox attenuation ----------
cox <- read.csv(file.path(OUT, "P0-1_specificity_cox.csv"))
cox$Model <- factor(cox$Model, levels = cox$Model)
pB <- ggplot(cox, aes(x = Model, y = HR, ymin = LCL, ymax = UCL)) +
  geom_hline(yintercept = 1, lty = 3, color = "grey40") +
  geom_pointrange(color = "#E64B35", size = 0.7, fatten = 3) +
  scale_y_log10(limits = c(0.5, 3.0)) +
  coord_flip() +
  labs(x = NULL, y = "Per-1-SD HR (log scale)", title = "B  Multivariable attenuation of the program HR") +
  theme_bw(base_size = 11) + theme(axis.text.y = element_text(size = 8.5))

## ---------- S22 C: VIF ----------
vif2 <- read.csv(file.path(OUT, "P0-1_specificity_vif.csv"))
vif4 <- read.csv(file.path(OUT, "P0-1_specificity_vif_model4.csv"))
vif2$Model <- "Model 2"; vif4$Model <- "Model 4"
vif <- rbind(vif2[, c("Variable", "VIF", "Model")], vif4[, c("Variable", "VIF", "Model")])
vif <- vif[vif$Variable == "mitoxy67", ]
pC <- ggplot(vif, aes(x = Model, y = VIF, fill = Model)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 10, lty = 2, color = "#E64B35") +
  annotate("text", x = 0.55, y = 10.5, label = "VIF = 10 (severe)", size = 3.2, color = "#E64B35") +
  geom_text(aes(label = round(VIF, 1)), vjust = -0.4, size = 4) +
  coord_cartesian(ylim = c(0, 22)) +
  scale_fill_manual(values = c("Model 2" = "#4DBBD5", "Model 4" = "#E64B35")) +
  labs(x = NULL, y = "VIF (mitoxy67)", title = "C  Variance inflation factor") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

## ---------- S24: purity-adjusted immune ----------
p0 <- read.csv(file.path(OUT, "P0-3_purity_adjusted_immune.csv"))
## 提取 CD4 与 CD274 的 raw/adjusted
cd4_raw <- p0$Value[p0$Metric == "CD4_raw_spearman"]
cd4_adj <- p0$Value[p0$Metric == "CD4_partial_rho"]
cd4_p   <- p0$Value[p0$Metric == "CD4_partial_P"]
dat <- data.frame(
  metric = rep(c("CD4 memory resting T", "CD274 transcript"), each = 2),
  type = rep(c("Raw", "Purity-adjusted"), 2),
  est = c(cd4_raw, cd4_adj, NA, NA))
## CD274 beta raw/adj 从回归输出（脚本日志: raw beta 0.614 P=2.5e-5; adj beta 0.579 P=2.4e-7 早期；现 CSV 有 CD4_reg_*）
cd274_raw <- 0.614; cd274_adj <- 0.579
dat$est[3] <- cd274_raw; dat$est[4] <- cd274_adj
dat$metric <- factor(dat$metric, levels = c("CD4 memory resting T", "CD274 transcript"))
pS <- ggplot(dat, aes(x = metric, y = est, fill = type)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("Raw" = "#4DBBD5", "Purity-adjusted" = "#E64B35")) +
  labs(x = NULL, y = "Association with program score\n(Spearman rho / standardized beta)",
       title = "Immune associations before vs after ESTIMATE-purity adjustment",
       subtitle = "CD4 partial rho = 0.259 (P = 4.1e-09); CD274 adj beta = 0.58 (P = 2.4e-07); median purity 0.745") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

ggsave(file.path(OUT, "FigS22_BC_panels.png"), pB + pC + plot_layout(ncol = 1, heights = c(1.2, 1)), width = 7, height = 7.5, dpi = 300)
ggsave(file.path(OUT, "FigS24_purity_immune.png"), pS, width = 7, height = 4.5, dpi = 300)
cat("S22 BC + S24 generated\n")
