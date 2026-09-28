## ============================================================
## FigS21 CANONICAL v2 (2026-09-27) — 修复拼图规范
## 修复：(1) Panel A 用 annotation_custom+rasterGrob 渲染（不再空白）
##       (2) 三面板统一对齐：A 全宽在上，B|C 等宽并排在下
##       (3) Panel B 去掉面板内文字标签（数字见表 S19），消除线-字重叠
##       (4) Panel C "VIF=10" 标签移至面板左上空白区，不再被裁切/压柱
## NPG 配色：#E64B35 / #4DBBD5
## ============================================================
suppressMessages({library(ggplot2); library(patchwork); library(png); library(grid)})
OUT_DIR <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"

## ---------- Panel A：相关矩阵热图（pheatmap 产物，raster 全幅嵌入） ----------
imgA <- readPNG(file.path(OUT_DIR, "P0-1_specificity_corheatmap.png"))
pA <- ggplot() +
  annotation_custom(rasterGrob(imgA, interpolate = TRUE),
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

## ---------- Panel B：多变量校正后 program HR 森林图（无面板内文字） ----------
cox <- read.csv(file.path(OUT_DIR, "P0-1_specificity_cox.csv"))
cox$Model <- factor(cox$Model, levels = rev(cox$Model))  # 1 在底部
pB <- ggplot(cox, aes(x = HR, y = Model)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL), height = 0.18,
                 linewidth = 0.9, color = "#E64B35") +
  geom_point(size = 2.5, color = "#E64B35") +
  scale_x_continuous(limits = c(0.55, 3.15), breaks = c(0.5, 1, 1.5, 2, 2.5, 3),
                     trans = "log10") +
  labs(x = "Per-1-SD HR (log scale)", y = NULL,
       title = "B  Multivariable attenuation of the program HR") +
  theme_bw(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", size = 11.5),
        axis.text.y = element_text(size = 8.8),
        plot.margin = margin(4, 6, 4, 4))

## ---------- Panel C：VIF 柱状图（阈值标签置于左上空白区） ----------
vif2 <- read.csv(file.path(OUT_DIR, "P0-1_specificity_vif.csv"))
vif4 <- read.csv(file.path(OUT_DIR, "P0-1_specificity_vif_model4.csv"))
vif_df <- data.frame(
  Model = c("Model 2\n(+OXPHOS+Glyco+\nApop+Prolif)", "Model 4\n(+Mitophagy+\nMitoDyn)"),
  VIF   = c(vif2$VIF[vif2$Variable == "mitoxy67"], vif4$VIF[vif4$Variable == "mitoxy67"]),
  col   = c("#4DBBD5", "#E64B35"))
vif_df$Model <- factor(vif_df$Model, levels = vif_df$Model)
pC <- ggplot(vif_df, aes(x = Model, y = VIF, fill = col)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "#E64B35", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f", VIF)), vjust = -0.5, size = 3.6, fontface = "bold") +
  annotate("text", x = 0.58, y = 21.2, label = "VIF = 10: severe collinearity threshold",
           hjust = 0, size = 2.9, color = "#E64B35", fontface = "italic") +
  scale_fill_identity() +
  scale_y_continuous(limits = c(0, 22.5), breaks = seq(0, 20, 5)) +
  labs(x = NULL, y = "VIF (mitoxy67)",
       title = "C  Variance inflation factor") +
  theme_bw(base_size = 10.5) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 11.5),
        axis.text.x = element_text(size = 8.4),
        plot.margin = margin(4, 6, 4, 4))

## ---------- 拼图：A 全宽在上，B|C 等宽并排在下 ----------
combined <- pA / (pB | pC) + plot_layout(heights = c(2.1, 1))

out_png <- file.path(OUT_DIR, "FigS21_specificity_benchmark_CANONICAL_20260927.png")
ggsave(out_png, combined, width = 10, height = 12.3, dpi = 300)
cat("FigS21 CANONICAL v2 saved:", out_png, "\n")
