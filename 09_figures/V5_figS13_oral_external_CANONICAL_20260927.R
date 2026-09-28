## ============================================================
## FigS13 CANONICAL v2 (2026-09-27) — 修复版
## 修复：(1) 数值列固定 x 位置，不再与误差线重叠
##       (2) 副标题缩短为两行，不再右缘截断
##       (3) 行序显式锁定（TCGA overall → oral → GSE41613 → GSE42743 overall/oral → GSE65858 oral）
## NPG 双色：TCGA = #E64B35，非 TCGA = #4DBBD5
## ============================================================
suppressMessages({library(ggplot2)})
OUT <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/20_figures_main"

d <- read.csv("D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation/FigureS13_oral_external_data_CANONICAL.csv")

## 显式行序（自上而下）+ 校验
order_cohorts <- c("TCGA overall (reference)", "TCGA oral subset",
                   "GSE41613 (100% oral)", "GSE42743 overall (96% oral)",
                   "GSE42743 oral subset", "GSE65858 oral subset")
key <- setNames(seq_along(order_cohorts), order_cohorts)
d$ord <- key[d$Cohort]
d <- d[order(d$ord), ]
stopifnot(identical(d$Cohort, order_cohorts))  # 行序校验（不符即报错停止）

labs_top_down <- c("TCGA overall (reference)",
                   "TCGA oral cavity subset",
                   "GSE41613 (100% oral)",
                   "GSE42743 overall (96% oral)",
                   "GSE42743 oral cavity subset",
                   "GSE65858 oral cavity subset")
d$ylab <- paste0(labs_top_down[d$ord], "\n(N=", d$N, "/", d$Events, " events)")
d$yf <- factor(d$ylab, levels = rev(d$ylab))  # factor 第一层在底部
d$is_tcga <- grepl("TCGA", d$Cohort)
d$sig <- ifelse(d$P < 0.05, "*", "")
d$txt <- sprintf("%.3f (%.2f\u2013%.2f)%s", d$HR, d$LCL, d$UCL, d$sig)

TXT_X <- 1.72   # 数值列固定 x（max UCL = 1.614 < 1.72，零重叠）
p <- ggplot(d, aes(x = HR, y = yf)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
  geom_errorbarh(aes(xmin = LCL, xmax = UCL, color = is_tcga),
                 height = 0.22, linewidth = 0.95) +
  geom_point(aes(color = is_tcga), size = 2.6) +
  geom_text(aes(x = TXT_X, y = yf, label = txt), hjust = 0, size = 3.0,
            color = "grey15") +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"), guide = "none") +
  scale_x_continuous(limits = c(0.35, 2.08), breaks = seq(0.5, 2.0, 0.5)) +
  labs(x = "Per-1-SD HR (95% CI), fixed 67-gene common program",
       y = NULL,
       title = "Oral-cavity external comparison (exploratory)",
       subtitle = "TCGA oral HR 1.34 (P=0.0019) not reproduced in oral-enriched external cohorts\n* Nominal P<0.05; estimates were directionally inconsistent across cohorts; anatomical matching did not restore transportability") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        axis.text.y = element_text(size = 9),
        plot.margin = margin(6, 10, 6, 6))

ggsave(file.path(OUT, "FigS13_oral_external_67g_CANONICAL_20260927.png"), p,
       width = 8, height = 4.6, dpi = 300)
cat("FigS13 fixed version saved\n")
print(d[, c("Cohort", "N", "Events", "HR", "LCL", "UCL", "P")], digits = 4, row.names = FALSE)
