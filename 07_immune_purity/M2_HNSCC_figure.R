## M2 HNSCC-only 亚组散点图：program score × Phenformin/Daporinad IC50
## 数据: M2_phenformin_daporinad_data.csv 中 UPPER_AERODIGESTIVE_TRACT 细胞系
suppressMessages({library(ggplot2)})
OUT_DIR <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"
d <- read.csv(file.path(OUT_DIR, "M2_phenformin_daporinad_data.csv"), check.names = FALSE)
hn <- d[grepl("UPPER_AERODIGESTIVE_TRACT", d$cell_line), ]
hn$Phenformin_IC50 <- as.numeric(hn$Phenformin_IC50)
hn$Daporinad_IC50 <- as.numeric(hn$Daporinad_IC50)

mk <- function(x, y, label, rho, pval) {
  dd <- data.frame(x = x, y = log10(y))
  dd <- dd[complete.cases(dd), ]
  p <- ggplot(dd, aes(x = x, y = y)) +
    geom_point(color = "#3C5488", alpha = 0.75, size = 2.2) +
    geom_smooth(method = "lm", se = TRUE, color = "#E64B35", linewidth = 0.8) +
    labs(x = "Mitoxyperiosis program score (67-gene, mean-Z)",
         y = label,
         title = sprintf("%s  HNSCC-only (Spearman rho = %.3f, P = %.3g)", label, rho, pval),
         subtitle = "CCLE UPPER_AERODIGESTIVE_TRACT cell lines") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(size = 10.5))
  p
}

rho_ph <- cor(hn$program_score, hn$Phenformin_IC50, method = "spearman", use = "complete.obs")
rho_da <- cor(hn$program_score, hn$Daporinad_IC50, method = "spearman", use = "complete.obs")
n_ph <- sum(!is.na(hn$Phenformin_IC50)); n_da <- sum(!is.na(hn$Daporinad_IC50))
p_ph <- 2 * (1 - pnorm(abs(rho_ph) * sqrt((n_ph - 2) / (1 - rho_ph^2))))
p_da <- 2 * (1 - pnorm(abs(rho_da) * sqrt((n_da - 2) / (1 - rho_da^2))))

p1 <- mk(hn$program_score, hn$Phenformin_IC50, "Phenformin IC50 (log10 nM)", rho_ph, p_ph)
p2 <- mk(hn$program_score, hn$Daporinad_IC50, "Daporinad IC50 (log10 nM)", rho_da, p_da)
ggsave(file.path(OUT_DIR, "M2_HNSCC_phenformin_scatter.png"), p1, width = 5.2, height = 4.2, dpi = 300)
ggsave(file.path(OUT_DIR, "M2_HNSCC_daporinad_scatter.png"), p2, width = 5.2, height = 4.2, dpi = 300)
cat(sprintf("HNSCC-only: Phenformin n=%d rho=%.3f P=%.4g | Daporinad n=%d rho=%.3f P=%.4g\n",
            n_ph, rho_ph, p_ph, n_da, rho_da, p_da))
