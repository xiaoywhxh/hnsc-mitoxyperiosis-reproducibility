# fix_fig_ucell.R -- 修复 fig_spatial_ucell_grid 的渲染缺陷（**不改动任何数值**）。
#
# 缺陷：step7_plots.R 第 57-58 行使用
#         aes(px, py, fill = ucell) + geom_point(shape = 15, ...)
#       ggplot2 中 shape 15 是"实心方块"，只认 colour、不认 fill；
#       因此 UCell 分数被静默忽略，全部 spot 画成默认黑色 → 图变成黑色剪影。
#
# 修法：fill -> colour（scale_fill_viridis_c -> scale_colour_viridis_c），shape 保持 15。
#       **只重新渲染图形；spot_scores.tsv / morans_I_by_sample.tsv 只读，不重算。**
#
# 输出：spatial_analysis/plots/fig_spatial_ucell_grid.{png,pdf}（覆盖）
#       旧版已备份为 fig_spatial_ucell_grid.DEFECTIVE_black.png

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT <- "D:/GSE208253_spatial/spatial_analysis"
IN   <- file.path(ROOT, "input")
SC   <- file.path(ROOT, "score")
MOR  <- file.path(ROOT, "moran")
PLOT <- file.path(ROOT, "plots")
dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)

spots  <- fread(file.path(IN, "spot_table_all.tsv"), sep = "\t", header = TRUE)
scores <- fread(file.path(SC, "spot_scores.tsv"),      sep = "\t", header = TRUE)
d <- merge(spots, scores, by = c("sample", "barcode"), all.x = TRUE)

d[, px := array_col + (array_row %% 2) * 0.5]
d[, py := -array_row]
samples <- paste0("s", 1:12)
d[, sample := factor(sample, levels = samples)]

sha <- function(p) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(file = p, algo = "sha256"))
  return(tools::md5sum(p)[[1]])
}
cat("[frozen-input sha256/md5]\n")
cat("  spot_scores.tsv        :", sha(file.path(SC, "spot_scores.tsv")), "\n")
cat("  morans_I_by_sample.tsv :", sha(file.path(MOR, "morans_I_by_sample.tsv")), "\n")

theme_sp <- theme_void(base_size = 8) +
  theme(plot.title    = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7.5, hjust = 0.5, colour = "grey30"),
        legend.key.height = unit(28, "pt"),
        legend.title = element_text(size = 7),
        legend.text  = element_text(size = 6.5),
        plot.margin  = margin(2, 2, 2, 2))

mor <- fread(file.path(MOR, "morans_I_by_sample.tsv"), sep = "\t", header = TRUE)
setnames(mor, "fdr_bh", "fdr", skip_absent = TRUE)

sub_lab <- function(s, stat_dt) {
  r <- stat_dt[sample == s]
  if (nrow(r) == 0) return(NA_character_)
  sprintf("I=%.3f  FDR=%.1e", r$morans_I[1], r$fdr[1])
}

# ---- 修正点：colour 而非 fill ----
plots_ucell <- lapply(samples, function(s) {
  ds <- d[sample == s]
  ggplot(ds, aes(px, py, colour = ucell)) +
    geom_point(shape = 15, size = 1.15) +
    scale_colour_viridis_c(option = "magma", name = "UCell",
                           limits = c(0, 1), oob = scales::squish) +
    coord_fixed() +
    labs(title = s, subtitle = sub_lab(s, mor)) +
    theme_sp
})

pA <- wrap_plots(plots_ucell, ncol = 4) +
  plot_annotation(
    title = "GSE208253 Visium - pre-specified 67-gene UCell score (Primary method)",
    subtitle = "Per-sample independent rank-based scoring; all 67 genes retained; 12 HPV-negative OSCC samples",
    theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30")))

ggsave(file.path(PLOT, "fig_spatial_ucell_grid.png"), pA, width = 13, height = 9.5, dpi = 320)
ggsave(file.path(PLOT, "fig_spatial_ucell_grid.pdf"), pA, width = 13, height = 9.5)
cat("[OK] regenerated fig_spatial_ucell_grid.{png,pdf}\n")
