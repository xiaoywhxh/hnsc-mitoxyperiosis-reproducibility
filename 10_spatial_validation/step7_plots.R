# Task #46 -- Spatial plots for the pre-specified 67-gene UCell validation on GSE208253.
#
# Plots produced (per-sample, using TRUE Visium array coordinates -- NOT pixel coords):
#   A. fig_spatial_ucell_grid.{png,pdf}   : 12-panel UCell score maps (Primary)
#   B. fig_spatial_region_grid.{png,pdf}  : 12-panel author region_4class overlay
#   C. fig_spatial_pathology_grid.{png,pdf}: 12-panel pathologist SCC vs non-SCC overlay
#   D. fig_score_distribution.{png,pdf}   : per-sample UCell violin + library-size violin
#   E. fig_moran_bar.{png,pdf}            : Moran's I per sample (with EI and FDR stars)
#   F. fig_effect_forest.{png,pdf}        : region + pathology effect-size forest (12 samples)
#
# No DEG / CellChat / trajectory / ML / drug prediction is performed here.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(cowplot)
})

IN   <- "D:/GSE208253_spatial/spatial_analysis/input"
SC   <- "D:/GSE208253_spatial/spatial_analysis/score"
MOR  <- "D:/GSE208253_spatial/spatial_analysis/moran"
EFF  <- "D:/GSE208253_spatial/spatial_analysis/effects"
PLOT <- "D:/GSE208253_spatial/spatial_analysis/plots"
dir.create(PLOT, showWarnings = FALSE, recursive = TRUE)

spots  <- fread(file.path(IN, "spot_table_all.tsv"), sep = "\t", header = TRUE)
scores <- fread(file.path(SC, "spot_scores.tsv"), sep = "\t", header = TRUE)
d <- merge(spots, scores, by = c("sample", "barcode"), all.x = TRUE)

# Visium hex array -> plotting cartesian coords.
# array_row = 0..77 ; array_col = 0..127 (odd columns only). Convert to a
# continuous hex layout: x = array_col + (array_row %% 2) * 0.5 ; y = array_row.
d[, px := array_col + (array_row %% 2) * 0.5]
d[, py := -array_row]                     # invert so row 0 is at the top
samples <- paste0("s", 1:12)
d[, sample := factor(sample, levels = samples)]

theme_sp <- theme_void(base_size = 8) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 7.5, hjust = 0.5, colour = "grey30"),
        legend.position = "right",
        legend.key.width = unit(7, "pt"),
        legend.key.height = unit(28, "pt"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6.5),
        plot.margin = margin(2, 2, 2, 2))

# ---------------------------------------------------------------- A. UCell grid
mor <- fread(file.path(MOR, "morans_I_by_sample.tsv"), sep = "\t", header = TRUE)
setnames(mor, "fdr_bh", "fdr", skip_absent = TRUE)

sub_lab <- function(s, stat_dt) {
  r <- stat_dt[sample == s]
  if (nrow(r) == 0) return(NA_character_)
  sprintf("I=%.3f  FDR=%.1e", r$morans_I[1], r$fdr[1])
}

plots_ucell <- lapply(samples, function(s) {
  ds <- d[sample == s]
  ggplot(ds, aes(px, py, fill = ucell)) +
    geom_point(shape = 15, size = 1.15) +
    scale_fill_viridis_c(option = "magma", name = "UCell",
                         limits = c(0, 1), oob = scales::squish) +
    coord_fixed() +
    labs(title = s, subtitle = sub_lab(s, mor)) +
    theme_sp
})
pA <- wrap_plots(plots_ucell, ncol = 4) +
  plot_annotation(
    title = "GSE208253 Visium - pre-specified 67-gene UCell score (Primary method)",
    subtitle = "Per-sample independent rank-based scoring; all 67 genes retained; 12 HPV-negative OSCC samples",
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30")))
ggsave(file.path(PLOT, "fig_spatial_ucell_grid.png"), pA, width = 13, height = 9.5, dpi = 320)
ggsave(file.path(PLOT, "fig_spatial_ucell_grid.pdf"), pA, width = 13, height = 9.5)

# ------------------------------------------------------------- B. region overlay
region_cols <- c(core = "#B2182B", edge = "#EF8A62",
                 transitory = "#67A9CF", nc = "#2166AC")
reg_lv <- c("core", "edge", "transitory", "nc")
d_reg <- d[!is.na(region_4class) & region_4class != ""]
d_reg[, region_4class := factor(region_4class, levels = reg_lv)]
# spots without author annotation (incl. the 477 missing in s3) drawn in grey underneath
d_un <- d[is.na(region_4class) | region_4class == ""]
plots_reg <- lapply(samples, function(s) {
  ds <- d_reg[sample == s]; du <- d_un[sample == s]
  ggplot() +
    geom_point(data = du, aes(px, py), shape = 15, size = 1.15, colour = "grey85") +
    geom_point(data = ds, aes(px, py, colour = region_4class), shape = 15, size = 1.15) +
    scale_colour_manual(values = region_cols, drop = FALSE, name = "region_4class") +
    coord_fixed() +
    labs(title = s, subtitle = sprintf("%d annotated spots", nrow(ds))) +
    theme_sp
})
pB <- wrap_plots(plots_reg, ncol = 4) +
  plot_annotation(
    title = "GSE208253 Visium - author region_4class annotation (Arora et al.)",
    subtitle = "Grey = in-tissue spots without author annotation (not imputed, not analysed)",
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9, colour = "grey30")))
ggsave(file.path(PLOT, "fig_spatial_region_grid.png"), pB, width = 13, height = 9.5, dpi = 320)
ggsave(file.path(PLOT, "fig_spatial_region_grid.pdf"), pB, width = 13, height = 9.5)

# ----------------------------------------------------- C. pathology SCC vs non-SCC
NON_SCC <- c("Lymphocyte Negative Stroma", "Lymphocyte Positive Stroma", "Muscle",
             "Glandular Stroma", "Non-cancerous Mucosa", "Lymphocyte Positive Muscles",
             "Artery/Vein")
EXCLUDED <- c("Artifact", "Cautery", "Fold", "Edge Effects", "Keratin")
d_path <- copy(d)
d_path[, path_class := fifelse(pathologist_anno_raw == "SCC", "SCC",
                        fifelse(pathologist_anno_raw %in% NON_SCC, "non-SCC",
                        fifelse(pathologist_anno_raw %in% EXCLUDED, "excluded", "other/NA")))]
d_path[, path_class := factor(path_class, levels = c("SCC", "non-SCC", "excluded", "other/NA"))]
path_cols <- c("SCC" = "#B2182B", "non-SCC" = "#2166AC",
               "excluded" = "grey55", "other/NA" = "grey85")
plots_path <- lapply(samples, function(s) {
  ds <- d_path[sample == s]
  # NB: counts must be computed OUTSIDE ggplot() -- `ds` is not in scope inside
  # the ggplot() call, so referencing ds$... there silently yields NA.
  n_scc <- sum(ds$path_class == "SCC"); n_non <- sum(ds$path_class == "non-SCC")
  ggplot(ds, aes(px, py, colour = path_class)) +
    geom_point(shape = 15, size = 1.15) +
    scale_colour_manual(values = path_cols, drop = FALSE, name = "pathologist") +
    coord_fixed() +
    labs(title = s, subtitle = sprintf("SCC=%d ; non-SCC=%d", n_scc, n_non)) +
    theme_sp
})
pC <- wrap_plots(plots_path, ncol = 4) +
  plot_annotation(
    title = "GSE208253 Visium - pathologist annotation (independent malignant-localisation check)",
    subtitle = paste0("non-SCC comparator PRE-DEFINED: ", paste(NON_SCC, collapse = "; "),
                      "  |  excluded: ", paste(EXCLUDED, collapse = "; ")),
    theme = theme(plot.title = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 7.5, colour = "grey30")))
ggsave(file.path(PLOT, "fig_spatial_pathology_grid.png"), pC, width = 13, height = 9.5, dpi = 320)
ggsave(file.path(PLOT, "fig_spatial_pathology_grid.pdf"), pC, width = 13, height = 9.5)

# ------------------------------------------------------- D. per-sample distributions
dl <- melt(d[!is.na(ucell)], id.vars = c("sample"), measure.vars = c("ucell"),
           variable.name = "metric", value.name = "value")
pD1 <- ggplot(dl, aes(sample, value)) +
  geom_violin(fill = "#F2A65A", colour = "grey30", alpha = .85, width = .85) +
  geom_boxplot(width = .12, outlier.size = .2) +
  labs(title = "Per-sample UCell score distribution",
       subtitle = "Primary method; 67 genes; per-sample independent ranking",
       x = NULL, y = "UCell score") + ylim(0, 1) +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold"))
pD2 <- ggplot(d, aes(sample, library_size)) +
  geom_violin(fill = "#7FB3D5", colour = "grey30", alpha = .85, width = .85) +
  geom_boxplot(width = .12, outlier.size = .2) +
  scale_y_log10() +
  labs(title = "Per-sample library size (total counts / spot)",
       subtitle = "Shown as a diagnostic: depth varies strongly between samples and regions",
       x = NULL, y = "UMI (log10)") +
  theme_bw(base_size = 9) + theme(plot.title = element_text(face = "bold"))
pD <- pD1 / pD2
ggsave(file.path(PLOT, "fig_score_distribution.png"), pD, width = 7.2, height = 7.6, dpi = 320)
ggsave(file.path(PLOT, "fig_score_distribution.pdf"), pD, width = 7.2, height = 7.6)

# ---------------------------------------------------------------- E. Moran's I bar
ei <- mor[, .(sample, morans_I, expected_I, fdr)]
ei[, sample := factor(sample, levels = samples)]     # numeric order s1..s12
el <- melt(ei, id.vars = c("sample", "fdr"), variable.name = "type", value.name = "I")
pE <- ggplot(el, aes(sample, I, fill = type)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = .3) +
  geom_text(data = ei, aes(sample, morans_I, label = ifelse(fdr < 0.05, "*", "ns")),
            inherit.aes = FALSE, vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(morans_I = "#4DBBD5", expected_I = "grey70"),
                    labels = c("Observed Moran's I", "E[I] = -1/(n-1)")) +
  labs(title = "Global Moran's I of the 67-gene UCell score, per sample",
       subtitle = "* FDR < 0.05 (BH within sample); hex adjacency, row-standardised W",
       x = NULL, y = "Moran's I", fill = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")
ggsave(file.path(PLOT, "fig_moran_bar.png"), pE, width = 7.2, height = 4.0, dpi = 320)
ggsave(file.path(PLOT, "fig_moran_bar.pdf"), pE, width = 7.2, height = 4.0)

# ----------------------------------------------------------------- F. effect forest
reff <- fread(file.path(EFF, "region_effect_sizes.tsv"), sep = "\t")
peff <- fread(file.path(EFF, "pathology_effect_sizes.tsv"), sep = "\t")
reff <- reff[!is.na(cliffs_delta)]
peff <- peff[!is.na(cliffs_delta)]

fr <- rbind(
  reff[comparison == "core_vs_nc",       .(sample, delta = cliffs_delta, grp = "core vs nc")],
  reff[comparison == "edge_vs_nc",       .(sample, delta = cliffs_delta, grp = "edge vs nc")],
  reff[comparison == "transitory_vs_nc", .(sample, delta = cliffs_delta, grp = "transitory vs nc")],
  peff[, .(sample, delta = cliffs_delta, grp = "SCC vs non-SCC")])
fr[, sample := factor(sample, levels = samples)]
pF <- ggplot(fr, aes(delta, sample, colour = delta > 0)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 1.8) +
  geom_segment(aes(x = 0, xend = delta, yend = sample), linewidth = .5) +
  facet_wrap(~grp, nrow = 1) +
  scale_colour_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"), guide = "none") +
  labs(title = "Per-sample effect sizes (Cliff's delta) - no pooled spot-level inference",
       subtitle = "Red = positive (program higher in group A); Blue = negative. Cross-sample summary = sign test + random-effects meta",
       x = "Cliff's delta", y = NULL) +
  theme_bw(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"))
ggsave(file.path(PLOT, "fig_effect_forest.png"), pF, width = 10.5, height = 4.2, dpi = 320)
ggsave(file.path(PLOT, "fig_effect_forest.pdf"), pF, width = 10.5, height = 4.2)

cat("Plots written to:", PLOT, "\n")
print(list.files(PLOT))
