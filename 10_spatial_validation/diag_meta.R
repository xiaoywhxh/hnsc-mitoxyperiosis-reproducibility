suppressPackageStartupMessages({library(data.table); library(metafor)})
eff <- fread("effects/region_effect_sizes.tsv", sep="\t", header=TRUE)
sub <- eff[comparison=="core_vs_nc" & !is.na(cliffs_delta)]
sub[, se := sqrt((n_A + n_B + 1) / (3 * n_A * n_B))]
print(sub[, .(sample, cliffs_delta, se)])
m <- rma(yi = cliffs_delta, sei = se, data = as.data.frame(sub), method = "REML")
print(m)
