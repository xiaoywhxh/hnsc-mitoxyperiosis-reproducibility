# Depth-sensitivity diagnostics (robustness). Does the core-vs-nc direction
# survive conditioning on sequencing depth / detected-gene count?
suppressPackageStartupMessages({ library(data.table); library(metafor) })
IN  <- "D:/GSE208253_spatial/spatial_analysis/input"
SC  <- "D:/GSE208253_spatial/spatial_analysis/score"
OUT <- "D:/GSE208253_spatial/spatial_analysis/robustness"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

spots <- fread(file.path(IN,"spot_table_all.tsv"), sep="\t", header=TRUE)
d <- merge(spots, fread(file.path(SC,"spot_scores.tsv"), sep="\t"), by=c("sample","barcode"))

cliff <- function(a,b){ if(!length(a)||!length(b)) return(NA_real_)
  r<-rank(c(a,b)); U<-sum(r[seq_along(a)])-length(a)*(length(a)+1)/2
  2*U/(length(a)*length(b))-1 }

# depth-adjusted UCell: residual of ucell ~ log(library) + log(n_gene_detected), per sample
d <- d[library_size > 0 & n_gene_detected > 0]
d[, llib := log(library_size)]
d[, lnd  := log(n_gene_detected)]
d[, ucell_adj := { fit <- lm(ucell ~ llib + lnd); resid(fit) }, by = sample]

rows <- list()
for (cn in c("core_vs_nc","edge_vs_nc","transitory_vs_nc")) {
  A <- sub("_vs_.*","",cn); B <- "nc"
  for (s in paste0("s",1:12)) {
    da <- d[sample==s & region_4class==A]; db <- d[sample==s & region_4class==B]
    if (nrow(da)<3||nrow(db)<3) next
    rows[[paste(cn,s)]] <- data.table(comparison=cn, sample=s,
      d_raw = cliff(da$ucell, db$ucell), d_adj = cliff(da$ucell_adj, db$ucell_adj))
  }
}
r <- rbindlist(rows)
fwrite(r, file.path(OUT,"depth_adjusted_region_direction.tsv"), sep="\t")
cat("=== Depth-adjusted UCell: direction per comparison ===\n")
print(r[, .(n=.N, n_pos_raw=sum(d_raw>0), n_pos_adj=sum(d_adj>0),
            mean_raw=round(mean(d_raw),3), mean_adj=round(mean(d_adj),3)), by=comparison])

# meta on adjusted
cat("\n=== RE meta on depth-adjusted delta ===\n")
for (cn in unique(r$comparison)) {
  sub <- as.data.frame(r[comparison==cn])
  sub$se <- sqrt((2 + 1) / (3 * 10 * 10)) * 2   # conservative ~0.1 placeholder removed; use empirical
  # empirical: use raw delta's observed spread vs n (approx). Instead run simple t-test on deltas.
  tt <- t.test(sub$d_adj)
  cat(sprintf("%-16s adj_delta_mean=%.3f  95%%CI[%.3f,%.3f]  t-test p=%.3g  k=%d\n",
              cn, mean(sub$d_adj), tt$conf.int[1], tt$conf.int[2], tt$p.value, nrow(sub)))
}

# also: is depth itself different between core and nc? (confounding check)
cat("\n=== depth difference core vs nc (median) ===\n")
dd <- d[region_4class %in% c("core","nc")]
print(dd[, .(med_lib_core = as.numeric(median(library_size[region_4class=="core"])),
             med_lib_nc   = as.numeric(median(library_size[region_4class=="nc"])),
             med_nd_core  = as.numeric(median(n_gene_detected[region_4class=="core"])),
             med_nd_nc    = as.numeric(median(n_gene_detected[region_4class=="nc"]))), by=sample])
