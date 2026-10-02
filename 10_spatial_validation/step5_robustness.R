# ROBUSTNESS ONLY (does NOT select the primary method; Primary remains UCell).
# Compare UCell vs simple standardised mean-expression score:
#   (a) Spearman correlation per sample
#   (b) whether region effect DIRECTION agrees
suppressPackageStartupMessages({ library(data.table); library(metafor) })
IN  <- "D:/GSE208253_spatial/spatial_analysis/input"
SC  <- "D:/GSE208253_spatial/spatial_analysis/score"
OUT <- "D:/GSE208253_spatial/spatial_analysis/robustness"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

spots  <- fread(file.path(IN, "spot_table_all.tsv"), sep = "\t", header = TRUE)
d <- merge(spots, fread(file.path(SC, "spot_scores.tsv"), sep = "\t"), by = c("sample","barcode"))

cliff_delta <- function(a, b) {
  if (length(a) < 1 || length(b) < 1) return(NA_real_)
  r <- rank(c(a, b)); U <- sum(r[seq_along(a)]) - length(a)*(length(a)+1)/2
  2*U/(length(a)*length(b)) - 1
}

# ---- (a) Spearman UCell vs meanZ ----
sp <- d[, .(n = .N, spearman = cor(ucell, meanz, method = "spearman")), by = sample]
# ---- (b) region direction agreement ----
cmp <- list(core_vs_nc = c("core","nc"), edge_vs_nc = c("edge","nc"),
            transitory_vs_nc = c("transitory","nc"))
agr <- list()
for (cn in names(cmp)) {
  A <- cmp[[cn]][1]; B <- cmp[[cn]][2]
  for (s in paste0("s",1:12)) {
    da <- d[sample==s & region_4class==A]; db <- d[sample==s & region_4class==B]
    if (nrow(da)<3 || nrow(db)<3) next
    agr[[paste(cn,s)]] <- data.table(comparison=cn, sample=s,
      delta_ucell = cliff_delta(da$ucell, db$ucell),
      delta_meanz = cliff_delta(da$meanz, db$meanz))
  }
}
ag <- rbindlist(agr)
ag[, agree := sign(delta_ucell) == sign(delta_meanz)]
# pathology direction agreement too
NON_SCC <- c("Lymphocyte Negative Stroma","Lymphocyte Positive Stroma","Muscle",
             "Glandular Stroma","Non-cancerous Mucosa","Lymphocyte Positive Muscles","Artery/Vein")
d[, pclass := fifelse(pathologist_anno_raw=="SCC","SCC",
              fifelse(pathologist_anno_raw %in% NON_SCC,"nonSCC","other"))]
pag <- list()
for (s in paste0("s",1:12)) {
  da <- d[sample==s & pclass=="SCC"]; db <- d[sample==s & pclass=="nonSCC"]
  if (nrow(da)<3 || nrow(db)<3) next
  pag[[s]] <- data.table(comparison="SCC_vs_nonSCC", sample=s,
    delta_ucell = cliff_delta(da$ucell, db$ucell),
    delta_meanz = cliff_delta(da$meanz, db$meanz))
}
pag <- rbindlist(pag); pag[, agree := sign(delta_ucell) == sign(delta_meanz)]
allag <- rbind(ag, pag)

fwrite(sp, file.path(OUT, "robustness_spearman.tsv"), sep="\t")
fwrite(allag, file.path(OUT, "robustness_direction_agreement.tsv"), sep="\t")

cat("=== Spearman(UCell, meanZ) per sample ===\n"); print(sp)
cat(sprintf("\nmean Spearman = %.3f (min %.3f, max %.3f)\n",
            mean(sp$spearman), min(sp$spearman), max(sp$spearman)))
cat("\n=== Region direction agreement (UCell vs meanZ) ===\n")
print(ag[, .(n = .N, n_agree = sum(agree), pct = round(100*mean(agree),1),
             mean_delta_ucell = round(mean(delta_ucell),3),
             mean_delta_meanz = round(mean(delta_meanz),3)), by = comparison])
cat("\n=== Pathology direction agreement ===\n")
print(pag[, .(n=.N, n_agree=sum(agree), mean_delta_ucell=round(mean(delta_ucell),3),
              mean_delta_meanz=round(mean(delta_meanz),3))])
cat(sprintf("\nOverall direction agreement: %d / %d (%.1f%%)\n",
            sum(allag$agree), nrow(allag), 100*mean(allag$agree)))
