suppressPackageStartupMessages(library(data.table))
spots <- fread("D:/GSE208253_spatial/spatial_analysis/input/spot_table_all.tsv", sep="\t", header=TRUE)
for (s in paste0("s", 1:12)) {
  dt <- spots[sample == s]
  lookup <- new.env(hash=TRUE)
  for (i in seq_len(nrow(dt))) assign(paste(dt$array_row[i], dt$array_col[i]), i, envir=lookup)
  offs <- rbind(c(0,-2),c(0,2),c(-1,-1),c(-1,1),c(1,-1),c(1,1))
  z <- 0
  for (i in seq_len(nrow(dt))) {
    r <- dt$array_row[i]; cc <- dt$array_col[i]; n <- 0
    for (k in seq_len(nrow(offs))) {
      v <- mget(paste(r+offs[k,1], cc+offs[k,2]), envir=lookup, ifnotfound=NA)[[1]]
      if (!is.na(v)) n <- n+1
    }
    if (n == 0) z <- z + 1
  }
  cat(sprintf("%s: isolated=%d of %d\n", s, z, nrow(dt)))
}
