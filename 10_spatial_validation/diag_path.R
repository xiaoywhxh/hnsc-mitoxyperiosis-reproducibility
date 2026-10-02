suppressPackageStartupMessages(library(data.table))
spots <- fread("input/spot_table_all.tsv", sep="\t", header=TRUE)
d <- spots[pathologist_anno_raw != ""]
print(d[, .N, by = pathologist_anno_raw][order(-N)])
cat("\n=== which categories are unambiguously NON-SCC tissue? ===\n")
cat("total annotated:", nrow(d), "\n")
