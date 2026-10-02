suppressPackageStartupMessages(library(data.table))
spots <- fread("input/spot_table_all.tsv", sep="\t", header=TRUE)
d <- merge(spots, fread("score/spot_scores.tsv", sep="\t"), by=c("sample","barcode"))
# distribution check + relation to library size / detected genes
s <- "s1"; dd <- d[sample==s]
cat("UCell vs n_gene_detected spearman:", cor(dd$ucell, dd$n_gene_detected, method="spearman"), "\n")
cat("meanZ vs n_gene_detected spearman:", cor(dd$meanz, dd$n_gene_detected, method="spearman"), "\n")
cat("UCell vs library_size spearman:", cor(dd$ucell, dd$library_size, method="spearman"), "\n")
cat("meanZ vs library_size spearman:", cor(dd$meanz, dd$library_size, method="spearman"), "\n")
cat("UCell range", range(dd$ucell), " meanZ range", range(dd$meanz), "\n")
# top genes per spot check: which genes drive meanZ
cat("\nUCell quartile vs meanZ median:/n")
q <- cut(dd$ucell, quantile(dd$ucell, 0:4/4), include.lowest=TRUE)
print(tapply(dd$meanz, q, median))
q2 <- cut(dd$meanz, quantile(dd$meanz, 0:4/4), include.lowest=TRUE)
cat("\nmeanZ quartile vs UCell median:/n")
print(tapply(dd$ucell, q2, median))
