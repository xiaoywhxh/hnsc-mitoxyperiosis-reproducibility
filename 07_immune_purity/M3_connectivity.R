## M3: LINCS HK1 敲低共识签名 × TCGA program 相关签名 连接性分析
## 数据: consensi-knockdown.tsv.bz2 (dhimmel/lincs, HK1 Entrez 3098 行) + TCGA-HNSC_full.txt
suppressMessages({
  library(data.table)
  library(org.Hs.eg.db)
})
OUT <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"

# ========== 1. TCGA program score + 全转录组 program 相关签名 ==========
cat(">>> 1. 读 TCGA 矩阵...\n")
expr <- fread("D:/HNSC_mitoxyperiosis_positron/data/TCGA_full/TCGA-HNSC_full.txt", header=TRUE)
genes <- expr$gene_id
mat <- as.matrix(expr[, -1, with=FALSE]); rownames(mat) <- genes
tumor <- grep("-01", colnames(mat)); mat <- mat[, tumor]
cat("肿瘤样本:", ncol(mat), " 基因:", nrow(mat), "\n")

common67 <- read.csv("D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation/program_common_genes_3OS.csv")$common_gene
g67 <- intersect(common67, rownames(mat))
cat("67-gene 覆盖:", length(g67), "/", length(common67), "\n")
z <- t(scale(t(mat[g67, , drop=FALSE])))
prog <- colMeans(z, na.rm=TRUE)

cat(">>> 2. 全转录组 × program Spearman (60k 基因, 约 1-2 分钟)...\n")
scores <- apply(mat, 1, function(x) {
  ok <- !is.na(x) & is.finite(x)
  if (sum(ok) < 50) return(NA)
  suppressWarnings(cor(x, prog, method="spearman", use="complete.obs"))
})
cat("program 相关签名基因数:", sum(!is.na(scores)), "\n")

# ========== 2. LINCS HK1 敲低共识签名 (Entrez 3098) ==========
cat(">>> 3. 读 LINCS consensi-knockdown HK1(3098) 签名...\n")
library(R.utils)
con <- bzfile("D:/HNSC_mitoxyperiosis_positron/data/DepMap/consensi-knockdown.tsv.bz2", "rt")
hdr <- strsplit(readLines(con, n=1), "\t")[[1]]
target <- NULL
repeat {
  line <- readLines(con, n=1)
  if (length(line) == 0) break
  parts <- strsplit(line, "\t")[[1]]
  if (parts[1] == "3098") { target <- parts; break }
}
close(con)
stopifnot(!is.null(target))
entrez_cols <- as.character(hdr[-1])
hk1_sig <- setNames(as.numeric(target[-1]), entrez_cols)
cat("HK1 敲低签名基因数:", length(hk1_sig), "\n")

# Entrez -> Symbol 映射
sym <- mapIds(org.Hs.eg.db, keys=names(hk1_sig), column="SYMBOL", keytype="ENTREZID")
keep <- !is.na(sym) & !duplicated(sym)
hk1_sig2 <- hk1_sig[keep]; names(hk1_sig2) <- sym[keep]
cat("映射后 HK1 签名基因数:", length(hk1_sig2), "\n")

# ========== 3. 连接性 ==========
common <- intersect(names(scores)[!is.na(scores)], names(hk1_sig2))
cat("共同基因:", length(common), "\n")
x <- scores[common]   # TCGA program 相关 ρ
y <- hk1_sig2[common] # HK1 敲低 z-score
r_conn <- cor(x, y, method="spearman")
n <- length(common)
t_stat <- r_conn * sqrt((n-2)/(1-r_conn^2))
p_conn <- 2 * (1 - pnorm(abs(t_stat)))
cat("\n========================================\n")
cat("连接性: HK1 敲低签名 × program 相关签名\n")
cat("n =", n, " Spearman rho =", round(r_conn, 4), " P =", format.pval(p_conn, digits=3), "\n")
cat("正相关 => program 高 ≈ HK1 敲低(低表达)转录状态\n")
cat("========================================\n")

# 反向对照: program 高 ≈ HK1 高表达? (y 取反)
r_rev <- cor(x, -y, method="spearman")
cat("反向(HK1 高表达) rho =", round(r_rev, 4), "\n")

# 保存
out_df <- data.frame(gene=common, program_rho=x, HK1_KD_z=y)
fwrite(out_df, file.path(OUT, "M3_connectivity_data.csv"))
cat("\n已保存 M3_connectivity_data.csv\n")
