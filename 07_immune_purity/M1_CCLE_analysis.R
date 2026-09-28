## ============================================================
## M1 (CCLE 版): 细胞系 program score + HNSCC 亚组分析
## 数据源: GSE36133 (CCLE Affymetrix 表达) + 细胞系注释
## 注: DepMap Chronos 依赖分数当前网络受限（figshare 403），
##     本模块为 CCLE 表达层面的虚拟敲低代理（PRKN/HK1 高低表达分组）
## ============================================================
suppressMessages({library(hgu133plus2.db); library(org.Hs.eg.db)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data/DepMap"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"

# 1. 解析 GSE36133 series matrix（探针级，注释到基因）
cat(">>> 1. 解析 GSE36133...\n")
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE36133_series_matrix.txt.gz"), open = "rt"), warn = FALSE)
start <- grep("!series_matrix_table_begin", lines)
hdr <- strsplit(lines[start + 1], "\t")[[1]]
samples <- hdr[-1]
# 探针行
probe_rows <- list()
for (i in (start + 2):length(lines)) {
  l <- lines[i]
  if (grepl("!series_matrix_table_end", l)) break
  parts <- strsplit(l, "\t")[[1]]
  probe_rows[[i - start - 1]] <- parts
}
cat("探针行数:", length(probe_rows), " 样本数:", length(samples), "\n")
probes <- sapply(probe_rows, function(p) gsub("\"", "", p[1]))
mat <- do.call(rbind, lapply(probe_rows, function(p) as.numeric(p[-1])))
rownames(mat) <- probes
colnames(mat) <- samples

# 2. 探针 → gene symbol（hgu133plus2.db，max median 聚合）
cat(">>> 2. 探针注释...\n")
sym <- mapIds(hgu133plus2.db, keys = rownames(mat), column = "SYMBOL", keytype = "PROBEID")
mat$gene <- sym[rownames(mat)]
mat <- mat[!is.na(mat$gene) & mat$gene != "", ]
# 聚合：每基因取中位数最大的探针
genes_uniq <- unique(mat$gene)
expr_g <- matrix(NA, nrow = length(genes_uniq), ncol = ncol(mat) - 1,
                 dimnames = list(genes_uniq, colnames(mat)[1:(ncol(mat)-1)]))
for (g in genes_uniq) {
  sub <- mat[mat$gene == g, 1:(ncol(mat)-1), drop = FALSE]
  med <- apply(sub, 1, median, na.rm = TRUE)
  expr_g[g, ] <- as.numeric(sub[which.max(med), ])
}
cat("基因级矩阵:", nrow(expr_g), "x", ncol(expr_g), "\n")

# 3. 细胞系注释（GSE36133 的 clinical 数据或标题解析）
# series matrix 的 !Sample_title 通常是细胞系名；临床特征行有 lineage
title_line <- grep("!Sample_title", lines)[1]
titles <- strsplit(lines[title_line], "\t")[[1]][-1]
titles <- gsub("\"", "", titles)
cat(">>> 3. 细胞系标题前5:", head(titles, 5), "\n")

# 4. program score（67-gene）
common67 <- read.csv("D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation/program_common_genes_3OS.csv")$common_gene
g67 <- intersect(common67, rownames(expr_g))
cat("CCLE 67-gene 覆盖:", length(g67), "/", length(common67), "\n")
z <- t(scale(t(expr_g[g67, ])))
prog <- colMeans(z, na.rm = TRUE)

# 5. PRKN/HK1 表达
prkn_expr <- expr_g["PRKN", ]
hk1_expr <- expr_g["HK1", ]

# 6. 输出
res <- data.frame(cell_line = names(prog), program_score = prog,
                  PRKN = prkn_expr[names(prog)], HK1 = hk1_expr[names(prog)])
write.csv(res, file.path(OUT_DIR, "M1_CCLE_program_scores.csv"), row.names = FALSE)
cat(">>> 已保存 M1_CCLE_program_scores.csv, n =", nrow(res), "\n")
cat(">>> program score 范围:", round(range(prog), 3), "\n")
cat(">>> PRKN 表达与 program 相关:", round(cor(prkn_expr, prog, method = "spearman", use = "complete.obs"), 3), "\n")
cat(">>> HK1 表达与 program 相关:", round(cor(hk1_expr, prog, method = "spearman", use = "complete.obs"), 3), "\n")
cat("DONE\n")
