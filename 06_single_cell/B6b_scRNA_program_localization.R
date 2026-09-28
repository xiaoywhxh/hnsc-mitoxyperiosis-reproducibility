## (b) 单细胞定位分析 — GSE103322 HNSCC scRNA-seq
## mitoxyperiosis-related program score 按细胞类型定位
## 方法: gene-level Z-score 均值（UCell/AUCell 未安装；均值 Z 为程序富集的稳健代理）

suppressMessages({
  library(ggplot2)
  library(dplyr)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/11_singlecell"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 读取表达矩阵（log 尺度，行 6+ 为基因）
# ============================================================
cat(">>> 读取 GSE103322 单细胞数据...\n")
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"), open = "rt"), warn = FALSE)

# 元数据行 1-6 (0-based: 行0=cell barcode, 行1-5=注释)
cell_barcodes <- strsplit(lines[0 + 1], "\t")[[1]][-1]
meta_line5 <- strsplit(lines[5 + 1], "\t")[[1]][-1]  # non-cancer cell type

# 细胞注释（来自 cell_annotations.Rds 更准确）
annot <- readRDS(file.path(DATA_DIR, "GSE103322_cell_annotations.Rds"))
cat("注释细胞数:", nrow(annot), "\n")
cat("cell_type 分布:\n")
t <- table(annot$cell_type)
for (i in seq_along(t)) cat("  ", names(t)[i], ":", t[i], "\n")

# 基因表达矩阵（行 6+）
gene_start <- 6  # 0-based 行6 = 'C9orf152'
n_cells <- length(cell_barcodes)
n_genes <- length(lines) - gene_start
cat("基因数:", n_genes, " 细胞数:", n_cells, "\n")

gene_names <- character(n_genes)
expr_mat <- matrix(NA, nrow = n_genes, ncol = n_cells)
for (i in 1:n_genes) {
  parts <- strsplit(lines[gene_start + i], "\t")[[1]]
  gene_names[i] <- gsub("'", "", parts[1])
  expr_mat[i, ] <- as.numeric(parts[2:(n_cells + 1)])
}
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- cell_barcodes
cat("表达矩阵:", nrow(expr_mat), "x", ncol(expr_mat), "\n")

# ============================================================
# 2. 基因集定义（canonical 73）
# ============================================================
core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
all73 <- unique(c(core_genes, mtor_genes, mito_genes, metab_genes))

modules <- list(
  Core = core_genes, mTOR = mtor_genes,
  MitoDynamics = mito_genes, MetaboImmune = metab_genes,
  All = all73
)

# 基因覆盖
found <- intersect(all73, rownames(expr_mat))
cat("\n73 基因在 scRNA 覆盖:", length(found), "/73\n")
cat("缺失:", setdiff(all73, rownames(expr_mat)), "\n")

# ============================================================
# 3. Program score（基因级 Z 均值）
# ============================================================
z_mat <- t(scale(t(expr_mat[found, ])))  # gene-level Z
program_score <- colMeans(z_mat, na.rm = TRUE)
names(program_score) <- colnames(expr_mat)

# 各模块
module_scores <- list()
for (m in names(modules)) {
  genes_m <- intersect(modules[[m]], rownames(expr_mat))
  module_scores[[m]] <- colMeans(z_mat[genes_m, , drop = FALSE], na.rm = TRUE)
}

# ============================================================
# 4. 按细胞类型定位
# ============================================================
# 匹配注释: 注释 cell_id 与表达列名格式一致
annot$cell_id_clean <- gsub("-", ".", annot$cell_id)
common <- intersect(names(program_score), annot$cell_id_clean)
cat("\n匹配细胞:", length(common), "\n")

score_df <- data.frame(
  Cell = names(program_score),
  Program_All = as.numeric(program_score),
  Core = as.numeric(module_scores$Core),
  mTOR = as.numeric(module_scores$mTOR),
  MitoDynamics = as.numeric(module_scores$MitoDynamics),
  MetaboImmune = as.numeric(module_scores$MetaboImmune),
  stringsAsFactors = FALSE
)
annot_sub <- annot[match(score_df$Cell, annot$cell_id_clean), ]
score_df$CellType <- annot_sub$cell_type
score_df <- score_df[!is.na(score_df$CellType), ]
cat("有注释细胞:", nrow(score_df), "\n")
cat("细胞类型分布:\n")
print(table(score_df$CellType))

# ============================================================
# 5. 统计比较 + 图
# ============================================================
# Kruskal-Wallis 检验
kw <- kruskal.test(Program_All ~ CellType, data = score_df)
cat(sprintf("\nKruskal-Wallis (overall program): P=%.3g\n", kw$p.value))

# 各细胞类型中位数
med <- tapply(score_df$Program_All, score_df$CellType, median)
cat("\nOverall program 各细胞类型中位数:\n")
for (i in seq_along(med)) cat(sprintf("  %s: %.4f\n", names(med)[i], med[i]))

# 保存
write.csv(score_df, file.path(OUT_DIR, "GSE103322_program_scores_by_celltype.csv"), row.names = FALSE)

# 图: overall program by cell type
p1 <- ggplot(score_df, aes(x = reorder(CellType, Program_All, median), y = Program_All, fill = CellType)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  scale_fill_manual(values = c("Malignant" = "#E64B35", "Fibroblast" = "#4DBBD5",
                               "T cell" = "#3C5488", "B cell" = "#F39B7F",
                               "Myeloid/Macro" = "#00A087", "Mast" = "#8491B4",
                               "Other" = "#B09C85")) +
  labs(title = "Mitoxyperiosis-related program score by cell type (GSE103322)",
       subtitle = "HNSCC single-cell atlas; overall program (mean Z)",
       x = "", y = "Program score (mean Z)") +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "Fig_scRNA_program_by_celltype.pdf"), p1, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "Fig_scRNA_program_by_celltype.png"), p1, width = 7, height = 5, dpi = 300)

# 模块热图（各细胞类型均值）
mod_means <- data.frame(CellType = score_df$CellType,
                        Core = score_df$Core, mTOR = score_df$mTOR,
                        Mito = score_df$MitoDynamics, MetaboImmune = score_df$MetaboImmune)
mod_summary <- aggregate(cbind(Core, mTOR, Mito, MetaboImmune) ~ CellType, data = mod_means, FUN = mean)
write.csv(mod_summary, file.path(OUT_DIR, "GSE103322_module_scores_by_celltype.csv"), row.names = FALSE)
cat("\n模块均值 by cell type:\n")
print(mod_summary)

# 热图
mod_long <- reshape2::melt(mod_summary, id.vars = "CellType", variable.name = "Module", value.name = "Mean")
p2 <- ggplot(mod_long, aes(x = Module, y = CellType, fill = Mean)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35", midpoint = 0) +
  labs(title = "Module-level program scores by cell type",
       x = "", y = "") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUT_DIR, "Fig_scRNA_module_heatmap.pdf"), p2, width = 6, height = 4)
ggsave(file.path(OUT_DIR, "Fig_scRNA_module_heatmap.png"), p2, width = 6, height = 4, dpi = 300)

cat("\n=== 单细胞定位分析 DONE ===")
