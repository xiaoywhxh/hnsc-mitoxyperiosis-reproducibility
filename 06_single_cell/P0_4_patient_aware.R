## P0-4: Single-cell patient-aware re-analysis (审稿意见 2 R1-M6 / R2-M5)
## GSE103322: cell-level KW P<0.001 存在 pseudoreplication → patient-level pseudobulk 复核
## 1) patient × celltype program score 矩阵 2) patient-level Malignant vs non-malignant 配对比较
## 3) mixed-model 近似（patient 为随机效应的 rank-based 检验）

suppressMessages({ library(data.table) })

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/11_singlecell"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. 加载 GSE103322 表达 + 注释
# ============================================================
annot <- readRDS(file.path(DATA_DIR, "GSE103322_cell_annotations.Rds"))  # 5902 cells
cat("Annotations:", nrow(annot), "cells,", length(unique(annot$patient)), "patients\n")

# 表达矩阵（行 6+ 基因, log 尺度; 前 6 行元数据）
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"), open = "rt"), warn = FALSE)
hdr <- strsplit(lines[1], "\t")[[1]]
meta_rows <- 6
cells <- hdr[-1]  # 5901
cat("Expression cells:", length(cells), "\n")

# 基因行
gene_lines <- lines[(meta_rows + 1):length(lines)]
genes <- sapply(gene_lines, function(l) gsub("'", "", strsplit(l, "\t")[[1]][1]))
expr_list <- lapply(gene_lines, function(l) as.numeric(strsplit(l, "\t")[[1]][-1]))
expr_mat <- do.call(rbind, expr_list)
rownames(expr_mat) <- genes
colnames(expr_mat) <- cells
cat("Expr dims:", dim(expr_mat), "\n")

# ============================================================
# 2. canonical 73 gene program score (cell-level mean-Z)
# ============================================================
core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
all73 <- unique(c(core_genes, mtor_genes, mito_genes, metab_genes))
found <- intersect(all73, rownames(expr_mat))
cat("73 genes detected in scRNA:", length(found), "/73; missing:", paste(setdiff(all73, found), collapse = ", "), "\n")

z <- t(scale(t(expr_mat[found, ])))
prog <- colMeans(z, na.rm = TRUE)
names(prog) <- colnames(expr_mat)

# ============================================================
# 3. 细胞注释对齐 + patient-level pseudobulk
# ============================================================
# annot 的 cell_id 与表达矩阵细胞名对齐（可能有格式差异）
annot$cell_id_clean <- gsub("-comb$", "", annot$cell_id)
expr_cells <- colnames(expr_mat)
# 尝试匹配
common <- intersect(annot$cell_id_clean, expr_cells)
cat("Matched cells:", length(common), "\n")

annot_m <- annot[match(common, annot$cell_id_clean), ]
prog_m <- prog[common]
annot_m$prog <- prog_m
annot_m$malignant <- ifelse(annot_m$cell_type == "Malignant", 1, 0)
annot_m$cell_group <- ifelse(annot_m$malignant == 1, "Malignant", "Non-malignant")

# ============================================================
# 4. Patient-level pseudobulk (mean program per patient × cellgroup)
# ============================================================
pb <- aggregate(prog ~ patient + cell_group, data = annot_m, FUN = mean)
cat("\nPatient × group pseudobulk matrix:\n")
print(table(annot_m$patient, annot_m$cell_group))

# 每患者 Malignant vs Non-malignant 配对比较（wilcoxon signed-rank / 符号检验）
wide <- reshape(pb, idvar = "patient", timevar = "cell_group", direction = "wide")
colnames(wide) <- gsub("prog\\.", "", colnames(wide))
# 列名含连字符时用位置访问
mal_col <- grep("Malignant", colnames(wide), value = TRUE)[1]
non_col <- grep("Non", colnames(wide), value = TRUE)[1]
cat("\nwide columns:", colnames(wide), "\n")
wide <- wide[complete.cases(wide), ]
cat("Patients with both groups:", nrow(wide), "\n")
if (nrow(wide) >= 3) {
  mal <- wide[[mal_col]]; non <- wide[[non_col]]
  d <- mal - non
  cat(sprintf("Paired difference: median %.4f, mean %.4f (SD %.4f)\n", median(d), mean(d), sd(d)))
  # Wilcoxon signed-rank
  wr <- wilcox.test(mal, non, paired = TRUE)
  cat(sprintf("Wilcoxon signed-rank P (patient-level paired): %.4g\n", wr$p.value))
  # 正方向比例
  cat(sprintf("Patients with Malignant > Non-malignant: %d / %d\n", sum(d > 0), length(d)))
}

# ============================================================
# 5. 输出
# ============================================================
write.csv(pb, file.path(OUT_DIR, "GSE103322_patient_pseudobulk.csv"), row.names = FALSE)
res <- data.frame(
  Metric = c("Cells_total", "Patients_total", "Genes_73_detected", "PRKN_detected",
             "Patient_pairs_both_groups",
             "Paired_median_diff_Malignant_minus_NonMalig",
             "Wilcoxon_signed_rank_P", "Patients_Malig_gt_NonMalig"),
  Value = c(nrow(annot_m), length(unique(annot_m$patient)), length(found),
            "PRKN" %in% found, nrow(wide),
            if (exists("d")) round(median(d), 4) else NA,
            if (exists("wr")) signif(wr$p.value, 3) else NA,
            if (exists("d")) sum(d > 0) else NA)
)
write.csv(res, file.path(OUT_DIR, "GSE103322_patient_level_summary.csv"), row.names = FALSE)
print(res)

cat("\n=== P0-4 DONE ===")
