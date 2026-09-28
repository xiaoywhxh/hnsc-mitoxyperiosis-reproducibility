## ============================================================
## M4: 线粒体应激代理分析（TCGA 501 + GSE103322 单细胞）
## 高/低 program 组比较: mtDNA 含量 / NRF2-ROS 防御 / UPRmt
## ============================================================
suppressMessages({library(survival)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"

# ---------- 应激代理基因集 ----------
mt_genes <- c('MT-ND1','MT-ND2','MT-ND4','MT-ND5','MT-CYB','MT-CO1','MT-CO2','MT-ATP6','MT-ATP8')
ros_genes <- c('SOD2','CAT','GPX1','GPX4','NQO1','HMOX1','SESN2','GCLC','GCLM','TXN','TXNRD1','PRDX1','PRDX3')
uprmt <- c('HSPD1','HSPE1','CLPP','LONP1','HSPA9','YME1L1','AFG3L2','SPG7','DNAJA3','DDIT3','ATF4','ATF5')
proxy_sets <- list(MitoContent_mtDNA = mt_genes, NRF2_ROS_defense = ros_genes, UPRmt = uprmt)

# 基因集分数（mean-Z）
score_set <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) < 3) return(rep(NA, ncol(mat)))
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

# ---------- TCGA 501 ----------
cat("===== TCGA 501 =====\n")
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]

common67 <- read.csv("D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation/program_common_genes_3OS.csv")$common_gene
expr67 <- log2(tumor_mat[common67, df$sample] + 1)
z67 <- t(scale(t(expr67)))
prog <- colMeans(z67, na.rm = TRUE)
df$prog <- prog[df$sample]
df$group <- ifelse(df$prog > median(df$prog), "High", "Low")

# 表达矩阵用 log2(TPM+1)
expr_all <- log2(tumor_mat[, df$sample] + 1)
res_tcga <- data.frame()
for (nm in names(proxy_sets)) {
  sc <- score_set(expr_all, proxy_sets[[nm]])
  df[[paste0("proxy_", nm)]] <- sc
  hi <- sc[df$group == "High"]; lo <- sc[df$group == "Low"]
  w <- wilcox.test(hi, lo)
  r <- cor(df$prog, sc, method = "spearman")
  res_tcga <- rbind(res_tcga, data.frame(
    Dataset = "TCGA", Proxy = nm, Genes = length(intersect(proxy_sets[[nm]], rownames(expr_all))),
    High_mean = mean(hi, na.rm = TRUE), Low_mean = mean(lo, na.rm = TRUE),
    Diff = mean(hi, na.rm = TRUE) - mean(lo, na.rm = TRUE),
    Wilcox_P = w$p.value, Spearman_rho = r
  ))
  cat(sprintf("  %s: High=%.3f Low=%.3f Wilcoxon P=%.3g rho=%.3f\n",
              nm, mean(hi, na.rm = TRUE), mean(lo, na.rm = TRUE), w$p.value, r))
}
res_tcga$FDR <- p.adjust(res_tcga$Wilcox_P, method = "BH")

# 单基因补充：SOD2 / SESN2 / DDIT3（关键应激基因）
key_genes <- c("SOD2", "SESN2", "DDIT3", "HSPD1", "LONP1")
cat("\n  关键单基因（High vs Low, Wilcoxon）:\n")
for (g in key_genes) {
  if (g %in% rownames(expr_all)) {
    hi <- expr_all[g, df$group == "High"]; lo <- expr_all[g, df$group == "Low"]
    w <- wilcox.test(hi, lo)
    cat(sprintf("    %s: P=%.3g (High %.2f vs Low %.2f)\n", g, w$p.value, mean(hi), mean(lo)))
  }
}

# ---------- GSE103322 单细胞 ----------
cat("\n===== GSE103322 单细胞（Malignant） =====\n")
# 解析表达文件（行 1=表头，行 2-6 元数据，行 7+ 基因；基因名带引号）
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"), open = "rt"), warn = FALSE)
meta_rows <- 6
hdr <- strsplit(lines[1], "\t")[[1]]
samples <- hdr[-1]
gene_lines <- lines[(meta_rows + 1):length(lines)]
genes_sc <- sapply(gene_lines, function(l) gsub("^'|'$|\"", "", strsplit(l, "\t")[[1]][1]))
expr_list <- lapply(gene_lines, function(l) as.numeric(strsplit(l, "\t")[[1]][-1]))
expr_sc <- do.call(rbind, expr_list)
rownames(expr_sc) <- genes_sc
colnames(expr_sc) <- samples
cat("scRNA 表达矩阵:", nrow(expr_sc), "x", ncol(expr_sc), "\n")

# 细胞注释（直接用 cell_id 交集，避免 match 错位）
annot <- readRDS(file.path(DATA_DIR, "GSE103322_cell_annotations.Rds"))
mal_cells <- annot$cell_id[annot$cancer_class == 1]
malignant_idx <- which(colnames(expr_sc) %in% mal_cells)
cat("Malignant 细胞:", length(malignant_idx), "\n")

# Malignant 内 program score（67-gene）与应激代理
g67_sc <- intersect(common67, rownames(expr_sc))
cat("scRNA 67-gene 覆盖:", length(g67_sc), "/", length(common67), "\n")
expr_mal <- expr_sc[, malignant_idx]
prog_mal <- score_set(expr_mal, g67_sc)
grp_mal <- ifelse(prog_mal > median(prog_mal, na.rm = TRUE), "High", "Low")

res_sc <- data.frame()
for (nm in names(proxy_sets)) {
  sc <- score_set(expr_mal, proxy_sets[[nm]])
  if (all(is.na(sc)) || sum(!is.na(sc)) < 10) {
    cat(sprintf("  %s: 覆盖不足/全 NA（跳过）\n", nm))
    next
  }
  hi <- na.omit(sc[grp_mal == "High"]); lo <- na.omit(sc[grp_mal == "Low"])
  w <- tryCatch(wilcox.test(hi, lo), error = function(e) NULL)
  r <- cor(prog_mal, sc, method = "spearman", use = "complete.obs")
  if (is.null(w)) {
    cat(sprintf("  %s: wilcox 失败\n", nm)); next
  }
  res_sc <- rbind(res_sc, data.frame(
    Dataset = "GSE103322_Malignant", Proxy = nm, Genes = length(intersect(proxy_sets[[nm]], rownames(expr_mal))),
    High_mean = mean(hi, na.rm = TRUE), Low_mean = mean(lo, na.rm = TRUE),
    Diff = mean(hi, na.rm = TRUE) - mean(lo, na.rm = TRUE),
    Wilcox_P = w$p.value, Spearman_rho = r
  ))
  cat(sprintf("  %s: High=%.3f Low=%.3f Wilcoxon P=%.3g rho=%.3f\n",
              nm, mean(hi, na.rm = TRUE), mean(lo, na.rm = TRUE), w$p.value, r))
}
res_sc$FDR <- p.adjust(res_sc$Wilcox_P, method = "BH")

# 输出
res_all <- rbind(res_tcga, res_sc)
write.csv(res_all, file.path(OUT_DIR, "M4_stress_proxy_summary.csv"), row.names = FALSE)
cat("\n已保存 M4_stress_proxy_summary.csv\n")
print(res_all, digits = 4, row.names = FALSE)
cat("\nDONE\n")
