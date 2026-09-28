## M4 修订: GSE103322 应激代理 patient-level 分析（审稿 5-0 Major 2 / 5-1 问题1）
## 每患者 malignant cells 聚合 program/NRF2/UPRmt 均值 → patient 为 n 做 Spearman；cell-level 仅可视化
suppressMessages({library(ggplot2)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/12_ici"

ros_genes <- c('SOD2','CAT','GPX1','GPX4','NQO1','HMOX1','SESN2','GCLC','GCLM','TXN','TXNRD1','PRDX1','PRDX3')
uprmt <- c('HSPD1','HSPE1','CLPP','LONP1','HSPA9','YME1L1','AFG3L2','SPG7','DNAJA3','DDIT3','ATF4','ATF5')

score_set <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  z <- t(scale(t(mat[genes, , drop = FALSE])))
  colMeans(z, na.rm = TRUE)
}

cat(">>> 读 GSE103322...\n")
lines <- readLines(gzfile(file.path(DATA_DIR, "GSE103322_HNSCC_all_data.txt.gz"), open = "rt"), warn = FALSE)
cell_barcodes <- strsplit(lines[1], "\t")[[1]][-1]
gene_start <- 6
n_cells <- length(cell_barcodes)
n_genes <- length(lines) - gene_start
gene_names <- character(n_genes)
expr_mat <- matrix(NA, nrow = n_genes, ncol = n_cells)
for (i in 1:n_genes) {
  parts <- strsplit(lines[gene_start + i], "\t")[[1]]
  gene_names[i] <- gsub("'", "", parts[1])
  expr_mat[i, ] <- as.numeric(parts[2:(n_cells + 1)])
}
rownames(expr_mat) <- gene_names
colnames(expr_mat) <- cell_barcodes

annot <- readRDS(file.path(DATA_DIR, "GSE103322_cell_annotations.Rds"))
mal <- annot$cell_id[annot$cell_type == "Malignant"]
expr_mal <- expr_mat[, mal]
mal_annot <- annot[annot$cell_type == "Malignant", ]
cat("Malignant 细胞:", ncol(expr_mal), " patient:", length(unique(mal_annot$patient)), "\n")

common67 <- read.csv("D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation/program_common_genes_3OS.csv")$common_gene

# per-cell 分数
prog_cell <- score_set(expr_mal, common67)
ros_cell  <- score_set(expr_mal, ros_genes)
uprmt_cell <- score_set(expr_mal, uprmt)

# patient 聚合（malignant cells 均值）
patients <- unique(mal_annot$patient)
df <- data.frame(patient = patients, prog = NA, NRF2 = NA, UPRmt = NA, n_cells = NA)
for (i in seq_along(patients)) {
  idx <- mal_annot$patient == patients[i]
  df$prog[i] <- mean(prog_cell[idx], na.rm = TRUE)
  df$NRF2[i] <- mean(ros_cell[idx], na.rm = TRUE)
  df$UPRmt[i] <- mean(uprmt_cell[idx], na.rm = TRUE)
  df$n_cells[i] <- sum(idx)
}
cat("\npatient-level 聚合（n =", nrow(df), "患者）:\n")
print(df[order(-df$prog), ], row.names = FALSE)

# patient-level Spearman
sp <- function(x, y) {
  ok <- complete.cases(x, y)
  r <- cor(x[ok], y[ok], method = "spearman")
  n <- sum(ok); t <- r * sqrt((n - 2) / (1 - r^2)); p <- 2 * (1 - pnorm(abs(t)))
  c(rho = r, P = p, n = n)
}
r1 <- sp(df$prog, df$NRF2); r2 <- sp(df$prog, df$UPRmt)
cat("\n========================================\n")
cat("patient-level（n =", nrow(df), "）:\n")
cat(sprintf("  program × NRF2/ROS:  rho=%.3f  P=%.3g\n", r1[1], r1[2]))
cat(sprintf("  program × UPRmt:     rho=%.3f  P=%.3g\n", r2[1], r2[2]))
cat("========================================\n")

# 图：patient-level
p1 <- ggplot(df, aes(prog, NRF2)) + geom_point(size = 2.5, color = "#3C5488") +
  geom_smooth(method = "lm", se = TRUE, color = "#E64B35") +
  labs(x = "Program score (patient-level malignant mean)", y = "NRF2/ROS-defense score",
       title = sprintf("patient-level (n=%d): rho=%.3f, P=%.3g", nrow(df), r1[1], r1[2])) +
  theme_bw(base_size = 11)
p2 <- ggplot(df, aes(prog, UPRmt)) + geom_point(size = 2.5, color = "#3C5488") +
  geom_smooth(method = "lm", se = TRUE, color = "#E64B35") +
  labs(x = "Program score (patient-level malignant mean)", y = "UPRmt score",
       title = sprintf("patient-level (n=%d): rho=%.3f, P=%.3g", nrow(df), r2[1], r2[2])) +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_DIR, "M4_patient_level_NRF2.png"), p1, width = 5, height = 4, dpi = 300)
ggsave(file.path(OUT_DIR, "M4_patient_level_UPRmt.png"), p2, width = 5, height = 4, dpi = 300)

write.csv(df, file.path(OUT_DIR, "M4_patient_level_stress.csv"), row.names = FALSE)
cat("已保存 M4_patient_level_stress.csv + 图\n")
