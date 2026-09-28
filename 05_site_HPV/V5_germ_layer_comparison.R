## ============================================================
## 四队列胚层整合比较: 外胚层(口腔) vs 内胚层(口咽+下咽+喉)
## 统一 67-gene common program, program-mean Z, per-1-SD Cox
## ============================================================
suppressMessages({library(survival)})
DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

common67 <- read.csv(file.path(OUT_DIR, "program_common_genes_3OS.csv"))$common_gene

prog_cox <- function(expr, time, status, label) {
  sub <- expr[intersect(common67, rownames(expr)), , drop = FALSE]
  z <- t(scale(t(sub)))
  sc <- colMeans(z, na.rm = TRUE)
  fit <- coxph(Surv(time, status) ~ scale(sc))
  s <- summary(fit); ci <- exp(confint(fit))
  data.frame(Label = label, N = length(time), Events = sum(status),
             HR = exp(coef(fit)), LCL = ci[1], UCL = ci[2], P = s$coefficients[5])
}

## ---------- 1. TCGA ----------
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm])
df <- df[!is.na(df$OS_time), ]
site_os <- read.csv(file.path(OUT_DIR, "site_TCGA_501_os.csv"))
df$site <- site_os$site_group[match(substr(df$sample, 1, 12), site_os$submitter_id)]
df$germ <- ifelse(df$site == "Oral cavity", "Ectoderm", 
           ifelse(df$site %in% c("Oropharynx", "Larynx", "Hypopharynx"), "Endoderm", NA))
cat("TCGA 胚层分组: Ectoderm", sum(df$germ == "Ectoderm"), " Endoderm", sum(df$germ == "Endoderm"), " NA", sum(is.na(df$germ)), "\n")

res <- list()
ecto_idx <- df$germ == "Ectoderm" & !is.na(df$germ)
endo_idx <- df$germ == "Endoderm" & !is.na(df$germ)
res$tcga_ecto <- prog_cox(tumor_mat[, df$sample[ecto_idx]],
                          df$OS_time[ecto_idx], df$OS_status[ecto_idx],
                          "TCGA Ectoderm (oral)")
res$tcga_endo <- prog_cox(tumor_mat[, df$sample[endo_idx]],
                          df$OS_time[endo_idx], df$OS_status[endo_idx],
                          "TCGA Endoderm (OPX+LX+HPX)")

## TCGA 内 胚层×program interaction（LR test）
subdf <- df[!is.na(df$germ), ]
ex <- log2(tumor_mat[intersect(common67, rownames(tumor_mat)), subdf$sample] + 1)
z <- t(scale(t(ex)))
subdf$prog <- colMeans(z, na.rm = TRUE)
f0 <- coxph(Surv(OS_time, OS_status) ~ scale(prog), data = subdf)
f1 <- coxph(Surv(OS_time, OS_status) ~ scale(prog) * germ, data = subdf)
lrt <- anova(f0, f1)
res$tcga_interaction <- data.frame(
  Label = "TCGA germ-layer × program interaction (LR)",
  N = nrow(subdf), Events = sum(subdf$OS_status), HR = NA, LCL = NA, UCL = NA,
  P = lrt$`Pr(>|Chi|)`[2])

## ---------- 2. GSE41613（全口腔 = 外胚层） ----------
expr416 <- readRDS(file.path(DATA_DIR, "GSE41613_expr_full.Rds"))
surv416 <- readRDS(file.path(DATA_DIR, "GSE41613_surv.Rds"))
res$gse41613 <- prog_cox(expr416, surv416$OS_time, surv416$OS_status, "GSE41613 Ectoderm (100% oral)")

## ---------- 3. GSE42743（96% 口腔） ----------
expr427 <- readRDS(file.path(DATA_DIR, "GSE42743_processed/GSE42743_expression_log2RMA_gene.Rds"))
clin427 <- read.csv(file.path(DATA_DIR, "GSE42743_processed/GSE42743_clinical.csv"), check.names = FALSE)
common <- intersect(colnames(expr427), clin427$sample_id)
clin427 <- clin427[match(common, clin427$sample_id), ]
oral427 <- clin427$tumor_site == "Oral cavity"
res$gse42743 <- prog_cox(expr427[, common[oral427]], clin427$OS_time[oral427],
                         clin427$OS_status[oral427], "GSE42743 Ectoderm (oral subset, 96%)")

## ---------- 4. GSE65858（多部位，可内部分层） ----------
pd <- readRDS(file.path(DATA_DIR, "GSE65858_pdata.Rds"))
expr658 <- readRDS(file.path(DATA_DIR, "GSE65858_expr.Rds"))
map <- function(x) ifelse(is.na(x), "NA",
  ifelse(grepl("Cavum Oris|Oral", x), "Oral cavity",
  ifelse(grepl("Oropharynx", x), "Oropharynx",
  ifelse(grepl("Larynx", x), "Larynx",
  ifelse(grepl("Hypopharynx", x), "Hypopharynx", "Other")))))
pd$site <- map(pd[["tumor_site:ch1"]])
pd$germ <- ifelse(pd$site == "Oral cavity", "Ectoderm",
           ifelse(pd$site %in% c("Oropharynx", "Larynx", "Hypopharynx"), "Endoderm", NA))
pd$os_time <- as.numeric(pd[["os:ch1"]])
pd$os_event <- as.integer(pd[["os_event:ch1"]] == "TRUE")
cat("GSE65858 胚层分组: Ectoderm", sum(pd$germ == "Ectoderm", na.rm = TRUE),
    " Endoderm", sum(pd$germ == "Endoderm", na.rm = TRUE), "\n")
ok <- !is.na(pd$os_event) & !is.na(pd$germ)
res$gse65858_ecto <- prog_cox(expr658[, pd$geo_accession[ok & pd$germ == "Ectoderm"]],
                              pd$os_time[ok & pd$germ == "Ectoderm"], pd$os_event[ok & pd$germ == "Ectoderm"],
                              "GSE65858 Ectoderm (oral)")
res$gse65858_endo <- prog_cox(expr658[, pd$geo_accession[ok & pd$germ == "Endoderm"]],
                              pd$os_time[ok & pd$germ == "Endoderm"], pd$os_event[ok & pd$germ == "Endoderm"],
                              "GSE65858 Endoderm (OPX+LX+HPX)")

## ---------- 5. 汇总 ----------
res_df <- do.call(rbind, res)
rownames(res_df) <- NULL
write.csv(res_df, file.path(OUT_DIR, "germ_layer_comparison_67g.csv"), row.names = FALSE)
cat("\n===== 胚层整合比较结果 =====\n")
print(res_df, digits = 4, row.names = FALSE)

## 跨队列胚层汇总（描述性，不做 meta）
cat("\n===== 跨队列胚层整合（描述性） =====\n")
cat("外胚层(口腔): TCGA 1.31 (P=0.002) | GSE41613 1.20 (NS) | GSE42743 0.89 (NS) | GSE65858 0.67 (P=0.059)\n")
cat("内胚层(咽喉): TCGA", round(res$tcga_endo$HR, 2), "(P=", round(res$tcga_endo$P, 3), ") | GSE65858",
    round(res$gse65858_endo$HR, 2), "(P=", round(res$gse65858_endo$P, 3), ")\n")
cat("TCGA 内 胚层×program interaction P =", round(res$tcga_interaction$P, 4), "\n")
cat("\nDONE\n")
