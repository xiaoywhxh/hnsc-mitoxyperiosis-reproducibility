## P0-D: site×program interaction + direction-signed survival (审稿意见 3 R2-M2/R3-M4/R1-M5)
## 1) program × site interaction P (TCGA 501) — interaction-first
## 2) direction-signed program score 的 survival association (unsigned vs signed HR)
## 3) 单细胞 patient cell counts（供 Supplementary）
## 4) stage/site missingness table

suppressMessages({ library(survival) })

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm], OS_status = clin$OS_status[cm],
                 age = clin$age[cm], gender = clin$gender[cm], stage = clin$ajcc_stage[cm],
                 stringsAsFactors = FALSE)
df <- df[!is.na(df$OS_time), ]
df$patient <- substr(df$sample, 1, 12)

# site（GDC）
site_map <- read.delim(file.path(OUT_DIR, "tcga_hnsc_primary_site.tsv"), check.names = FALSE, quote = "")
site_map$patient <- substr(site_map$submitter_id, 1, 12)
sm <- match(df$patient, site_map$patient)
ps <- tolower(site_map$primary_site[sm])
df$site <- ifelse(grepl("tongue|mouth|gum|palate|lip|floor", ps), "Oral cavity",
           ifelse(grepl("larynx|hypopharyn", ps), "Larynx/Hypopharynx",
           ifelse(grepl("tonsil|oropharyn|base of tongue", ps), "Oropharynx", "Other")))

# program score
core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
all73 <- unique(c(core_genes, mtor_genes, mito_genes, metab_genes))
expr <- log2(tumor_mat[intersect(all73, rownames(tumor_mat)), df$sample] + 1)
z <- t(scale(t(expr)))
prog <- colMeans(z, na.rm = TRUE)
df$prog <- prog[df$sample]

# ============================================================
# 1. site × program interaction (interaction-first)
# ============================================================
cat("\n===== 1. program × site interaction (TCGA 501) =====\n")
# 用 Oral cavity 作为 reference（最大组），LRT 比较含/不含交互项
df$site <- factor(df$site, levels = c("Oral cavity", "Larynx/Hypopharynx", "Oropharynx"))
fit_main <- coxph(Surv(OS_time, OS_status) ~ scale(prog) + site, data = df)
fit_int  <- coxph(Surv(OS_time, OS_status) ~ scale(prog) * site, data = df)
lr <- anova(fit_main, fit_int, test = "Chisq")
cat("LRT interaction P:", signif(lr$`Pr(>|Chi|)`[2], 3), "\n")
s_int <- summary(fit_int)
print(s_int$coefficients)
write.csv(data.frame(test = "program x site interaction (LRT, 3-site, Oral cavity ref)",
                     P = lr$`Pr(>|Chi|)`[2]),
          file.path(OUT_DIR, "program_site_interaction.csv"), row.names = FALSE)

# 用 Larynx 作 reference 再算（敏感性）
df$site2 <- relevel(df$site, ref = "Larynx/Hypopharynx")
fit_int2 <- coxph(Surv(OS_time, OS_status) ~ scale(prog) * site2, data = df)
cat("Alt reference interaction LRT:\n")
print(anova(coxph(Surv(OS_time, OS_status) ~ scale(prog) + site2, data = df),
            fit_int2, test = "Chisq"))

# ============================================================
# 2. direction-signed score survival (unsigned vs signed)
# ============================================================
cat("\n===== 2. Direction-signed program score survival =====\n")
neg_genes <- c("BCL2","BCL2L1","MCL1","PTEN","TSC1","TSC2","DEPTOR")
dir_w <- setNames(rep(1, 73), all73)
dir_w[neg_genes] <- -1
prog_signed <- colMeans(dir_w[all73] * z, na.rm = TRUE)
fit_u <- coxph(Surv(OS_time, OS_status) ~ scale(prog), data = df)
fit_s <- coxph(Surv(OS_time, OS_status) ~ scale(prog_signed), data = df)
cat(sprintf("Unsigned: HR=%.3f (%.3f-%.3f) P=%.4f\n",
            exp(coef(fit_u)), exp(confint(fit_u))[1], exp(confint(fit_u))[2],
            summary(fit_u)$coefficients[5]))
cat(sprintf("Signed:   HR=%.3f (%.3f-%.3f) P=%.4f\n",
            exp(coef(fit_s)), exp(confint(fit_s))[1], exp(confint(fit_s))[2],
            summary(fit_s)$coefficients[5]))
write.csv(data.frame(Score = c("unsigned", "direction-signed"),
                     HR = c(exp(coef(fit_u)), exp(coef(fit_s))),
                     LCL = c(exp(confint(fit_u))[1], exp(confint(fit_s))[1]),
                     UCL = c(exp(confint(fit_u))[2], exp(confint(fit_s))[2]),
                     P = c(summary(fit_u)$coefficients[5], summary(fit_s)$coefficients[5])),
          file.path(OUT_DIR, "program_signed_survival.csv"), row.names = FALSE)

# ============================================================
# 3. stage/site missingness
# ============================================================
cat("\n===== 3. Missingness =====\n")
cat("Total 501:\n")
cat("  stage NA:", sum(is.na(df$stage)), "\n")
cat("  site NA:", sum(is.na(df$site) | df$site == "Other"), "\n")
cat("  complete (age+sex+stage):", sum(!is.na(df$age) & !is.na(df$stage)), "\n")
cat("  complete (+site):", sum(!is.na(df$age) & !is.na(df$stage) & df$site != "Other"), "\n")
print(table(df$stage, useNA = "ifany"))
print(table(df$site, useNA = "ifany"))
miss <- data.frame(Variable = c("age", "sex", "stage", "site(3-class)", "complete age+sex+stage", "complete +site"),
                   N_complete = c(sum(!is.na(df$age)), sum(!is.na(df$gender)),
                                  sum(!is.na(df$stage)), sum(df$site != "Other"),
                                  sum(!is.na(df$age) & !is.na(df$stage)),
                                  sum(!is.na(df$age) & !is.na(df$stage) & df$site != "Other")),
                   N_missing = c(sum(is.na(df$age)), sum(is.na(df$gender)),
                                 sum(is.na(df$stage)), sum(df$site == "Other" | is.na(df$site)),
                                 sum(is.na(df$age) | is.na(df$stage)),
                                 sum(is.na(df$age) | is.na(df$stage) | df$site == "Other")))
write.csv(miss, file.path(OUT_DIR, "missingness_table.csv"), row.names = FALSE)
print(miss)

# ============================================================
# 4. 单细胞 patient cell counts（供 Supplementary）
# ============================================================
cat("\n===== 4. Single-cell patient counts =====\n")
pb <- read.csv("D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/11_singlecell/GSE103322_patient_pseudobulk.csv")
tab <- table(pb$patient, pb$cell_group)
print(tab)
write.csv(as.data.frame.matrix(tab),
          "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/11_singlecell/GSE103322_patient_cell_counts.csv")

cat("\n=== P0-D DONE ===")
