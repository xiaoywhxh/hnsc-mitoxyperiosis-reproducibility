## P1-supp: Program-level multivariable Cox + site-stratified (审稿意见 2 R2-M2/M3)
## 1) overall program Cox 调整 age + sex + stage (+ site)
## 2) site-stratified program association（TCGA 内部，exploratory）
## 3) 输出 Fig2 数据（clinical-adjusted program forest）

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
df <- data.frame(sample = colnames(tumor_mat), OS_time = clin$OS_time[cm],
                 OS_status = clin$OS_status[cm], age = clin$age[cm],
                 gender = clin$gender[cm], stage = clin$ajcc_stage[cm],
                 diag = clin$primary_diagnosis[cm], stringsAsFactors = FALSE)
df <- df[!is.na(df$OS_time), ]
cat("501 cohort:", nrow(df), "\n")

# ============================================================
# 0. Program score (mean-Z overall)
# ============================================================
core_genes <- c("PRKN","VDAC1","VDAC2","VDAC3","BAX","BAK1","BID","BBC3","PMAIP1","BCL2","BCL2L1","MCL1","TSPO","AIFM1")
mtor_genes <- c("MTOR","RICTOR","RPTOR","MLST8","MAPKAP1","PRR5","PRR5L","DEPTOR","TSC1","TSC2","RHEB","AKT1","AKT2","AKT3","PTEN","PIK3CA","PIK3CB","PIK3CD","RRAGA","RRAGB","RRAGC","RRAGD")
mito_genes <- c("DNM1L","FIS1","MFF","MIEF1","MIEF2","MFN1","MFN2","OPA1","MARCHF5","PINK1","SLC25A3","SLC25A5","SLC25A6","CYCS")
metab_genes <- c("HK1","HK2","PFKFB3","PKM","LDHA","LDHB","IDH1","IDH2","MDH1","MDH2","CS","GLS","GLS2","GOT1","GOT2","NLRP3","CASP1","CASP4","CASP5","STING1","CGAS","MAOB","ENDOG")
all73 <- unique(c(core_genes, mtor_genes, mito_genes, metab_genes))
expr <- log2(tumor_mat[intersect(all73, rownames(tumor_mat)), df$sample] + 1)
z <- t(scale(t(expr)))
prog_all <- colMeans(z, na.rm = TRUE)
df$prog <- prog_all[df$sample]

# ============================================================
# 1. Site 分类（GDC primary_site + primary_diagnosis）
# ============================================================
site_map <- read.delim(file.path(OUT_DIR, "tcga_hnsc_primary_site.tsv"), check.names = FALSE, quote = "")
# GDC submitter_id 是 case barcode (15位 TCGA-XX-XXXX)
site_map$patient <- substr(site_map$submitter_id, 1, 12)
df$patient <- substr(df$sample, 1, 12)
sm <- match(df$patient, site_map$patient)
ps <- tolower(site_map$primary_site[sm])
df$site <- ifelse(grepl("tongue|mouth|gum|palate|lip|floor", ps), "Oral cavity",
           ifelse(grepl("larynx|hypopharyn", ps), "Larynx/Hypopharynx",
           ifelse(grepl("tonsil|oropharyn|base of tongue", ps), "Oropharynx",
           "Other")))
cat("\nSite distribution (GDC primary_site):\n"); print(table(df$site, useNA = "ifany"))

# ============================================================
# 2. Multivariable program Cox: model A (age+sex+stage) / B (+site)
# ============================================================
df$gender2 <- factor(df$gender)
df$stage2 <- ifelse(grepl("III|IV", df$stage), "III/IV", ifelse(grepl("I|II", df$stage), "I/II", NA))

fit_u <- coxph(Surv(OS_time, OS_status) ~ scale(prog), data = df)
fit_a <- coxph(Surv(OS_time, OS_status) ~ scale(prog) + age + gender2 + stage2, data = df)
df_c <- df[!is.na(df$stage2), ]
fit_b <- coxph(Surv(OS_time, OS_status) ~ scale(prog) + age + gender2 + stage2 + site,
               data = df_c)

report <- function(fit, label, n) {
  s <- summary(fit)
  cat(sprintf("%-35s HR=%.3f (%.3f-%.3f) P=%.4g  N=%d\n", label,
              exp(coef(fit)["scale(prog)"]), exp(confint(fit)["scale(prog)", 1]),
              exp(confint(fit)["scale(prog)", 2]), s$coefficients["scale(prog)", 5], n))
}
cat("\n===== Program (overall) survival association, adjusted =====\n")
report(fit_u, "Univariable", nrow(df))
report(fit_a, "Adjusted age+sex+stage", nrow(df))
report(fit_b, "Adjusted age+sex+stage+site", nrow(df_c))

# ============================================================
# 3. Site-stratified program Cox (exploratory)
# ============================================================
cat("\n===== Site-stratified program HR (univariable, exploratory) =====\n")
site_res <- data.frame()
for (s in sort(unique(df$site))) {
  sub <- df[df$site == s, ]
  if (nrow(sub) < 30 || sum(sub$OS_status) < 10) next
  fit <- coxph(Surv(OS_time, OS_status) ~ scale(prog), data = sub)
  ss <- summary(fit)
  site_res <- rbind(site_res, data.frame(
    Site = s, N = nrow(sub), Events = sum(sub$OS_status),
    HR = exp(coef(fit)), LCL = exp(confint(fit))[1], UCL = exp(confint(fit))[2],
    P = ss$coefficients[5]))
  cat(sprintf("  %-20s N=%d E=%d HR=%.3f (%.3f-%.3f) P=%.4g\n", s, nrow(sub),
              sum(sub$OS_status), exp(coef(fit)), exp(confint(fit))[1],
              exp(confint(fit))[2], ss$coefficients[5]))
}
write.csv(site_res, file.path(OUT_DIR, "program_site_stratified.csv"), row.names = FALSE)

# ============================================================
# 4. 输出调整森林数据（Fig2）
# ============================================================
adj_res <- data.frame(
  Model = c("Univariable", "Adjusted: age + sex + stage", "Adjusted: age + sex + stage + site"),
  HR = c(exp(coef(fit_u)["scale(prog)"]), exp(coef(fit_a)["scale(prog)"]), exp(coef(fit_b)["scale(prog)"])),
  LCL = c(exp(confint(fit_u)["scale(prog)", 1]), exp(confint(fit_a)["scale(prog)", 1]), exp(confint(fit_b)["scale(prog)", 1])),
  UCL = c(exp(confint(fit_u)["scale(prog)", 2]), exp(confint(fit_a)["scale(prog)", 2]), exp(confint(fit_b)["scale(prog)", 2])),
  P = c(summary(fit_u)$coefficients["scale(prog)", 5], summary(fit_a)$coefficients["scale(prog)", 5], summary(fit_b)$coefficients["scale(prog)", 5]),
  N = c(nrow(df), nrow(df), nrow(df_c))
)
write.csv(adj_res, file.path(OUT_DIR, "program_adjusted_forest.csv"), row.names = FALSE)
print(adj_res)

cat("\n=== P1-supp DONE ===")
