## ============================================================
## 四队列癌种(解剖部位)比例分析 + 亚组分析可行性评估
## 日期: 2026-08-31
## 队列: TCGA / GSE41613 / GSE65858 / GSE42743
## ============================================================
suppressMessages({
  library(survival)
})

DATA_DIR <- "D:/HNSC_mitoxyperiosis_positron/data"
OUT_DIR  <- "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/08_external_validation"

## ---------- 1. TCGA: GDC primary_site 归类 ----------
site_tab <- read.delim(file.path(OUT_DIR, "tcga_hnsc_primary_site.tsv"),
                       header = TRUE, stringsAsFactors = FALSE)
cat("TCGA site 行数:", nrow(site_tab), "\n")

# 解剖学分类（base of tongue 属口咽部 C01.9，非口腔）
site_map <- function(x) {
  tolower_x <- tolower(x)
  ifelse(grepl("floor of mouth|mouth|gum|palate|lip|cheek|buccal|retromolar|alveolar|other and unspecified parts of tongue", tolower_x),
         "Oral cavity",
  ifelse(grepl("tonsil|base of tongue|oropharynx|lingual tonsil", tolower_x),
         "Oropharynx",
  ifelse(grepl("larynx", tolower_x),
         "Larynx",
  ifelse(grepl("hypopharynx", tolower_x),
         "Hypopharynx",
  ifelse(grepl("nasal|sinus|nasopharynx", tolower_x),
         "Nasal/Paranasal",
         "Other")))))
}
site_tab$site_group <- site_map(site_tab$primary_site)
tcga_site <- table(site_tab$site_group)
cat("\n=== TCGA 部位分组（GDC primary_site, n=", nrow(site_tab), "）===\n", sep="")
print(sort(tcga_site, decreasing = TRUE))

# 详细: primary_site -> group 映射检查
cat("\n=== TCGA primary_site 明细 ===\n")
print(table(site_tab$primary_site, site_tab$site_group))

# 501 例有 OS 的样本部位
tumor_mat <- readRDS(file.path(DATA_DIR, "TCGA_full", "TCGA_HNSC_tumor_full_expr.Rds"))
clin <- readRDS(file.path(DATA_DIR, "TCGA_HNSC_clinical.Rds"))
clin$OS_time <- ifelse(!is.na(clin$days_to_death), clin$days_to_death,
                       ifelse(!is.na(clin$days_to_last_followup), clin$days_to_last_followup, NA))
clin$OS_status <- ifelse(clin$vital_status == "Dead", 1, 0)
pid <- substr(colnames(tumor_mat), 1, 12)
cm <- match(pid, clin$bcr_patient_barcode)
os_ok <- !is.na(clin$OS_time[cm])
os_pat <- pid[os_ok]

site_os <- site_tab[site_tab$submitter_id %in% os_pat, ]
cat("\n=== TCGA 有 OS 样本部位分布（n=", nrow(site_os), "）===\n", sep="")
t_site_os <- sort(table(site_os$site_group), decreasing=TRUE)
print(t_site_os)
cat("占比:\n")
print(round(100 * t_site_os / sum(t_site_os), 1))

# 事件数
os_dead <- clin$OS_status[match(site_os$submitter_id, clin$bcr_patient_barcode)]
t_dead <- table(site_os$site_group, os_dead)
cat("\n=== TCGA 各部位事件数 ===\n")
print(t_dead)

# 保存 TCGA site
write.csv(data.frame(submitter_id = site_os$submitter_id,
                     primary_site = site_os$primary_site,
                     site_group = site_os$site_group,
                     OS_dead = os_dead),
          file.path(OUT_DIR, "site_TCGA_501_os.csv"), row.names = FALSE)

## ---------- 2. GSE65858: tumor_site ----------
pd <- readRDS(file.path(DATA_DIR, "GSE65858_pdata.Rds"))
ts <- pd[["tumor_site:ch1"]]
cat("\n=== GSE65858 部位分布（n=", length(ts), "）===\n", sep="")
print(sort(table(ts, useNA="ifany"), decreasing=TRUE))

# 部位分类统一
gse65858_map <- function(x) {
  ifelse(is.na(x), "NA",
  ifelse(grepl("Cavum Oris|Oral", x), "Oral cavity",
  ifelse(grepl("Oropharynx", x), "Oropharynx",
  ifelse(grepl("Larynx", x), "Larynx",
  ifelse(grepl("Hypopharynx", x), "Hypopharynx", "Other")))))
}
pd$site_group <- gse65858_map(ts)
t658 <- table(pd$site_group)
print(t658)
cat("占比:\n"); print(round(100 * t658 / sum(t658), 1))

# 事件数(OS)
os658 <- pd[["os_event:ch1"]]
t658_ev <- table(pd$site_group, os658)
print(t658_ev)
write.csv(data.frame(geo = pd$geo_accession, site_group = pd$site_group, os_event = os658),
          file.path(OUT_DIR, "site_GSE65858.csv"), row.names = FALSE)

## ---------- 3. GSE42743: tumor_site ----------
clin427 <- read.csv(file.path(DATA_DIR, "GSE42743_processed/GSE42743_clinical.csv"),
                    check.names = FALSE)
cat("\n=== GSE42743 部位分布（n=", nrow(clin427), "）===\n", sep="")
print(sort(table(clin427$tumor_site), decreasing=TRUE))
write.csv(data.frame(sample = clin427$sample_id, site_group = clin427$tumor_site,
                     OS_status = clin427$OS_status),
          file.path(OUT_DIR, "site_GSE42743.csv"), row.names = FALSE)

## ---------- 4. GSE41613: 待 SOFT 数据（占位） ----------
cat("\n=== GSE41613: 待 series matrix 完整下载后补充 ===\n")

## ---------- 5. 汇总表 ----------
cat("\n\n========== 汇总：四队列部位比例 ==========\n")
cat("TCGA(有OS 501): Oral cavity=", t_site_os["Oral cavity"], ", Oropharynx=", t_site_os["Oropharynx"],
    ", Larynx=", t_site_os["Larynx"], ", Hypopharynx=", t_site_os["Hypopharynx"],
    ", Other=", sum(t_site_os[c("Nasal/Paranasal","Other")], na.rm=TRUE), "\n", sep="")
cat("GSE65858(270): Oral=", t658["Oral cavity"], ", Oropharynx=", t658["Oropharynx"],
    ", Larynx=", t658["Larynx"], ", Hypopharynx=", t658["Hypopharynx"], ", NA=", t658["NA"], "\n", sep="")
cat("GSE42743(74): Oral=", table(clin427$tumor_site)["Oral cavity"],
    ", Oropharynx=", table(clin427$tumor_site)["Oropharynx"], "\n", sep="")

cat("\nDONE\n")
