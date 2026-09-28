## B3: Tier A/B/C gene evidence table
## 依据审稿意见 1.1/1.2：将 73 基因重建为三级证据体系
## Tier A = direct mechanistic (原始 Cell 论文实验确立的机制基因)
## Tier B = mechanistically proximal (线粒体通透性/线粒体-质膜互作/氧化膜损伤/mTORC2/线粒体应激)
## Tier C = extended context (糖酵解/TCA/cGAS-STING/inflammasome/一般凋亡/mitophagy/代谢免疫)

# Canonical 73 unique genes per manuscript (2026-08-30)
# Note: 数据一致性发现 — Rds 基因集与论文描述的模块构成不一致，
#       本表以论文 Methods 2.2 的 canonical 73 genes 为准

core_genes <- c("PRKN", "VDAC1", "VDAC2", "VDAC3", "BAX", "BAK1", "BID",
                "BBC3", "PMAIP1", "BCL2", "BCL2L1", "MCL1", "TSPO", "AIFM1")
mtor_genes <- c("MTOR", "RICTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "PRR5L",
                "DEPTOR", "TSC1", "TSC2", "RHEB", "AKT1", "AKT2", "AKT3",
                "PTEN", "PIK3CA", "PIK3CB", "PIK3CD", "RRAGA", "RRAGB", "RRAGC", "RRAGD")
mito_genes <- c("DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "MFN1", "MFN2", "OPA1",
                "MARCHF5", "PINK1", "SLC25A3", "SLC25A5", "SLC25A6", "CYCS")
metab_genes <- c("HK1", "HK2", "PFKFB3", "PKM", "LDHA", "LDHB", "IDH1", "IDH2",
                 "MDH1", "MDH2", "CS", "GLS", "GLS2", "GOT1", "GOT2",
                 "NLRP3", "CASP1", "CASP4", "CASP5", "STING1", "CGAS", "MAOB", "ENDOG")

all73 <- c(core_genes, mtor_genes, mito_genes, metab_genes)
stopifnot(length(all73) == 73, length(unique(all73)) == 73)

module_of <- function(g) {
  if (g %in% core_genes) return("Core")
  if (g %in% mtor_genes) return("mTOR")
  if (g %in% mito_genes) return("MitoDynamics")
  if (g %in% metab_genes) return("MetaboImmune")
  return(NA)
}

# ============================================================
# Tier assignment + evidence annotations
# ============================================================
# 依据: Wang Y et al. Cell 2025 (PMID 41317732) 原始机制:
# 免疫代谢应激 → BAX/BAK1/BID 依赖氧化应激 → 持续线粒体-质膜接触
# → mTORC2 抑制肌动蛋白骨架重塑 → 线粒体滞留外周 → 局部膜氧化损伤 → 裂解

tier_A <- c("BAX", "BAK1", "BID", "MTOR", "RICTOR", "MAPKAP1")  # 直接机制基因(实验确立)
# BAX/BAK1/BID: 原始研究依赖的氧化应激介导者; MTOR/RICTOR/MAPKAP1: mTORC2 催化/特异组分

tier_B <- c("VDAC1", "VDAC2", "VDAC3", "TSPO", "DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2",
            "MFN1", "MFN2", "OPA1", "MARCHF5", "PINK1", "SLC25A3", "SLC25A5", "SLC25A6",
            "CYCS", "AIFM1", "ENDOG", "PRKN", "MLST8", "PRR5", "PRR5L", "DEPTOR",
            "TSC1", "TSC2", "RHEB", "AKT1", "AKT2", "AKT3", "PTEN", "PIK3CA", "PIK3CB", "PIK3CD")
# 线粒体通透性(VDAC/TSPO)、动力学(DNM1L/MFN/OPA1...)、转运(SLC25)、线粒体释放死亡因子(CYCS/AIFM1/ENDOG)、
# mitophagy(PRKN/PINK1)、mTORC2 其余组分(MLST8/PRR5/PRR5L/DEPTOR)、mTOR 上游(TSC/RHEB/AKT/PI3K)

tier_C <- c("BBC3", "PMAIP1", "BCL2", "BCL2L1", "MCL1", "RPTOR",
            "RRAGA", "RRAGB", "RRAGC", "RRAGD",
            "HK1", "HK2", "PFKFB3", "PKM", "LDHA", "LDHB", "IDH1", "IDH2",
            "MDH1", "MDH2", "CS", "GLS", "GLS2", "GOT1", "GOT2",
            "NLRP3", "CASP1", "CASP4", "CASP5", "STING1", "CGAS", "MAOB")
# BCL2 家族其他(BBC3/PMAIP1/BCL2/BCL2L1/MCL1)、mTORC1 组分(RPTOR/raptor)、Rag GTPase(RRAGA-D)、
# 糖酵解/TCA、炎症小体、cGAS-STING、MAOB

tier_of <- function(g) {
  if (g %in% tier_A) return("A")
  if (g %in% tier_B) return("B")
  if (g %in% tier_C) return("C")
  return(NA)
}

# Verify completeness
stopifnot(all(all73 %in% c(tier_A, tier_B, tier_C)))
cat("Tier counts: A=", length(tier_A), " B=", length(tier_B), " C=", length(tier_C), "\n")

# ============================================================
# Evidence annotations per gene (role + evidence source + direction)
# ============================================================
# Manual curated annotations (concise)
annot <- list(
  # Tier A
  BAX = "MOMP 执行者；原始研究氧化应激介导核心", BAK1 = "MOMP 执行者；原始研究氧化应激介导核心",
  BID = "BH3-only 激活子；线粒体膜损伤启动", MTOR = "mTORC1/2 催化亚基；mTORC2 轴核心",
  RICTOR = "mTORC2 特异组分；抑制肌动蛋白重塑的直接机制", MAPKAP1 = "mTORC2 特异组分(SIN1)",
  # Tier B
  VDAC1 = "线粒体外膜通透性通道；mPTP 组分", VDAC2 = "线粒体外膜通透性通道", VDAC3 = "线粒体外膜通透性通道",
  TSPO = "外膜转运蛋白；线粒体应激响应", DNM1L = "线粒体分裂(Drp1)", FIS1 = "线粒体分裂受体",
  MFF = "线粒体分裂因子", MIEF1 = "线粒体分裂/融合调控(MiD1)", MIEF2 = "线粒体分裂/融合调控(MiD2)",
  MFN1 = "线粒体融合", MFN2 = "线粒体融合；ER-线粒体接触", OPA1 = "线粒体嵴重塑/融合",
  MARCHF5 = "线粒体 E3 泛素连接酶(MITOL)", PINK1 = "mitophagy 传感器", PRKN = "mitophagy 执行者(Parkin)",
  SLC25A3 = "线粒体磷酸转运", SLC25A5 = "ANT2 腺苷酸转运；mPTP 候选组分",
  SLC25A6 = "ANT3 腺苷酸转运", CYCS = "细胞色素 c；线粒体死亡信号释放",
  AIFM1 = "凋亡诱导因子；线粒体释放", ENDOG = "线粒体内切酶 G；caspase 非依赖死亡",
  MLST8 = "mTORC1/2 共同组分(GβL)", PRR5 = "mTORC2 组分(Protor-1)", PRR5L = "mTORC2 组分(Protor-2)",
  DEPTOR = "mTORC2 抑制性组分", TSC1 = "mTOR 上游抑制子(hamartin)", TSC2 = "mTOR 上游抑制子(tuberin)",
  RHEB = "mTORC1 激活子；与 mTORC2 通路交叉", AKT1 = "PI3K-AKT-mTOR 轴", AKT2 = "PI3K-AKT-mTOR 轴",
  AKT3 = "PI3K-AKT-mTOR 轴", PTEN = "PI3K 通路抑制子", PIK3CA = "PI3K 催化亚基",
  PIK3CB = "PI3K 催化亚基", PIK3CD = "PI3K 催化亚基",
  # Tier C
  BBC3 = "BH3-only 凋亡启动子(PUMA)", PMAIP1 = "BH3-only 凋亡启动子(NOXA)", BCL2 = "抗凋亡 BCL2 家族",
  BCL2L1 = "抗凋亡 BCL2 家族(Bcl-xL)", MCL1 = "抗凋亡 BCL2 家族", RPTOR = "mTORC1 特异组分(raptor)；非 mTORC2",
  RRAGA = "Rag GTPase；mTORC1 氨基酸感知", RRAGB = "Rag GTPase；mTORC1 氨基酸感知",
  RRAGC = "Rag GTPase；mTORC1 氨基酸感知", RRAGD = "Rag GTPase；mTORC1 氨基酸感知",
  HK1 = "糖酵解第一步；线粒体外膜己糖激酶", HK2 = "糖酵解；线粒体结合己糖激酶",
  PFKFB3 = "糖酵解调控", PKM = "糖酵解丙酮酸激酶", LDHA = "乳酸脱氢酶", LDHB = "乳酸脱氢酶",
  IDH1 = "TCA/IDH", IDH2 = "TCA/IDH", MDH1 = "TCA 苹果酸脱氢酶", MDH2 = "TCA 苹果酸脱氢酶",
  CS = "TCA 柠檬酸合酶", GLS = "谷氨酰胺分解", GLS2 = "谷氨酰胺分解",
  GOT1 = "谷氨酰胺代谢/转氨", GOT2 = "谷氨酰胺代谢/转氨",
  NLRP3 = "炎症小体感受器", CASP1 = "炎症小体效应酶", CASP4 = "非经典炎症小体",
  CASP5 = "非经典炎症小体", STING1 = "cGAS-STING 先天免疫", CGAS = "胞质 DNA 感受器",
  MAOB = "单胺氧化酶；ROS 生成与代谢-免疫交叉"
)

evid <- list(
  # Evidence type: D=direct perturbation in original study, P=pathway member/literature, R=review-based
  BAX = "D", BAK1 = "D", BID = "D", MTOR = "D", RICTOR = "D", MAPKAP1 = "D",
  VDAC1 = "P", VDAC2 = "P", VDAC3 = "P", TSPO = "P", DNM1L = "P", FIS1 = "P",
  MFF = "P", MIEF1 = "P", MIEF2 = "P", MFN1 = "P", MFN2 = "P", OPA1 = "P",
  MARCHF5 = "P", PINK1 = "P", PRKN = "P", SLC25A3 = "P", SLC25A5 = "P", SLC25A6 = "P",
  CYCS = "P", AIFM1 = "P", ENDOG = "P", MLST8 = "P", PRR5 = "P", PRR5L = "P",
  DEPTOR = "P", TSC1 = "P", TSC2 = "P", RHEB = "P", AKT1 = "P", AKT2 = "P",
  AKT3 = "P", PTEN = "P", PIK3CA = "P", PIK3CB = "P", PIK3CD = "P",
  BBC3 = "P", PMAIP1 = "P", BCL2 = "P", BCL2L1 = "P", MCL1 = "P", RPTOR = "P",
  RRAGA = "P", RRAGB = "P", RRAGC = "P", RRAGD = "P",
  HK1 = "R", HK2 = "R", PFKFB3 = "R", PKM = "R", LDHA = "R", LDHB = "R",
  IDH1 = "R", IDH2 = "R", MDH1 = "R", MDH2 = "R", CS = "R", GLS = "R", GLS2 = "R",
  GOT1 = "R", GOT2 = "R", NLRP3 = "R", CASP1 = "R", CASP4 = "R", CASP5 = "R",
  STING1 = "R", CGAS = "R", MAOB = "R"
)

# Direction expectation in mitoxyperiosis (pro-death/positive = promotes mitoxyperiosis biology)
direction <- list(
  BAX = "+", BAK1 = "+", BID = "+", MTOR = "+", RICTOR = "+", MAPKAP1 = "+",
  VDAC1 = "+", VDAC2 = "+", VDAC3 = "+", TSPO = "+", DNM1L = "+", FIS1 = "+",
  MFF = "+", MIEF1 = "+", MIEF2 = "+", MFN1 = "±", MFN2 = "±", OPA1 = "±",
  MARCHF5 = "-", PINK1 = "-", PRKN = "-", SLC25A3 = "+", SLC25A5 = "+", SLC25A6 = "+",
  CYCS = "+", AIFM1 = "+", ENDOG = "+", MLST8 = "+", PRR5 = "+", PRR5L = "+",
  DEPTOR = "-", TSC1 = "-", TSC2 = "-", RHEB = "+", AKT1 = "+", AKT2 = "+",
  AKT3 = "+", PTEN = "-", PIK3CA = "+", PIK3CB = "+", PIK3CD = "+",
  BBC3 = "+", PMAIP1 = "+", BCL2 = "-", BCL2L1 = "-", MCL1 = "-", RPTOR = "±",
  RRAGA = "±", RRAGB = "±", RRAGC = "±", RRAGD = "±",
  HK1 = "+", HK2 = "+", PFKFB3 = "+", PKM = "+", LDHA = "+", LDHB = "±",
  IDH1 = "±", IDH2 = "±", MDH1 = "±", MDH2 = "±", CS = "±", GLS = "+", GLS2 = "±",
  GOT1 = "±", GOT2 = "±", NLRP3 = "+", CASP1 = "+", CASP4 = "+", CASP5 = "+",
  STING1 = "+", CGAS = "+", MAOB = "+"
)

# ============================================================
# Build master table
# ============================================================
df <- data.frame(
  Gene = all73,
  Module = sapply(all73, module_of),
  Tier = sapply(all73, tier_of),
  Annotation_CN = sapply(all73, function(g) annot[[g]]),
  Evidence_Type = sapply(all73, function(g) evid[[g]]),
  Expected_Direction = sapply(all73, function(g) direction[[g]]),
  stringsAsFactors = FALSE
)
df$Evidence_Label <- ifelse(df$Evidence_Type == "D",
  "Direct (original Cell study perturbation)",
  ifelse(df$Evidence_Type == "P", "Pathway member / literature",
         "Extended context (review-based)"))

# Sort: Tier then module
df <- df[order(df$Tier, df$Module, df$Gene), ]
rownames(df) <- NULL

write.csv(df, "D:/HNSC_mitoxyperiosis_positron/V2_Reanalysis/02_gene_program/mitoxyperiosis_gene_evidence_master.csv", row.names = FALSE)

# Summary
cat("\n=== Tier x Module summary ===\n")
print(table(df$Tier, df$Module))

cat("\n=== Tier A genes ===\n")
print(df[df$Tier == "A", c("Gene", "Module", "Annotation_CN")])

cat("\n=== DONE ===")
