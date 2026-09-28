###############################################################################
# Mitoxyperiosis 73-gene signature definition
# Cell 2025.11 first report — 4 functional modules
# Date: 2026-07-20
#
# This file defines all gene sets directly (no RDS dependency).
# Source this file to get gene_sets_mito object.
###############################################################################

# Core module (14 genes): mitochondrial outer membrane permeability regulators
Mitoxyperiosis_Core <- c(
  "PRKN", "VDAC1", "VDAC2", "VDAC3", "BAX", "BAK1", "BID",
  "BBC3", "PMAIP1", "BCL2", "BCL2L1", "MCL1", "TSPO", "AIFM1"
)

# mTOR signaling module (22 genes): growth-survival decision network
# NOTE: RRAGA/B/C/D are Rag GTPases (HGNC official symbols).
#   Old Rds used non-standard "RAGA/B/C/D" which silently failed intersect().
#   Phase 1 incorrectly replaced them with PIK3R1/3/5+STK11 (opposite direction).
#   Restored to RRAGA/B/C/D on 2026-07-25 per root cause analysis.
Mitoxyperiosis_mTOR <- c(
  "MTOR", "RICTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "PRR5L", "DEPTOR",
  "TSC1", "TSC2", "RHEB", "AKT1", "AKT2", "AKT3",
  "PTEN", "PIK3CA", "PIK3CB", "PIK3CD", "RRAGA", "RRAGB", "RRAGC", "RRAGD"
)

# Mitochondrial dynamics module (20 genes): fission-fusion-mitophagy
Mitoxyperiosis_MitoDynamics <- c(
  "DNM1L", "FIS1", "MFF", "MIEF1", "MIEF2", "MFN1", "MFN2", "OPA1",
  "MARCHF5", "PINK1", "BCL2", "BAX", "SLC25A3", "SLC25A5", "SLC25A6",
  "VDAC1", "VDAC2", "VDAC3", "TSPO", "CYCS"
)

# Metabolic-immune crosstalk module (26 genes): metabolism → inflammasome → immune
Mitoxyperiosis_MetaboImmune <- c(
  "HK1", "HK2", "PFKFB3", "PKM", "LDHA", "LDHB", "IDH1", "IDH2",
  "MDH1", "MDH2", "CS", "GLS", "GLS2", "GOT1", "GOT2",
  "NLRP3", "CASP1", "CASP4", "CASP5", "STING1", "CGAS",
  "MAOB", "AIFM1", "ENDOG", "BBC3", "PMAIP1"
)

# All genes (73 unique, after deduplication)
Mitoxyperiosis_All <- unique(c(
  Mitoxyperiosis_Core, Mitoxyperiosis_mTOR,
  Mitoxyperiosis_MitoDynamics, Mitoxyperiosis_MetaboImmune
))

# Package as list (same format as gene_sets_mitoxyperiosis.Rds)
gene_sets_mito <- list(
  Mitoxyperiosis_Core = Mitoxyperiosis_Core,
  Mitoxyperiosis_mTOR = Mitoxyperiosis_mTOR,
  Mitoxyperiosis_MitoDynamics = Mitoxyperiosis_MitoDynamics,
  Mitoxyperiosis_MetaboImmune = Mitoxyperiosis_MetaboImmune,
  Mitoxyperiosis_All = Mitoxyperiosis_All
)

cat("Mitoxyperiosis基因集定义完成!\n")
for (gs in names(gene_sets_mito)) {
  cat("  ", gs, ":", length(gene_sets_mito[[gs]]), "genes\n")
}
