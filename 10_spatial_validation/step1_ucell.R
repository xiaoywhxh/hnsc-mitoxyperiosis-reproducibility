###############################################################################
# STEP 1 (corrected v2): 67-gene UCell score, FULL TRANSCRIPTOME background.
# Uses pre-exported full matrices (.npz) via reticulate-free approach:
# exports were written with scipy.sparse.save_npz -> we re-read with Matrix::readMM
# Instead, we re-derive from the full CSC through a small python pre-step that
# writes per-spot present-gene ranks directly. Here we simply recompute in R by
# loading the npz through the 'npz' helper (saved as .mtx instead).
###############################################################################
suppressPackageStartupMessages({ library(Matrix); library(data.table) })
IN  <- "D:/GSE208253_spatial/spatial_analysis/input"
OUT <- "D:/GSE208253_spatial/spatial_analysis/score"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
GENES <- readLines("D:/GSE208253_spatial/_67gene_extract/67gene_list.txt"); GENES <- GENES[nzchar(GENES)]
MAXRANK <- 1500

res_all <- list()
for (k in 1:12) {
  s <- paste0("s", k)
  full <- readMM(file.path(IN, paste0(s, "_full.mtx")))   # genes x spots
  full <- as(full, "dgCMatrix")
  spots <- readLines(file.path(IN, paste0(s, "_spots.tsv")))
  genes <- readLines(file.path(IN, paste0(s, "_genes_full.tsv")))
  rownames(full) <- genes; colnames(full) <- spots
  n_g <- nrow(full); n_c <- ncol(full)
  gi <- match(GENES, genes); stopifnot(!any(is.na(gi)))
  U <- numeric(n_c); det <- integer(n_c)
  # colSums on sparse is fast; do ranking per spot
  for (col in seq_len(n_c)) {
    x <- full[, col]
    nzidx <- which(x@i >= 0)  # not valid; use explicit
    # dense-ish: get nonzero positions
    xi <- x@i + 1L; xv <- x@x
    if (length(xi) == 0) { U[col] <- 0; next }
    ord <- order(-xv, xi)
    rk <- integer(length(xi)); rk[ord] <- seq_along(xi)
    ranks_desc <- (n_g - length(xi)) + rk
    pos <- match(gi, xi)
    keep <- which(!is.na(pos))
    det[col] <- length(keep)
    pr <- ranks_desc[pos[keep]]
    pr <- pr[pr <= MAXRANK]
    if (length(pr) == 0) { U[col] <- 0; next }
    n <- length(pr); S <- sum(pr)
    U[col] <- 1 - (S - n*(n+1)/2) / (n*(2*MAXRANK - n + 1)/2 - n*(n+1)/2)
  }
  res_all[[s]] <- data.table(sample = s, barcode = spots, ucell = U, n_gene_detected = det)
  cat(sprintf("%s: n=%d mean=%.4f sd=%.4f median=%.4f\n", s, n_c, mean(U), sd(U), median(U)))
}
fwrite(rbindlist(res_all), file.path(OUT, "spot_ucell_score.tsv"), sep = "\t")
cat("DONE\n")
