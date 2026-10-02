# Moran's I for spatial score using Visium hexagonal adjacency.
# Implements the standard Moran's I with row-standardised weights (W),
# analytically-derived expectation/variance under randomisation (normality),
# then cross-checks on samples without isolated spots via spdep::moran.test.
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(spdep) })
IN  <- "D:/GSE208253_spatial/spatial_analysis/input"
SC  <- "D:/GSE208253_spatial/spatial_analysis/score"
OUT <- "D:/GSE208253_spatial/spatial_analysis/moran"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

spots  <- fread(file.path(IN, "spot_table_all.tsv"), sep = "\t", header = TRUE)
scores <- fread(file.path(SC, "spot_scores.tsv"), sep = "\t", header = TRUE)
spots  <- merge(spots, scores, by = c("sample", "barcode"), all.x = TRUE)

build_nb <- function(dt) {
  n <- nrow(dt)
  lookup <- new.env(hash = TRUE)
  for (i in seq_len(n)) assign(paste(dt$array_row[i], dt$array_col[i]), i, envir = lookup)
  offs <- rbind(c(0,-2), c(0,2), c(-1,-1), c(-1,1), c(1,-1), c(1,1))
  ii <- integer(0); jj <- integer(0)
  for (i in seq_len(n)) {
    r <- dt$array_row[i]; cc <- dt$array_col[i]
    for (k in seq_len(nrow(offs))) {
      v <- mget(paste(r + offs[k,1], cc + offs[k,2]), envir = lookup, ifnotfound = NA)[[1]]
      if (!is.na(v)) { ii <- c(ii, i); jj <- c(jj, v) }
    }
  }
  sparseMatrix(i = ii, j = jj, dims = c(n, n))
}

morans_I_manual <- function(x, W) {
  # W: symmetric binary adjacency (n x n). Row-standardise.
  n <- length(x)
  rs <- Matrix::rowSums(W)
  rs[rs == 0] <- 1                      # isolated spots: zero row
  Ws <- Diagonal(x = 1/rs) %*% W        # row-standardised
  xz <- x - mean(x)
  S0 <- sum(Ws)                         # = number of non-isolated rows
  num <- as.numeric(t(xz) %*% (Ws %*% xz))
  den <- as.numeric(sum(xz^2))
  I <- (n / S0) * (num / den)
  # expectation and variance under normality assumption (row-standardised W)
  S1 <- 0.5 * sum((Ws + Matrix::t(Ws))^2)
  S2 <- sum((Matrix::rowSums(Ws) + Matrix::colSums(Ws))^2)
  k <- (sum(xz^4) / n) / (sum(xz^2) / n)^2
  EI <- -1 / (n - 1)
  A <- n * ((n^2 - 3*n + 3) * S1 - n * S2 + 3 * S0^2)
  B <- k * ((n^2 - n) * S1 - 2*n*S2 + 6 * S0^2)
  C <- (n - 1) * (n - 2) * (n - 3) * S0^2
  VI <- (A - B) / C - EI^2
  z <- (I - EI) / sqrt(VI)
  p <- 2 * pnorm(-abs(z))
  list(I = I, EI = EI, VI = VI, z = z, p = p, kbar = mean(Matrix::rowSums(W)))
}

res <- list()
for (s in paste0("s", 1:12)) {
  dt <- spots[sample == s][order(array_row, array_col)]
  W  <- build_nb(as.data.frame(dt))
  m  <- morans_I_manual(dt$ucell, W)
  # cross-check with spdep if no isolated spots
  chk <- NA_real_
  if (all(Matrix::rowSums(W) > 0)) {
    nb <- lapply(seq_len(nrow(dt)), function(i) which(W[i, ] != 0))
    attr(nb, "class") <- "nb"; attr(nb, "region.id") <- as.character(seq_len(nrow(dt)))
    lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
    chk <- unname(moran.test(dt$ucell, lw, zero.policy = TRUE,
                             alternative = "two.sided")$estimate[["Moran I statistic"]])
  }
  res[[s]] <- data.table(sample = s, n_spots = nrow(dt),
                         mean_n_neighbours = m$kbar,
                         morans_I = m$I, expected_I = m$EI, variance = m$VI,
                         sd = sqrt(m$VI), z = m$z, p_value = m$p,
                         spdep_check_I = chk)
  cat(sprintf("%s: n=%d kbar=%.2f  I=%.4f  E[I]=%.4f  sd=%.4f  z=%.2f  p=%.3g  spdep=%.4f\n",
              s, nrow(dt), m$kbar, m$I, m$EI, sqrt(m$VI), m$z, m$p,
              ifelse(is.na(chk), NA, chk)))
}

res <- rbindlist(res)
res[, fdr_bh := p.adjust(p_value, method = "BH")]
fwrite(res, file.path(OUT, "morans_I_by_sample.tsv"), sep = "\t")

pos <- sum(res$morans_I > 0)
cat(sprintf("\n=== Moran's I direction: %d / 12 samples > 0 ; %d / 12 FDR<0.05 ===\n",
            pos, sum(res$fdr_bh < 0.05)))
cat(sprintf("I range: %.4f .. %.4f ; median %.4f\n",
            min(res$morans_I), max(res$morans_I), median(res$morans_I)))
print(res[, .(sample, n_spots, morans_I = round(morans_I,4), z = round(z,2),
              p_value = signif(p_value,3), fdr_bh = signif(fdr_bh,3),
              spdep = round(spdep_check_I,4))])
cat("\nWritten:", file.path(OUT, "morans_I_by_sample.tsv"), "\n")
