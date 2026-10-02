"""
UCell-style rank-based signature score (Mannen et al. 2022, CSBJ).

Definition implemented (matches UCell::calculate_Uscore with maxRank=1500):
  For each cell/spot j:
    1. Rank ALL genes by expression (descending). Genes with equal expression
       share the average rank; zero-expression genes occupy the lowest ranks.
    2. Take the signature genes that are expressed (>0) and whose descending
       rank <= maxRank. Let n = their count, S = sum of their ranks.
    3. U = 1 - (S - Smin) / (Smax - Smin)
         Smin = n(n+1)/2
         Smax = n(2*maxRank - n + 1)/2
    Bounded to [0, 1]; higher = higher signature expression.
  Computed INDEPENDENTLY per sample (no cross-sample depth normalisation).
"""
import numpy as np


def ucell_score_per_spot(mat_csc, sig_rows, max_rank=1500):
    """
    mat_csc : scipy.sparse CSC matrix (n_genes x n_spots), raw counts
    sig_rows: 1-D array of gene row indices for the signature
    returns : (U float64[n_spots], n_detected int[n_spots])
    """
    n_genes, n_spots = mat_csc.shape
    indptr, indices, data = mat_csc.indptr, mat_csc.indices, mat_csc.data
    sig_set = np.asarray(sig_rows)
    U = np.zeros(n_spots, dtype=np.float64)
    ndet = np.zeros(n_spots, dtype=np.int32)

    for j in range(n_spots):
        st, en = indptr[j], indptr[j + 1]
        if en == st:
            continue
        rows = indices[st:en]            # gene row indices with count>0
        vals = data[st:en]
        m = rows.size
        # descending rank among present genes (ties -> average)
        order = np.argsort(-vals, kind="stable")
        rk = np.empty(m, dtype=np.float64)
        # handle ties by average rank
        sv = vals[order]
        r = np.empty(m, dtype=np.float64)
        i = 0
        while i < m:
            k = i
            while k + 1 < m and sv[k + 1] == sv[i]:
                k += 1
            avg = (i + k) / 2.0 + 1.0    # 1-based average rank
            r[i:k + 1] = avg
            i = k + 1
        rk[order] = r
        # Global descending rank: expressed genes occupy the top m ranks (1..m);
        # zero-expression genes share ranks m+1..n_genes. So a present gene's
        # global descending rank IS its rank among present genes.
        ranks_desc = rk
        # signature genes present
        pos = np.searchsorted(rows, sig_set)          # rows is sorted? CSC indices sorted within column -> yes
        ok = (pos < m) & (rows[np.clip(pos, 0, m - 1)] == sig_set)
        sel = ranks_desc[pos[ok]]
        ndet[j] = int(ok.sum())
        sel = sel[sel <= max_rank]
        n = sel.size
        if n == 0:
            U[j] = 0.0
            continue
        S = sel.sum()
        Smin = n * (n + 1) / 2.0
        Smax = n * (2 * max_rank - n + 1) / 2.0
        U[j] = 1.0 - (S - Smin) / (Smax - Smin)
    return U, ndet


def mean_expression_score(counts_csc_log, sig_rows, gene_mean, gene_sd):
    """
    Robustness-only comparator: gene-wise z-scored mean expression across signature.
    counts_csc_log : log1p(CPM) CSC (n_genes x n_spots)
    gene_mean/sd   : arrays over the SIGNATURE genes (in sig_rows order)
    returns        : mean z across signature genes per spot
    """
    sub = counts_csc_log[sig_rows, :].toarray()       # n_sig x n_spots
    z = (sub - gene_mean[:, None]) / gene_sd[:, None]
    return z.mean(axis=0)
