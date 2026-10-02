"""
STEP 1: compute 67-gene UCell score (Primary) and mean-expression score
(robustness only) for each of the 12 samples, independently.
"""
import h5py, numpy as np, os, sys, json
from scipy import sparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ucell import ucell_score_per_spot, mean_expression_score

BASE = 'D:/GSE208253_spatial/GSE208253_processed'
OUT  = 'D:/GSE208253_spatial/spatial_analysis/score'
os.makedirs(OUT, exist_ok=True)

GENES = [g.strip() for g in open('D:/GSE208253_spatial/_67gene_extract/67gene_list.txt', encoding='utf-8') if g.strip()]
SAMPLES = [(f's{i}', f'GSM{6339630+i}') for i in range(1, 13)]

all_rows = []
summary = []
for sN, gsm in SAMPLES:
    h5p = f'{BASE}/{gsm}_{sN}/filtered_feature_bc_matrix.h5'
    with h5py.File(h5p, 'r') as f:
        indptr  = f['matrix/indptr'][()]
        indices = f['matrix/indices'][()]
        data    = f['matrix/data'][()]
        barcodes = [b.decode() for b in f['matrix/barcodes'][()]]
        feat_names = np.array([x.decode() for x in f['matrix/features/name'][()]])
    n_spots = len(barcodes)
    # sort indices within each column (CSC guarantees sorted, but enforce)
    full = sparse.csc_matrix((data, indices, indptr), shape=(len(feat_names), n_spots))
    full.sort_indices()
    sig_rows = np.array([np.where(feat_names == g)[0][0] for g in GENES])

    # --- Primary: UCell ---
    U, ndet = ucell_score_per_spot(full, sig_rows, max_rank=1500)

    # --- Robustness: mean z-scored log1p(CPM) of 67 genes ---
    lib = np.asarray(full.sum(axis=0)).ravel()
    lib[lib == 0] = 1
    # CSC: data is ordered by column; expand per-entry column index
    col_of_entry = np.repeat(np.arange(n_spots), np.diff(full.indptr))
    csc = full.astype(np.float32)
    csc.data = np.log1p(csc.data * (1e6 / lib)[col_of_entry])
    sub = csc[sig_rows, :].toarray()
    gmean = sub.mean(axis=1); gsd = sub.std(axis=1); gsd[gsd == 0] = 1.0
    meanZ = mean_expression_score(csc, sig_rows, gmean, gsd)

    for i, bc in enumerate(barcodes):
        all_rows.append(f'{sN}\t{bc}\t{U[i]:.6f}\t{meanZ[i]:.6f}\t{ndet[i]}\t{int(lib[i])}')
    summary.append(dict(sample=sN, n_spots=n_spots,
                        ucell_mean=float(U.mean()), ucell_sd=float(U.std()),
                        ucell_median=float(np.median(U)), ucell_p10=float(np.percentile(U,10)),
                        ucell_p90=float(np.percentile(U,90)),
                        meanz_mean=float(meanZ.mean()), meanz_sd=float(meanZ.std())))
    print(f'{sN}: n={n_spots} UCell mean={U.mean():.4f} sd={U.std():.4f} median={np.median(U):.4f} '
          f'p10={np.percentile(U,10):.4f} p90={np.percentile(U,90):.4f} | meanZ mean={meanZ.mean():.4f}')

with open(f'{OUT}/spot_scores.tsv', 'w', encoding='utf-8') as f:
    f.write('sample\tbarcode\tucell\tmeanz\tn_gene_detected\tlibrary_size\n')
    f.write('\n'.join(all_rows) + '\n')
json.dump(summary, open(f'{OUT}/_score_summary.json', 'w'), indent=1)
print('\nWritten:', f'{OUT}/spot_scores.tsv')
