# Rebuild with FULL transcriptome background (36,601 genes) for correct UCell ranking.
# UCell ranks each gene within the WHOLE transcriptome per cell, then takes the
# signature genes' ranks. Using only 67 genes as background saturates the score.
import h5py, numpy as np, os, json
from scipy import sparse
from scipy.io import mmwrite

BASE = 'D:/GSE208253_spatial/GSE208253_processed'
OUT  = 'D:/GSE208253_spatial/spatial_analysis/input'
GENES = [g.strip() for g in open('D:/GSE208253_spatial/_67gene_extract/67gene_list.txt', encoding='utf-8') if g.strip()]
SAMPLES = [(f's{i}', f'GSM{6339630+i}') for i in range(1, 13)]

for sN, gsm in SAMPLES:
    h5p = f'{BASE}/{gsm}_{sN}/filtered_feature_bc_matrix.h5'
    with h5py.File(h5p, 'r') as f:
        indptr  = f['matrix/indptr'][()]
        indices = f['matrix/indices'][()]
        data    = f['matrix/data'][()]
        barcodes = [b.decode() for b in f['matrix/barcodes'][()]]
        feat_names = np.array([x.decode() for x in f['matrix/features/name'][()]])
    n_spots = len(barcodes)
    # full matrix (CSC) as genes x spots CSR for fast column access -> keep CSC, use it
    full = sparse.csc_matrix((data, indices, indptr), shape=(len(feat_names), n_spots))
    # save full matrix as .npz (compressed sparse) to rank against all genes later
    sparse.save_npz(f'{OUT}/{sN}_full.npz', full)
    print(f'{sN}: full matrix {full.shape}, nnz={full.nnz}')
print('DONE')
