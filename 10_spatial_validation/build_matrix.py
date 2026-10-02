# Build per-sample full (dense-in-logic) 67-gene count matrix from GSE208253 h5 files.
# Output: one .mtx (Matrix Market, genes x spots) + genes.tsv + spots.tsv per sample.
import h5py, numpy as np, os, gzip, csv, json
from scipy import sparse
from scipy.io import mmwrite

BASE = 'D:/GSE208253_spatial/GSE208253_processed'
OUT  = 'D:/GSE208253_spatial/spatial_analysis/input'
os.makedirs(OUT, exist_ok=True)

GENES = [g.strip() for g in open('D:/GSE208253_spatial/_67gene_extract/67gene_list.txt', encoding='utf-8') if g.strip()]
assert len(GENES) == 67, len(GENES)

SAMPLES = [(f's{i}', f'GSM{6339630+i}') for i in range(1, 13)]

manifest = []
for sN, gsm in SAMPLES:
    h5p = f'{BASE}/{gsm}_{sN}/filtered_feature_bc_matrix.h5'
    with h5py.File(h5p, 'r') as f:
        shape = tuple(f['matrix/shape'][()])          # [n_genes, n_spots]
        indptr  = f['matrix/indptr'][()]
        indices = f['matrix/indices'][()]
        data    = f['matrix/data'][()]
        barcodes = [b.decode() for b in f['matrix/barcodes'][()]]
        feat_names = [x.decode() for x in f['matrix/features/name'][()]]
    n_genes, n_spots = int(shape[0]), int(shape[1])
    assert n_spots == len(barcodes), (n_spots, len(barcodes))
    # CSC: indptr indexed by spot(col); indices = gene rows
    gene_row = {g: i for i, g in enumerate(feat_names)}
    rows = [gene_row[g] for g in GENES]           # 67 gene row indices
    rowset = set(rows)
    gpos = {r: k for k, r in enumerate(rows)}
    # extract 67-gene submatrix (67 x n_spots), keep all spots (zeros included implicitly)
    sub = sparse.lil_matrix((67, n_spots), dtype=np.int32)
    for col in range(n_spots):
        st, en = int(indptr[col]), int(indptr[col+1])
        for r_, v_ in zip(indices[st:en], data[st:en]):
            if r_ in rowset:
                sub[gpos[r_], col] = int(v_)
    sub = sub.tocsr()
    mmwrite(f'{OUT}/{sN}_67gene.mtx', sub)
    with open(f'{OUT}/{sN}_genes.tsv','w',encoding='utf-8') as fh:
        fh.write('\n'.join(GENES))
    with open(f'{OUT}/{sN}_spots.tsv','w',encoding='utf-8') as fh:
        fh.write('\n'.join(barcodes))
    nnz = int(sub.nnz)
    manifest.append(dict(sample=sN, gsm=gsm, n_genes=67, n_spots=n_spots,
                         nnz_67gene=nnz, sparsity=round(1-nnz/(67*n_spots),4)))
    print(f'{sN} {gsm}: spots={n_spots} nnz={nnz} sparsity={1-nnz/(67*n_spots):.4f}')

json.dump(manifest, open(f'{OUT}/_matrix_manifest.json','w'), indent=1)
print('DONE ->', OUT)
