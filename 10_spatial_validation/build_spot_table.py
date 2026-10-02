# Assemble analysis-ready inputs: per-sample spot table with counts(67) handled by R,
# plus coordinates + author annotation joined on (sample, barcode).
import gzip, csv, json, os
from collections import defaultdict

BASE = 'D:/GSE208253_spatial/GSE208253_processed'
OUT  = 'D:/GSE208253_spatial/spatial_analysis/input'
os.makedirs(OUT, exist_ok=True)
SAMPLES = [(f's{i}', f'GSM{6339630+i}') for i in range(1, 13)]

# author annotation keyed by (sample, barcode)  [GEO barcode]
ann4 = {}; annraw = {}
with open('D:/GSE208253_spatial/spatial_prep/author_region_annotation.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        ann4[(r['GEO_sample'], r['GEO_barcode'])] = r['region_4class']
        annraw[(r['GEO_sample'], r['GEO_barcode'])] = r['pathologist_anno_raw']

rows_all = []
for sN, gsm in SAMPLES:
    p = f'{BASE}/{gsm}_{sN}/spatial/tissue_positions_list.csv.gz'
    with gzip.open(p, 'rt') as f:
        for line in f:
            parts = line.rstrip('\n').split(',')
            bc, it, ar, ac, pr, pc = parts[0], int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])
            if it != 1:      # keep only in-tissue
                continue
            rows_all.append([sN, gsm, bc, ar, ac, pr, pc,
                             ann4.get((sN, bc), ''), annraw.get((sN, bc), '')])

hdr = ['sample','GSM','barcode','array_row','array_col','pxl_row','pxl_col',
       'region_4class','pathologist_anno_raw']
with open(f'{OUT}/spot_table_all.tsv','w',encoding='utf-8',newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(hdr); w.writerows(rows_all)

# summary
per = defaultdict(lambda: dict(n=0, ann=0))
for r in rows_all:
    per[r[0]]['n'] += 1
    if r[7]: per[r[0]]['ann'] += 1
print('sample  in_tissue  annotated')
for s in sorted(per, key=lambda x:int(x[1:])):
    print(f'  {s:<6} {per[s]["n"]:>8} {per[s]["ann"]:>10}')
print('TOTAL in_tissue', sum(v['n'] for v in per.values()),
      ' annotated', sum(v['ann'] for v in per.values()))
json.dump({k: v for k, v in per.items()}, open(f'{OUT}/_spot_table_summary.json','w'), indent=1)
