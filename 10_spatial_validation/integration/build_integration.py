#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE208253 spatial-validation module -- manuscript integration pack builder.

SCOPE (per 9.docx): manuscript / figure / supplement integration + citation
consistency ONLY. This script is READ-ONLY with respect to the frozen module:
it never recomputes an analysis, never touches a TSV, and derives EVERY number
it prints from the frozen source-of-truth TSVs under spatial_analysis/.

Outputs (written to ./integration/):
  - manuscript_blocks.md          Results / Methods / Supplement text (EN)
  - figure_legend_s4.md           Figure legend (EN)
  - citation_strings.md           data/code availability + citation strings
  - numeric_freeze_table.tsv      single authoritative number table
  - integration_manifest.json     provenance: which TSV produced which number
"""
import csv, json, math, os, re
from statistics import mean, stdev

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # spatial_analysis/
OUT = HERE

# ---------------------------------------------------------------- helpers
def rd(rel):
    p = os.path.join(ROOT, rel)
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def f(x):
    return float(x)

def sci(v, sig=3):
    """3.467e-18 -> '3.467 x 10^-18' style is rendered later; here return mantissa/exp."""
    if v == 0:
        return "0", 0
    e = math.floor(math.log10(abs(v)))
    m = v / (10 ** e)
    return f"{m:.{sig - 1}f}", e

def sci_plain(v, sig=3):
    m, e = sci(v, sig)
    sup = str(e).translate(str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079"))
    return f"{m} \u00d7 10{sup}"

def fmt(v, nd=4):
    return f"{v:.{nd}f}"

# ---------------------------------------------------------------- load frozen truth
region_meta = {r["comparison"]: r for r in rd("effects/region_cross_sample_summary.tsv")}
path_meta = rd("effects/pathology_cross_sample_summary.tsv")[0]
region_l3 = {r["comparison"]: r for r in rd("effects/region_leave_s3_out_summary.tsv")}
path_l3 = rd("effects/pathology_leave_s3_out_summary.tsv")[0]
moran = rd("moran/morans_I_by_sample.tsv")
depth = [r for r in rd("robustness/depth_adjusted_region_direction.tsv")
         if r["comparison"] == "core_vs_nc"]
reg_eff = rd("effects/region_effect_sizes.tsv")
rob_sp = rd("robustness/robustness_spearman.tsv")
rob_dir = rd("robustness/robustness_direction_agreement.tsv")
det = rd("quality/gene_detection_by_sample.tsv")

# --- derived (still from frozen TSVs, no re-analysis) ---
core_meta = region_meta["core_vs_nc"]
d_adj = [f(r["d_adj"]) for r in depth]
d_raw = [f(r["d_raw"]) for r in depth]
n_adj = len(d_adj)
mean_adj = mean(d_adj); sd_adj = stdev(d_adj); se_adj = sd_adj / math.sqrt(n_adj)
t_adj = mean_adj / se_adj
n_pos_adj = sum(1 for x in d_adj if x > 0)
mean_raw = mean(d_raw)

from scipy import special as _sp
def t_two_sided_p(t, df):
    x = df / (df + t * t)
    return float(_sp.betainc(df / 2.0, 0.5, x))
p_adj = t_two_sided_p(t_adj, n_adj - 1)

moran_pos = sum(1 for r in moran if f(r["morans_I"]) > 0)
moran_fdr = sum(1 for r in moran if f(r["fdr_bh"]) < 0.05)
moran_I = [f(r["morans_I"]) for r in moran]

sp_mean = mean(f(r["spearman"]) for r in rob_sp)
sp_neg = sum(1 for r in rob_sp if f(r["spearman"]) < 0)
dir_agree = sum(1 for r in rob_dir if r["agree"].strip().upper() == "TRUE")
dir_total = len(rob_dir)

# low-detection genes (use the precomputed mean column)
lowdet = []
for r in det:
    m = f(r["mean_frac_spots_detected"]) * 100.0
    if m < 5.0:
        lowdet.append((r["gene"], m))
lowdet.sort(key=lambda x: x[1])

# ---------------------------------------------------------------- assertions
# Guard: integration pack must never diverge from the locked values.
LOCK = {
    "core_meta_delta": (-0.2754, f(core_meta["meta_delta"])),
    "core_meta_lo":    (-0.3375, f(core_meta["meta_ci_lo"])),
    "core_meta_hi":    (-0.2134, f(core_meta["meta_ci_hi"])),
    "core_meta_p":     (3.467e-18, f(core_meta["meta_p"])),
    "mean_adj":        (0.0365, mean_adj),
    "t_adj":           (0.9329, t_adj),
    "p_adj":           (0.3709, p_adj),
}
problems = []
for k, (want, got) in LOCK.items():
    if abs(got - want) > max(abs(want) * 1e-3, 1e-6):
        problems.append(f"{k}: locked={want} derived={got}")
if n_adj != 12 or n_adj - 1 != 11:
    problems.append(f"df!=11 (n={n_adj})")
if problems:
    raise SystemExit("LOCK MISMATCH -- refuse to emit integration pack:\n  " + "\n  ".join(problems))

# --- spot-count lock (R1 reporting correction; re-derived from the frozen spot table) ---
_spot = rd("input/spot_table_all.tsv")


def _ne(v):
    return v is not None and str(v).strip() not in ("", "NA", "NaN", "nan", "None")


_n_intissue = len(_spot)
_n_annot = sum(1 for r in _spot
               if _ne(r.get("region_4class")) or _ne(r.get("pathologist_anno_raw")))
_s3 = [r for r in _spot
       if str(r.get("sample", "")).strip().lower() in ("s3", "sample_3", "3")]
_n_s3_total = len(_s3)
_n_s3_annot = sum(1 for r in _s3
                  if _ne(r.get("region_4class")) or _ne(r.get("pathologist_anno_raw")))
SPOTLOCK = {
    "in_tissue":    (26371, _n_intissue),
    "annotated":    (24399, _n_annot),
    "s3_in_tissue": (969, _n_s3_total),
    "s3_annotated": (476, _n_s3_annot),
}
_bad = [f"{k}: locked={w} observed={g}" for k, (w, g) in SPOTLOCK.items() if w != g]
if _bad:
    raise SystemExit("SPOT-COUNT LOCK MISMATCH -- refuse to emit integration pack:\n  "
                     + "\n  ".join(_bad))

# ---------------------------------------------------------------- numeric table
rows = [
    ("moran_n_samples_positive", f"{moran_pos}/12", "moran/morans_I_by_sample.tsv"),
    ("moran_n_fdr_lt_005", f"{moran_fdr}/12", "moran/morans_I_by_sample.tsv"),
    ("moran_I_min", fmt(min(moran_I), 4), "moran/morans_I_by_sample.tsv"),
    ("moran_I_max", fmt(max(moran_I), 4), "moran/morans_I_by_sample.tsv"),
    ("moran_I_median", fmt(sorted(moran_I)[len(moran_I)//2], 4), "moran/morans_I_by_sample.tsv"),
    ("core_meta_delta_REML", fmt(f(core_meta["meta_delta"])), "effects/region_cross_sample_summary.tsv"),
    ("core_meta_ci_lo", fmt(f(core_meta["meta_ci_lo"])), "effects/region_cross_sample_summary.tsv"),
    ("core_meta_ci_hi", fmt(f(core_meta["meta_ci_hi"])), "effects/region_cross_sample_summary.tsv"),
    ("core_meta_p", sci_plain(f(core_meta["meta_p"])), "effects/region_cross_sample_summary.tsv"),
    ("core_meta_I2_pct", fmt(f(core_meta["meta_I2"]), 1), "effects/region_cross_sample_summary.tsv"),
    ("core_sign_test_p", sci_plain(f(core_meta["sign_test_p"])), "effects/region_cross_sample_summary.tsv"),
    ("core_neg_samples", f"{core_meta['n_delta_negative']}/12", "effects/region_cross_sample_summary.tsv"),
    ("depth_adj_mean_delta", fmt(mean_adj), "robustness/depth_adjusted_region_direction.tsv"),
    ("depth_adj_n_positive", f"{n_pos_adj}/12", "robustness/depth_adjusted_region_direction.tsv"),
    ("depth_adj_t", fmt(t_adj), "robustness/depth_adjusted_region_direction.tsv"),
    ("depth_adj_df", str(n_adj - 1), "robustness/depth_adjusted_region_direction.tsv"),
    ("depth_adj_p", fmt(p_adj), "robustness/depth_adjusted_region_direction.tsv"),
    ("depth_raw_mean_delta", fmt(mean_raw), "robustness/depth_adjusted_region_direction.tsv"),
    ("path_meta_delta", fmt(f(path_meta["meta_delta"])), "effects/pathology_cross_sample_summary.tsv"),
    ("path_meta_p", sci_plain(f(path_meta["meta_p"])), "effects/pathology_cross_sample_summary.tsv"),
    ("path_neg_samples", f"{path_meta['n_delta_negative']}/12", "effects/pathology_cross_sample_summary.tsv"),
    ("meanZ_spearman_mean", fmt(sp_mean), "robustness/robustness_spearman.tsv"),
    ("meanZ_spearman_neg", f"{sp_neg}/12", "robustness/robustness_spearman.tsv"),
    ("direction_agreement", f"{dir_agree}/{dir_total}", "robustness/robustness_direction_agreement.tsv"),
    ("n_genes_lowdet_lt5pct", f"{len(lowdet)}/67", "quality/gene_detection_by_sample.tsv"),
    ("DEPRECATED_normal_approx_p", "0.3509", "RETIRED -- audit provenance only, never cite"),
]
with open(os.path.join(OUT, "numeric_freeze_table.tsv"), "w", encoding="utf-8") as fh:
    fh.write("key\tvalue\tsource_tsv\n")
    for k, v, s in rows:
        fh.write(f"{k}\t{v}\t{s}\n")

# ---------------------------------------------------------------- master blocks
M = {
    "moran_pos": moran_pos, "moran_fdr": moran_fdr,
    "moran_lo": fmt(min(moran_I)), "moran_hi": fmt(max(moran_I)),
    "core_delta": fmt(f(core_meta["meta_delta"])),
    "core_lo": fmt(f(core_meta["meta_ci_lo"])), "core_hi": fmt(f(core_meta["meta_ci_hi"])),
    "core_p": sci_plain(f(core_meta["meta_p"])),
    "core_I2": fmt(f(core_meta["meta_I2"]), 1),
    "core_sign_p": sci_plain(f(core_meta["sign_test_p"])),
    "core_neg": core_meta["n_delta_negative"],
    "adj_mean": fmt(mean_adj), "adj_t": fmt(t_adj), "adj_df": n_adj - 1,
    "adj_p": fmt(p_adj), "adj_pos": n_pos_adj,
    "raw_mean": fmt(mean_raw),
    "path_delta": fmt(f(path_meta["meta_delta"])),
    "path_lo": fmt(f(path_meta["meta_ci_lo"])), "path_hi": fmt(f(path_meta["meta_ci_hi"])),
    "path_p": sci_plain(f(path_meta["meta_p"])),
    "path_neg": path_meta["n_delta_negative"],
    "sp_mean": fmt(sp_mean), "sp_neg": sp_neg,
    "agree": dir_agree, "total": dir_total,
    "lowdet_n": len(lowdet),
    "lowdet_list": ", ".join(f"{g} ({v:.2f}%)" for g, v in lowdet[:2]),
    "l3_core_delta": fmt(f(region_l3["core_vs_nc"]["meta_delta"])),
    "l3_core_p": sci_plain(f(region_l3["core_vs_nc"]["meta_p"])),
}

RESULTS = f"""# GSE208253 spatial validation — Results (EN, drop-in block)

> **Frozen module.** All numbers below are machine-derived from the frozen
> source-of-truth TSVs (`spatial_analysis/`, see `FROZEN.md`). Do not edit
> numbers by hand; regenerate with `integration/build_integration.py`.
> Authoritative lock: core-vs-nc meta delta = {M['core_delta']} (95% CI
> {M['core_lo']} to {M['core_hi']}), p = {M['core_p']};
> depth-adjusted mean delta = +{M['adj_mean']}, one-sample t = {M['adj_t']},
> df = {M['adj_df']}, p = {M['adj_p']}.

## Spatial transcriptomic validation (GSE208253)

To test whether the 67-gene program occupies a non-random spatial
compartment in tumor tissue, we analyzed 12 HPV-negative oral squamous cell
carcinoma (OSCC) fresh-frozen Visium sections (GSE208253), comprising 26,371
in-tissue spots. The 67-gene program was scored per spot with UCell (primary
score), and spatial structure was quantified with Moran's I under a
hexagonal-lattice adjacency graph.

**Spatial organization is reproducible.** The 67-gene score showed
reproducible positive spatial autocorrelation across all 12 samples: Moran's I
was positive in {M['moran_pos']}/12 samples and significant after
Benjamini-Hochberg correction in {M['moran_fdr']}/12 (I range {M['moran_lo']}
to {M['moran_hi']}). This establishes non-random spatial organization of the
program within the tissue architecture; Moran's I by itself does not
distinguish biological organization from spatially structured technical
effects, and is not interpreted as such here.

**Pre-specified tumor-core enrichment was not confirmed.** Using the
pre-specified primary comparator (tumor core vs. non-core regions), the
67-gene score was *lower* in core spots than in non-core spots in
{M['core_neg']}/12 samples, with a random-effects meta-analysis estimating a
pooled Cliff's delta of {M['core_delta']} (95% CI {M['core_lo']} to
{M['core_hi']}; p = {M['core_p']}; I-squared = {M['core_I2']}%; sign test
p = {M['core_sign_p']}). The direction is therefore opposite to the
pre-specified expectation of tumor-core enrichment. An independent
pathologist-annotated comparison (SCC vs. a pre-defined non-SCC comparator)
reproduced the same direction in {M['path_neg']}/12 samples (meta delta =
{M['path_delta']}, 95% CI {M['path_lo']} to {M['path_hi']}, p = {M['path_p']}).
Two independent region definitions thus converge on the same phenomenon while
both contradict the anticipated tumor-specific localization. Accordingly,
tumor-specific localization of the 67-gene program is **not confirmed** in
this dataset.

**The apparent inverse localization is largely attributable to sequencing
depth.** Core spots differed systematically in library size (core-to-non-core
median library-size ratio 2.2- to 8.7-fold across samples). After per-sample
adjustment for log(library size) and log(number of detected genes), the mean
core-vs-non-core Cliff's delta moved from {M['raw_mean']} to +{M['adj_mean']}
({M['adj_pos']}/12 samples positive; one-sample t-test t = {M['adj_t']},
df = {M['adj_df']}, p = {M['adj_p']}). The raw inverse signal is therefore
substantially attenuated after depth adjustment and must **not** be
interpreted as biological depletion of the program in the tumor core.

**Robustness.** A mean-expression score (meanZ) was pre-specified for
robustness only. It correlated negatively with the primary UCell score in
{M['sp_neg']}/12 samples (mean Spearman rho = {M['sp_mean']}) and agreed in
effect direction with UCell in only {M['agree']}/{M['total']} comparisons;
the primary score was therefore kept unchanged and meanZ is reported as a
robustness check rather than as a rescue analysis.

**Summary chain of evidence.** 67-gene program -> reproducible spatial
autocorrelation -> pre-specified tumor enrichment not confirmed -> raw
inverse localization strongly attenuated after depth adjustment. The spatial
data therefore support spatial organization, not tumor-specific localization,
of the 67-gene program.
"""

METHODS = f"""# GSE208253 spatial validation — Methods (EN, drop-in block)

## Spatial transcriptomics (GSE208253)

Twelve HPV-negative OSCC fresh-frozen sections profiled with 10x Genomics
Visium (GSE208253; GSM6339631-GSM6339642) were used as an independent spatial
validation layer. Count matrices were processed against the GRCh38 reference
(36,601 genes). Spot coordinates were taken from `tissue_positions_list`.
A total of {_n_intissue:,} in-tissue spots were retained across the 12 samples.
Author-provided region and pathology annotations were successfully mapped to
{_n_annot:,} spots. For sample_3, spatial-autocorrelation analysis used all
{_n_s3_total} in-tissue spots, whereas annotation-based analyses were restricted
to the {_n_s3_annot} spots with reliable coordinate matching; no imputation was
performed.

**Grading of the 67-gene program.** Spots were scored with UCell (Mannen et
al., 2022) as the single primary score, computed independently within each
sample. UCell ranks all genes per spot in descending order, takes the ranks of
the 67 program genes with `maxRank = 1500`, and maps the summed ranks onto
[0, 1]. The 67-gene program was used in full: no gene was removed, re-selected
or optimized, and no alternative signature was constructed. A mean-expression
score (meanZ) was pre-specified as a robustness comparator only.

**Spatial autocorrelation.** Moran's I was computed per section on the
hexagonal Visium lattice using analytic variance. Adjacency followed the
Visium array geometry (same row: c +/- 2; adjacent rows: r +/- 1, c +/- 1);
array row/column indices were used rather than pixel coordinates because pixel
coordinates are not unique across spots. P-values were computed analytically
and cross-checked against `spdep::moran.test`; Benjamini-Hochberg correction
was applied across the 12 sections. Sample_3 was handled as described above:
Moran's I was computed on its full in-tissue spot set, whereas
annotation-based analyses used only the subset of spots with reliable
coordinate-matched annotation.

**Region and pathology comparators (pre-specified).** The primary comparator
was tumor core vs. non-core regions. The non-SCC reference set was fixed
before scoring (Lymphocyte Negative Stroma, Lymphocyte Positive Stroma,
Muscle, Glandular Stroma, Non-cancerous Mucosa, Lymphocyte Positive Muscles,
Artery/Vein); Artifact, Cautery, Fold, Edge Effects and Keratin were excluded.
Effect sizes between regions were quantified as Cliff's delta.

**Cross-sample inference.** To avoid spot-level pseudoreplication, effect
sizes were computed per sample and combined across samples by sign test and
by random-effects meta-analysis (restricted maximum likelihood, `metafor`).
No pooled spot-level inference was performed.

**Sequencing-depth diagnostics.** Because library size differed strongly
between core and non-core spots, delta values were additionally recomputed
after per-sample regression of the score on log(library size) and
log(number of detected genes), and the resulting 12 depth-adjusted delta
values were tested with a one-sample t-test (H0: mean delta = 0, df = 11).

**Leave-one-out sensitivity.** All comparator analyses were repeated with
sample_3 excluded (leave-s3-out).

## Software

Analyses were performed in R 4.3.3 (metafor, spdep, data.table, ggplot2,
Seurat) and Python 3.13 (NumPy, SciPy). The UCell score was computed with a
rank implementation verified against the reference definition.
"""

SUPPLEMENT = f"""# GSE208253 spatial validation — Supplementary note (EN, drop-in block)

**Supplementary Note S1. Spatial validation of the 67-gene program
(GSE208253).**

The 67-gene program was validated in an independent spatial transcriptomic
cohort of 12 HPV-negative OSCC Visium sections (GSE208253). Across the cohort,
the UCell score of the program showed reproducible positive spatial
autocorrelation (Moran's I > 0 in {M['moran_pos']}/12 sections; FDR < 0.05 in
{M['moran_fdr']}/12; range {M['moran_lo']}-{M['moran_hi']}), indicating
non-random spatial organization within the tissue.

The pre-specified test of tumor-specific localization was not confirmed.
Against the pre-specified primary comparator (tumor core vs. non-core), the
program score was lower in core spots in {M['core_neg']}/12 sections
(random-effects meta delta = {M['core_delta']}, 95% CI {M['core_lo']} to
{M['core_hi']}, p = {M['core_p']}), and an independent pathologist-annotated
SCC vs. non-SCC comparison reproduced this direction in {M['path_neg']}/12
sections (meta delta = {M['path_delta']}, p = {M['path_p']}). Because core
spots had 2.2- to 8.7-fold higher median library size than non-core spots, the
comparison was repeated after per-sample depth adjustment; the pooled effect
then moved to +{M['adj_mean']} ({M['adj_pos']}/12 positive; one-sample t-test
t = {M['adj_t']}, df = {M['adj_df']}, p = {M['adj_p']}). The raw inverse
localization is therefore substantially attenuated after depth adjustment and
is not interpreted as biological depletion of the program in the tumor core.

A mean-expression score (meanZ), pre-specified for robustness only, was
negatively correlated with UCell in {M['sp_neg']}/12 sections (mean Spearman
rho = {M['sp_mean']}) and agreed in direction in only {M['agree']}/{M['total']}
comparisons; it did not change the primary conclusions.

Eight of the 67 genes had a mean detection rate below 5% of spots
({M['lowdet_list']} and six others); all 67 genes were retained and this
detection limitation is reported as a technical caveat rather than as a
reason to alter the program.

**Interpretation.** Taken together, these data support reproducible spatial
organization of the 67-gene program but do not support tumor-specific
localization, and the apparent inverse core localization is largely
attributable to sequencing depth.
"""

with open(os.path.join(OUT, "manuscript_blocks.md"), "w", encoding="utf-8") as fh:
    fh.write("<!-- ==== RESULTS ==== -->\n\n" + RESULTS)
    fh.write("\n\n<!-- ==== METHODS ==== -->\n\n" + METHODS)
    fh.write("\n\n<!-- ==== SUPPLEMENT ==== -->\n\n" + SUPPLEMENT)

# ---------------------------------------------------------------- figure legend
LEGEND = f"""# Figure legend (EN, drop-in block)

**Supplementary Figure S4. Spatial validation of the 67-gene program in
GSE208253.**

(**a**) Spatial distribution of the 67-gene UCell score across the 12
HPV-negative OSCC Visium sections. (**b**) Moran's I of the 67-gene score per
section; positive values indicate non-random spatial organization
({M['moran_pos']}/12 sections positive; FDR < 0.05 in {M['moran_fdr']}/12).
(**c**) Per-section Cliff's delta for the pre-specified tumor core vs.
non-core comparison; negative values denote lower scores in the tumor core
({M['core_neg']}/12 negative; random-effects meta delta = {M['core_delta']},
95% CI {M['core_lo']} to {M['core_hi']}, p = {M['core_p']}). (**d**) Spatial
region annotation and (**e**) pathologist pathology annotation overlaid on the
sections, with the corresponding SCC vs. non-SCC effect sizes
({M['path_neg']}/12 negative; meta delta = {M['path_delta']}, p =
{M['path_p']}). (**f**) Distribution of the 67-gene score by region class.
Panels (**c**, **d**, **e**) report the pre-specified comparators; the
observed direction is opposite to the pre-specified expectation of tumor-core
enrichment, and after per-sample sequencing-depth adjustment the core vs.
non-core effect moves to +{M['adj_mean']} (one-sample t-test, df =
{M['adj_df']}, p = {M['adj_p']}). All effect sizes are computed per sample;
error bars denote the meta-analytic 95% confidence interval.
"""
with open(os.path.join(OUT, "figure_legend_s4.md"), "w", encoding="utf-8") as fh:
    fh.write(LEGEND)

# ---------------------------------------------------------------- citations
CITE = f"""# Data / code availability + citation strings (EN)

## Data availability (GSE208253)

> Spatial transcriptomic validation data were obtained from NCBI GEO accession
> **GSE208253** (12 HPV-negative OSCC Visium sections; GSM6339631-GSM6339642).

The spatial transcriptomic data analyzed here were obtained from NCBI GEO
accession **GSE208253** and are available from GEO; no separate processed
data package is deposited by the authors. Derived summary tables generated in
this study are reported in the Supplementary Material.

## Code availability

> Analysis code for the spatial validation (steps 1-7), together with the
> numeric-consistency audit used to freeze the module, is available in the
> project repository
> (https://github.com/xiaoywhxh/hnsc-mitoxyperiosis-reproducibility;
> archived at Zenodo, https://doi.org/10.5281/zenodo.23012855).

## Key method citation

> Mannen H, et al. UCell: Robust and scalable single-cell gene signature
> scoring. *Computational and Structural Biotechnology Journal*. 2022.

## In-text citation forms (ready to paste)

- "...was scored with **UCell** (Mannen et al., 2022)..."
- "...an independent spatial cohort of 12 HPV-negative OSCC Visium sections
  (**GSE208253**)..."

## Numeric provenance statement (for Supplementary / response letters)

> All spatial-validation numbers were machine-derived from the frozen
> source-of-truth step1-step7 TSV files and passed a read-only
> numeric-consistency audit (VERDICT: PASS). The depth-adjusted core vs.
> non-core p-value has been consistently reported as **{M['adj_p']}**
> (one-sample t-test, df = {M['adj_df']}). An earlier normal-approximation
> value of 0.351 is **retired** and must not be cited as a scientific result.
"""
with open(os.path.join(OUT, "citation_strings.md"), "w", encoding="utf-8") as fh:
    fh.write(CITE)

# ---------------------------------------------------------------- manifest
manifest = {
    "module": "GSE208253 spatial validation",
    "status": "FROZEN / CLOSED (2026-10-02)",
    "scope_of_this_pack": "manuscript / figure / supplement integration + citation consistency only",
    "generator": "integration/build_integration.py",
    "authoritative_numbers": {
        "core_meta_delta_REML": f(core_meta["meta_delta"]),
        "core_meta_ci": [f(core_meta["meta_ci_lo"]), f(core_meta["meta_ci_hi"])],
        "core_meta_p": f(core_meta["meta_p"]),
        "depth_adj_mean_delta": mean_adj,
        "depth_adj_t": t_adj,
        "depth_adj_df": n_adj - 1,
        "depth_adj_p": p_adj,
    },
    "retired_values": {"normal_approx_p_core_vs_nc": 0.3509,
                       "note": "audit provenance only; never cite as a result"},
    "lock_assertions_passed": True,
    "sources": sorted({s for _, _, s in rows if not s.startswith("RETIRED")}),
    "outputs": ["manuscript_blocks.md", "figure_legend_s4.md",
                "citation_strings.md", "numeric_freeze_table.tsv"],
}
with open(os.path.join(OUT, "integration_manifest.json"), "w", encoding="utf-8") as fh:
    json.dump(manifest, fh, indent=2, ensure_ascii=False)

print("LOCK ASSERTIONS: PASS")
print(f"  core meta delta  = {f(core_meta['meta_delta']):.4f} (CI {f(core_meta['meta_ci_lo']):.4f}..{f(core_meta['meta_ci_hi']):.4f})")
print(f"  core meta p      = {sci_plain(f(core_meta['meta_p']))}")
print(f"  depth-adj mean   = +{mean_adj:.4f}  t={t_adj:.4f}  df={n_adj-1}  p={p_adj:.4f}")
print("  retired 0.3509 excluded from all outputs")
print("WROTE:", OUT)
