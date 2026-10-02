# 10_spatial_validation — GSE208253 spatial validation module

**Status: FROZEN / CLOSED (2026-10-02)** — see `FROZEN.md`.

## What this module is

An orthogonal, pre-specified spatial transcriptomic test of the locked 67-gene
mitoxyperiosis-related program in 12 HPV-negative oral squamous cell carcinoma
(OSCC) Visium sections (NCBI GEO accession **GSE208253**), analysed with UCell as
the single primary score.

## Files

| Path | Content |
|---|---|
| `step1_ucell.R`, `step1_compute_scores.py`, `ucell.py` | UCell scoring (rank implementation verified against the reference definition) |
| `build_matrix.py`, `build_full_matrix.py`, `build_spot_table.py` | input assembly from the GEO per-sample matrices |
| `step2_moran.R` | Moran's I on the hexagonal Visium lattice (analytic variance, BH-FDR across sections) |
| `step3_region_effects.R` | pre-specified author region comparator (tumour core vs non-core), Cliff's delta |
| `step4_pathology.R` | independent pathologist-annotated SCC vs non-SCC comparator |
| `step5_robustness.R` | meanZ robustness score and direction agreement |
| `step6_depth_check.R` | sequencing-depth diagnostics (per-sample adjustment; one-sample t-test) |
| `step7_plots.R` | figure generation |
| `fix_fig_ucell.R` | render fix for the UCell grid panel (see below) |
| `audit_numeric_consistency.py` | read-only numeric-consistency audit (VERDICT: PASS) |
| `integration/` | manuscript integration pack (drop-in Results/Methods/Supplement blocks, figure legend, citation strings, numeric freeze table, generator with LOCK assertions, citation-consistency checker) |
| `moran/`, `effects/`, `robustness/`, `quality/` | frozen source-of-truth result tables cited by `FROZEN.md` |
| `plots/` | publication figures (PNG; PDF equivalents omitted) |
| `report/GSE208253_validation_master_tables.md` | master tables (rendered HTML report omitted) |

## Reproducing

Inputs are the public GSE208253 per-sample 10x Visium matrices and the
author-provided region/pathology annotations. Rebuild the local input layer with
`build_spot_table.py` / `build_matrix.py` / `build_full_matrix.py`, then run
`step1`–`step7` in order. `integration/build_integration.py` regenerates the
manuscript blocks; it carries LOCK assertions on both the frozen statistics and
the spot counts and refuses to emit if either drifts.

## Excluded from this repository

`input/` (per-spot matrices), `score/spot_scores.tsv` (per-spot scores) and the
rendered HTML report are derived artefacts that can be regenerated from the
public source data with the scripts above; they are omitted to keep the
repository code-only, consistent with the other modules.

## Note on `fix_fig_ucell.R`

`step7_plots.R` originally mapped the UCell score to `fill` while drawing
`geom_point(shape = 15)`. In ggplot2, shape 15 is a solid point that consumes
`colour`, not `fill`, so the score mapping was silently ignored and every spot was
drawn in the default black. `fix_fig_ucell.R` re-renders that single panel with
`colour`/`scale_colour_viridis_c`. **This is a rendering fix only: the underlying
scores and statistics are unchanged** (the frozen `score/spot_scores.tsv` and
`moran/morans_I_by_sample.tsv` checksums are identical before and after).

## Key frozen results

| Quantity | Value |
|---|---|
| Moran's I positive / FDR<0.05 | 12/12 sections / 12/12 sections (I range 0.0571-0.2304) |
| tumour core vs non-core, meta delta | -0.2754 (95% CI -0.3375 to -0.2134), P = 3.47e-18 |
| SCC vs non-SCC, meta delta | -0.1327 (10/12 negative), P = 7.22e-6 |
| depth-adjusted core vs non-core mean delta | +0.0365 (8/12 positive), one-sample t = 0.9329, df = 11, P = 0.3709 |
| meanZ vs UCell (robustness) | 12/12 negative (mean Spearman rho = -0.2442) |

The 0.3509 normal-approximation value is retired and retained only as audit
provenance; it must not be cited as a scientific result.
