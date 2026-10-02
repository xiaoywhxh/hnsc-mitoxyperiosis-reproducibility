<!-- ==== RESULTS ==== -->

# GSE208253 spatial validation — Results (EN, drop-in block)

> **Frozen module.** All numbers below are machine-derived from the frozen
> source-of-truth TSVs (`spatial_analysis/`, see `FROZEN.md`). Do not edit
> numbers by hand; regenerate with `integration/build_integration.py`.
> Authoritative lock: core-vs-nc meta delta = -0.2754 (95% CI
> -0.3375 to -0.2134), p = 3.47 × 10⁻¹⁸;
> depth-adjusted mean delta = +0.0365, one-sample t = 0.9329,
> df = 11, p = 0.3709.

## Spatial transcriptomic validation (GSE208253)

To test whether the 67-gene program occupies a non-random spatial
compartment in tumor tissue, we analyzed 12 HPV-negative oral squamous cell
carcinoma (OSCC) fresh-frozen Visium sections (GSE208253), comprising 26,371
in-tissue spots. The 67-gene program was scored per spot with UCell (primary
score), and spatial structure was quantified with Moran's I under a
hexagonal-lattice adjacency graph.

**Spatial organization is reproducible.** The 67-gene score showed
reproducible positive spatial autocorrelation across all 12 samples: Moran's I
was positive in 12/12 samples and significant after
Benjamini-Hochberg correction in 12/12 (I range 0.0571
to 0.2304). This establishes non-random spatial organization of the
program within the tissue architecture; Moran's I by itself does not
distinguish biological organization from spatially structured technical
effects, and is not interpreted as such here.

**Pre-specified tumor-core enrichment was not confirmed.** Using the
pre-specified primary comparator (tumor core vs. non-core regions), the
67-gene score was *lower* in core spots than in non-core spots in
12/12 samples, with a random-effects meta-analysis estimating a
pooled Cliff's delta of -0.2754 (95% CI -0.3375 to
-0.2134; p = 3.47 × 10⁻¹⁸; I-squared = 77.0%; sign test
p = 4.88 × 10⁻⁴). The direction is therefore opposite to the
pre-specified expectation of tumor-core enrichment. An independent
pathologist-annotated comparison (SCC vs. a pre-defined non-SCC comparator)
reproduced the same direction in 10/12 samples (meta delta =
-0.1327, 95% CI -0.1907 to -0.0747, p = 7.22 × 10⁻⁶).
Two independent region definitions thus converge on the same phenomenon while
both contradict the anticipated tumor-specific localization. Accordingly,
tumor-specific localization of the 67-gene program is **not confirmed** in
this dataset.

**The apparent inverse localization is largely attributable to sequencing
depth.** Core spots differed systematically in library size (core-to-non-core
median library-size ratio 2.2- to 8.7-fold across samples). After per-sample
adjustment for log(library size) and log(number of detected genes), the mean
core-vs-non-core Cliff's delta moved from -0.2869 to +0.0365
(8/12 samples positive; one-sample t-test t = 0.9329,
df = 11, p = 0.3709). The raw inverse signal is therefore
substantially attenuated after depth adjustment and must **not** be
interpreted as biological depletion of the program in the tumor core.

**Robustness.** A mean-expression score (meanZ) was pre-specified for
robustness only. It correlated negatively with the primary UCell score in
12/12 samples (mean Spearman rho = -0.2442) and agreed in
effect direction with UCell in only 5/48 comparisons;
the primary score was therefore kept unchanged and meanZ is reported as a
robustness check rather than as a rescue analysis.

**Summary chain of evidence.** 67-gene program -> reproducible spatial
autocorrelation -> pre-specified tumor enrichment not confirmed -> raw
inverse localization strongly attenuated after depth adjustment. The spatial
data therefore support spatial organization, not tumor-specific localization,
of the 67-gene program.


<!-- ==== METHODS ==== -->

# GSE208253 spatial validation — Methods (EN, drop-in block)

## Spatial transcriptomics (GSE208253)

Twelve HPV-negative OSCC fresh-frozen sections profiled with 10x Genomics
Visium (GSE208253; GSM6339631-GSM6339642) were used as an independent spatial
validation layer. Count matrices were processed against the GRCh38 reference
(36,601 genes). Spot coordinates were taken from `tissue_positions_list`.
A total of 26,371 in-tissue spots were retained across the 12 samples.
Author-provided region and pathology annotations were successfully mapped to
24,399 spots. For sample_3, spatial-autocorrelation analysis used all
969 in-tissue spots, whereas annotation-based analyses were restricted
to the 476 spots with reliable coordinate matching; no imputation was
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


<!-- ==== SUPPLEMENT ==== -->

# GSE208253 spatial validation — Supplementary note (EN, drop-in block)

**Supplementary Note S1. Spatial validation of the 67-gene program
(GSE208253).**

The 67-gene program was validated in an independent spatial transcriptomic
cohort of 12 HPV-negative OSCC Visium sections (GSE208253). Across the cohort,
the UCell score of the program showed reproducible positive spatial
autocorrelation (Moran's I > 0 in 12/12 sections; FDR < 0.05 in
12/12; range 0.0571-0.2304), indicating
non-random spatial organization within the tissue.

The pre-specified test of tumor-specific localization was not confirmed.
Against the pre-specified primary comparator (tumor core vs. non-core), the
program score was lower in core spots in 12/12 sections
(random-effects meta delta = -0.2754, 95% CI -0.3375 to
-0.2134, p = 3.47 × 10⁻¹⁸), and an independent pathologist-annotated
SCC vs. non-SCC comparison reproduced this direction in 10/12
sections (meta delta = -0.1327, p = 7.22 × 10⁻⁶). Because core
spots had 2.2- to 8.7-fold higher median library size than non-core spots, the
comparison was repeated after per-sample depth adjustment; the pooled effect
then moved to +0.0365 (8/12 positive; one-sample t-test
t = 0.9329, df = 11, p = 0.3709). The raw inverse
localization is therefore substantially attenuated after depth adjustment and
is not interpreted as biological depletion of the program in the tumor core.

A mean-expression score (meanZ), pre-specified for robustness only, was
negatively correlated with UCell in 12/12 sections (mean Spearman
rho = -0.2442) and agreed in direction in only 5/48
comparisons; it did not change the primary conclusions.

Eight of the 67 genes had a mean detection rate below 5% of spots
(CASP5 (0.49%), GLS2 (0.60%) and six others); all 67 genes were retained and this
detection limitation is reported as a technical caveat rather than as a
reason to alter the program.

**Interpretation.** Taken together, these data support reproducible spatial
organization of the 67-gene program but do not support tumor-specific
localization, and the apparent inverse core localization is largely
attributable to sequencing depth.
