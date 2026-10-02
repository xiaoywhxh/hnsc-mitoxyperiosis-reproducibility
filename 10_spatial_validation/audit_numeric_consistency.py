# -*- coding: utf-8 -*-
"""
read-only / summary-only NUMERIC CONSISTENCY AUDIT
==================================================
Authorised once, before freezing the GSE208253 module (per 8.docx).

What it does
------------
1. RECOMPUTES every headline number directly from the raw step1-step7 TSVs
   (and, where the TSV only holds per-sample rows, re-derives the cross-sample
   summary with the ORIGINAL statistical definition).
2. PARSES the delivered documents (HTML report, master-tables MD) and pulls out
   every headline number that appears there.
3. COMPARES the two, flagging MATCH / MISMATCH / NOT-FOUND and, for each
   mismatch, states which statistic definition each side used.

It WRITES NOTHING except this audit's own report.  It does not touch step1-7
outputs, does not re-select genes, does not change region definitions.

Run:  python audit_numeric_consistency.py
"""
import csv, os, re, json, math, sys

BASE = "D:/GSE208253_spatial/spatial_analysis"
REPORT_HTML = os.path.join(BASE, "report", "GSE208253_spatial_validation_report.html")
REPORT_MD   = os.path.join(BASE, "report", "GSE208253_validation_master_tables.md")
OUT_DIR     = os.path.join(BASE, "audit")
os.makedirs(OUT_DIR, exist_ok=True)

R = []          # audit rows
def rec(item, expected, found, status, note=""):
    R.append(dict(item=item, recomputed=expected, in_document=found,
                  status=status, note=note))

def rd(p, delim="\t"):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delim))

def num(x):
    try: return float(str(x).replace(",", ""))
    except Exception: return float("nan")

# ============================================================================
# PART A -- RECOMPUTE FROM RAW TSV (source of truth)
# ============================================================================
# --- A1. Moran's I -----------------------------------------------------------
moran = rd(os.path.join(BASE, "moran", "morans_I_by_sample.tsv"))
n_moran = len(moran)
n_moran_pos = sum(1 for r in moran if num(r["morans_I"]) > 0)
n_moran_sig = sum(1 for r in moran if num(r["fdr_bh"]) < 0.05)
moran_min = min(num(r["morans_I"]) for r in moran)
moran_max = max(num(r["morans_I"]) for r in moran)
morans = sorted(num(r["morans_I"]) for r in moran)
moran_med = (morans[5] + morans[6]) / 2 if n_moran % 2 == 0 else morans[n_moran//2]
s10 = [r for r in moran if r["sample"] == "s10"][0]
s10_spdep = s10.get("spdep_check_I", "").strip()

rec("Moran's I: n samples computed", n_moran, None, "INFO")
rec("Moran's I: n with I>0", n_moran_pos, None, "INFO")
rec("Moran's I: n with FDR<0.05", n_moran_sig, None, "INFO")
rec("Moran's I: min", round(moran_min, 4), None, "INFO")
rec("Moran's I: max", round(moran_max, 4), None, "INFO")
rec("Moran's I: median", round(moran_med, 4), None, "INFO")
rec("Moran's I: s10 spdep cross-check identical",
    abs(num(s10["morans_I"]) - num(s10_spdep)) < 1e-9 if s10_spdep else None,
    None, "INFO", f"manual={num(s10['morans_I']):.9f} spdep={s10_spdep}")

# --- A2. region effect sizes + cross-sample summary --------------------------
reff = rd(os.path.join(BASE, "effects", "region_effect_sizes.tsv"))
rsum = rd(os.path.join(BASE, "effects", "region_cross_sample_summary.tsv"))
rl3  = rd(os.path.join(BASE, "effects", "region_leave_s3_out_summary.tsv"))

def per_comp(eff, cmp_name, side="A"):
    rows = [r for r in eff if r["comparison"] == cmp_name and str(r["cliffs_delta"]).strip()]
    return rows

core = per_comp(reff, "core_vs_nc")
rec("core_vs_nc: n samples with delta", len(core), None, "INFO")
rec("core_vs_nc: n delta>0", sum(1 for r in core if num(r["cliffs_delta"]) > 0), None, "INFO")
rec("core_vs_nc: n delta<0", sum(1 for r in core if num(r["cliffs_delta"]) < 0), None, "INFO")
rec("core_vs_nc: delta min", round(min(num(r["cliffs_delta"]) for r in core), 4), None, "INFO")
rec("core_vs_nc: delta max", round(max(num(r["cliffs_delta"]) for r in core), 4), None, "INFO")

core_sum = [r for r in rsum if r["comparison"] == "core_vs_nc"][0]
rec("core_vs_nc: meta delta (REML)", round(num(core_sum["meta_delta"]), 4), None, "INFO")
rec("core_vs_nc: meta p", num(core_sum["meta_p"]), None, "INFO")
rec("core_vs_nc: meta CI lo", round(num(core_sum["meta_ci_lo"]), 4), None, "INFO")
rec("core_vs_nc: meta CI hi", round(num(core_sum["meta_ci_hi"]), 4), None, "INFO")
rec("core_vs_nc: meta I2", round(num(core_sum["meta_I2"]), 1), None, "INFO")
rec("core_vs_nc: sign test p", num(core_sum["sign_test_p"]), None, "INFO")

core_l3 = [r for r in rl3 if r["comparison"] == "core_vs_nc"][0]
rec("core_vs_nc leave-s3-out: n", core_l3["n_samples"], None, "INFO")
rec("core_vs_nc leave-s3-out: n delta>0", core_l3["n_delta_positive"], None, "INFO")
rec("core_vs_nc leave-s3-out: meta delta", round(num(core_l3["meta_delta"]), 4), None, "INFO")
rec("core_vs_nc leave-s3-out: meta p", num(core_l3["meta_p"]), None, "INFO")

for c in ["edge_vs_nc", "transitory_vs_nc"]:
    cs = [r for r in rsum if r["comparison"] == c][0]
    rows = per_comp(reff, c)
    rec(f"{c}: n delta>0", sum(1 for r in rows if num(r["cliffs_delta"]) > 0), None, "INFO")
    rec(f"{c}: n delta<0", sum(1 for r in rows if num(r["cliffs_delta"]) < 0), None, "INFO")
    rec(f"{c}: sign test p", num(cs["sign_test_p"]), None, "INFO")
    rec(f"{c}: meta delta", round(num(cs["meta_delta"]), 4), None, "INFO")
    rec(f"{c}: meta p", num(cs["meta_p"]), None, "INFO")
    rec(f"{c}: meta I2", round(num(cs["meta_I2"]), 1), None, "INFO")

# --- A3. pathology effect sizes ---------------------------------------------
peff = rd(os.path.join(BASE, "effects", "pathology_effect_sizes.tsv"))
psum = rd(os.path.join(BASE, "effects", "pathology_cross_sample_summary.tsv"))
pl3  = rd(os.path.join(BASE, "effects", "pathology_leave_s3_out_summary.tsv"))
ppos = sum(1 for r in peff if num(r["cliffs_delta"]) > 0)
pneg = sum(1 for r in peff if num(r["cliffs_delta"]) < 0)
psum0 = psum[0]
rec("pathology SCC vs non-SCC: n delta>0", ppos, None, "INFO")
rec("pathology SCC vs non-SCC: n delta<0", pneg, None, "INFO")
rec("pathology: meta delta", round(num(psum0["meta_delta"]), 4), None, "INFO")
rec("pathology: meta p", num(psum0["meta_p"]), None, "INFO")
rec("pathology: meta CI lo", round(num(psum0["meta_ci_lo"]), 4), None, "INFO")
rec("pathology: meta CI hi", round(num(psum0["meta_ci_hi"]), 4), None, "INFO")
rec("pathology: meta I2", round(num(psum0["meta_I2"]), 1), None, "INFO")
rec("pathology: sign test p", num(psum0["sign_test_p"]), None, "INFO")
rec("pathology leave-s3-out: n delta>0", pl3[0]["n_delta_positive"], None, "INFO")
rec("pathology leave-s3-out: meta delta", round(num(pl3[0]["meta_delta"]), 4), None, "INFO")
rec("pathology leave-s3-out: meta p", num(pl3[0]["meta_p"]), None, "INFO")

# s9 / s10 SCC counts (the plot-NA bug check)
for s in ["s9", "s10"]:
    for r in peff:
        if r["sample"] == s:
            rec(f"{s}: n SCC spots", r["n_SCC"], None, "INFO")
            rec(f"{s}: n non-SCC spots", r["n_nonSCC"], None, "INFO")

# --- A4. robustness (UCell vs meanZ) ----------------------------------------
rspear = rd(os.path.join(BASE, "robustness", "robustness_spearman.tsv"))
rdir   = rd(os.path.join(BASE, "robustness", "robustness_direction_agreement.tsv"))
sp_neg = sum(1 for r in rspear if num(r["spearman"]) < 0)
sp_vals = [num(r["spearman"]) for r in rspear]
sp_mean = sum(sp_vals) / len(sp_vals)
dir_agree = sum(1 for r in rdir if r["agree"] == "TRUE")
dir_tot   = len(rdir)
rec("Spearman: n samples", len(rspear), None, "INFO")
rec("Spearman: n negative", sp_neg, None, "INFO")
rec("Spearman: mean", round(sp_mean, 4), None, "INFO")
rec("Spearman: min", round(min(sp_vals), 4), None, "INFO")
rec("Spearman: max", round(max(sp_vals), 4), None, "INFO")
rec("direction agreement: agree", dir_agree, None, "INFO")
rec("direction agreement: total", dir_tot, None, "INFO")
rec("direction agreement: pct", round(100*dir_agree/dir_tot, 1), None, "INFO")

# --- A5. depth-adjusted (THE CONTESTED NUMBER) ------------------------------
rdepth = rd(os.path.join(BASE, "robustness", "depth_adjusted_region_direction.tsv"))
depth_stats = {}
for c in ["core_vs_nc", "edge_vs_nc", "transitory_vs_nc"]:
    adj = [num(r["d_adj"]) for r in rdepth if r["comparison"] == c]
    raw = [num(r["d_raw"]) for r in rdepth if r["comparison"] == c]
    n = len(adj); m = sum(adj)/n
    sd = math.sqrt(sum((x-m)**2 for x in adj)/(n-1))
    se = sd/math.sqrt(n)
    t  = m/se
    # ---- TWO CANDIDATE p's, deliberately computed BOTH ways ----
    try:
        from scipy import stats as _st
        p_t    = float(_st.t.sf(abs(t), n-1) * 2)      # correct: t dist, df=11
    except Exception:
        p_t = float("nan")
    p_norm = 2*(1 - 0.5*(1 + math.erf(abs(t)/math.sqrt(2))))  # WRONG: normal approx
    depth_stats[c] = dict(n=n, mean_raw=sum(raw)/n, mean_adj=m, sd=sd, se=se, t=t,
                          p_t=p_t, p_norm=p_norm,
                          n_pos=sum(1 for x in adj if x > 0))
    rec(f"depth-adj {c}: n", n, None, "INFO")
    rec(f"depth-adj {c}: mean raw delta", round(sum(raw)/n, 4), None, "INFO")
    rec(f"depth-adj {c}: mean adjusted delta", round(m, 6), None, "INFO")
    rec(f"depth-adj {c}: n adj delta>0", sum(1 for x in adj if x > 0), None, "INFO")
    rec(f"depth-adj {c}: t statistic", round(t, 4), None, "INFO")
    rec(f"depth-adj {c}: p (t-dist df=11, CORRECT)", round(p_t, 6), None, "DEFINITION")
    rec(f"depth-adj {c}: p (normal-approx, WRONG)", round(p_norm, 6), None, "DEFINITION")

# --- A6. spot counts / detection counts -------------------------------------
spot = rd(os.path.join(BASE, "input", "spot_table_all.tsv"))
n_spot = len(spot)
ann = [r for r in spot if r["region_4class"].strip() not in ("", "NA")]
n_ann = len(ann)
from collections import Counter
by_s = Counter(r["sample"] for r in spot)
ann_by_s = Counter(r["sample"] for r in ann)
rec("total in-tissue spots", n_spot, None, "INFO")
rec("total annotated spots", n_ann, None, "INFO")
rec("per-sample in-tissue spot counts", dict(sorted(by_s.items(), key=lambda kv: int(kv[0][1:]))), None, "INFO")
rec("s3 annotated spots", ann_by_s.get("s3", 0), None, "INFO")
rec("s3 in-tissue spots", by_s.get("s3", 0), None, "INFO")

gdet = rd(os.path.join(BASE, "quality", "gene_detection_by_sample.tsv"))
n_genes_det = len(gdet)
n_lt5  = sum(1 for r in gdet if num(r["mean_frac_spots_detected"]) < 0.05)
n_lt1  = sum(1 for r in gdet if num(r["mean_frac_spots_detected"]) < 0.01)
casp5 = [r for r in gdet if r["gene"] == "CASP5"][0]
gls2  = [r for r in gdet if r["gene"] == "GLS2"][0]
rec("n genes in detection table", n_genes_det, None, "INFO")
rec("n genes mean detection <5%", n_lt5, None, "INFO")
rec("n genes mean detection <1%", n_lt1, None, "INFO")
rec("CASP5 mean detection pct", round(num(casp5["mean_frac_spots_detected"])*100, 2), None, "INFO")
rec("GLS2 mean detection pct", round(num(gls2["mean_frac_spots_detected"])*100, 2), None, "INFO")

# ============================================================================
# PART B -- EXTRACT NUMBERS FROM DELIVERED DOCUMENTS
# ============================================================================
html_txt = re.sub(r"<[^>]+>", " ", open(REPORT_HTML, encoding="utf-8").read())
html_txt = html_txt.replace("&nbsp;", " ").replace("&lt;", "<").replace("&gt;", ">") \
                   .replace("&amp;", "&").replace("&quot;", '"')
html_txt = re.sub(r"\s+", " ", html_txt)
md_txt = re.sub(r"\s+", " ", open(REPORT_MD, encoding="utf-8").read())

def find_pcts(txt):
    """collect all numeric tokens that look like a p-value"""
    pats = set()
    for m in re.finditer(r"(\d+(?:\.\d+)?)\s*[×x]\s*10\s*[⁻\-]\s*(\d+)", txt):
        pats.add(m.group(0))
    for m in re.finditer(r"p\s*=\s*([0-9]*\.?[0-9]+(?:e[+-]?\d+)?)", txt, re.I):
        pats.add(m.group(1))
    return pats

html_p = find_pcts(html_txt)
md_p   = find_pcts(md_txt)

# --- B1. the contested depth-adjusted p -------------------------------------
# Normalise both docs: strip only REAL HTML tags (markdown tables contain bare
# '<' / '>' in headers like "δ>0", which a naive <[^>]+> regex would swallow
# along with the numbers between them), map unicode superscripts, collapse ws.
SUP = {"⁻":"-","⁰":"0","¹":"1","²":"2","³":"3","⁴":"4","⁵":"5","⁶":"6","⁷":"7","⁸":"8","⁹":"9"}
HTML_TAG = re.compile(r"</?(?:div|span|p|table|tr|td|th|thead|tbody|b|i|em|strong|code|figure|figcaption|img|h[1-6]|ul|ol|li|a|br|html|head|body|style|title|meta|svg|sup|sub)\b[^>]*>", re.I)
def norm(t):
    t = HTML_TAG.sub(" ", t)
    for a,b in [("&nbsp;"," "),("&lt;","<"),("&gt;",">"),("&amp;","&"),("&quot;",'"'),("&#39;","'")]:
        t = t.replace(a,b)
    for k,v in SUP.items(): t = t.replace(k,v)
    t = t.replace("×","x").replace("−","-").replace("–","-")
    return re.sub(r"\s+"," ", t)

Hn, Mn = norm(open(REPORT_HTML, encoding="utf-8").read()), norm(open(REPORT_MD, encoding="utf-8").read())
html_txt, md_txt = Hn, Mn

# Strip out "correction / deprecation notes" before checking that a legacy wrong
# value is ABSENT -- a note explaining the fix legitimately quotes the old value.
def strip_deprecation_notes(t):
    # drop any sentence/clause that mentions 更正/audit/废弃/legacy/审计
    parts = re.split(r"(?<=[。；;])", t)
    keep = [p for p in parts if not re.search(r"更正|审计|废弃|legacy|audit|原为|曾误用|此前版本", p)]
    return " ".join(keep)

Hn_clean, Mn_clean = strip_deprecation_notes(Hn), strip_deprecation_notes(Mn)

core_d = depth_stats["core_vs_nc"]

# canonical renderings we now expect to find (post-fix)
def contains_any(t, *cs): return [c for c in cs if c in t]

correct_p_t   = f"{core_d['p_t']:.4f}"          # 0.3709
wrong_p_norm  = f"{core_d['p_norm']:.4f}"       # 0.3509
legacy_351    = "0.351"

rec("HTML depth-adj p uses correct t value (%s)" % correct_p_t,
    correct_p_t, contains_any(Hn, correct_p_t, correct_p_t[:5]), 
    "MATCH" if (correct_p_t in Hn or correct_p_t[:5] in Hn) else "MISS")
rec("HTML depth-adj p: legacy WRONG value 0.351 absent (ignoring fix-notes)",
    "absent", contains_any(Hn_clean, legacy_351), "MATCH" if legacy_351 not in Hn_clean else "MISMATCH")
rec("MD depth-adj p uses correct t value (%s)" % correct_p_t,
    correct_p_t, contains_any(Mn, correct_p_t, correct_p_t[:5]),
    "MATCH" if (correct_p_t in Mn or correct_p_t[:5] in Mn) else "MISS")
rec("MD depth-adj p: legacy WRONG value 0.351 absent (ignoring fix-notes)",
    "absent", contains_any(Mn_clean, legacy_351), "MATCH" if legacy_351 not in Mn_clean else "MISMATCH")

# --- B2. meta p for core_vs_nc (must be 3.47e-18 now) -----------------------
core_meta = num(core_sum["meta_p"])            # 3.467e-18
rec("HTML core meta p = correct 3.5 x 10^-18 (not 10^-9)",
    "3.5 x 10-18", contains_any(Hn, "3.5 x 10-18", "3.47e-18"),
    "MATCH" if ("3.5 x 10-18" in Hn or "3.47e-18" in Hn) else "MISS")
rec("HTML core meta p: legacy WRONG exponent 10^-9 absent (ignoring fix-notes)",
    "absent", contains_any(Hn_clean, "3.5 x 10-9", "3.5e-9"),
    "MATCH" if ("3.5 x 10-9" not in Hn_clean and "3.5e-9" not in Hn_clean) else "MISMATCH")
rec("MD core meta p = correct 3.47e-18",
    "3.47e-18", contains_any(Mn, "3.47e-18", "3.47 x 10-18", "3.5e-18"),
    "MATCH" if any(x in Mn for x in ["3.47e-18", "3.47 x 10-18"]) else "MISS")
rec("MD core meta p: legacy WRONG 3.5e-9 absent (ignoring fix-notes)",
    "absent", contains_any(Mn_clean, "3.5e-9", "3.5e-09"),
    "MATCH" if ("3.5e-9" not in Mn_clean and "3.5e-09" not in Mn_clean) else "MISMATCH")

# --- B3. other meta p values must still be present & correct ----------------
# MD is a "master table" that reports meta delta + meta p in the summary rows of
# sections 2/3; we require the value be present in EITHER canonical form.
for label, tv, cands in [
    ("edge_vs_nc meta p",        num([r for r in rsum if r["comparison"]=="edge_vs_nc"][0]["meta_p"]),       ["1.1 x 10-4","1.12e-4","1.121e-04","1.1e-4"]),
    ("transitory_vs_nc meta p",  num([r for r in rsum if r["comparison"]=="transitory_vs_nc"][0]["meta_p"]), ["1.9 x 10-6","1.85e-6","1.854e-06","1.9e-6"]),
    ("pathology meta p",         num(psum[0]["meta_p"]),                                                     ["7.2 x 10-6","7.22e-6","7.216e-06","7.2e-6"]),
    ("core leave-s3 meta p",     num([r for r in rl3 if r["comparison"]=="core_vs_nc"][0]["meta_p"]),        ["1.1 x 10-16","1.06e-16","1.057e-16","1.1e-16"]),
    ("pathology leave-s3 meta p",num(pl3[0]["meta_p"]),                                                      ["3.7 x 10-5","3.70e-5","3.702e-05","3.7e-5"]),
]:
    hit_h = contains_any(Hn, *cands); hit_m = contains_any(Mn, *cands)
    rec(f"HTML {label} present & correct", f"{tv:.3e}", hit_h,
        "MATCH" if hit_h else "MISS")
    rec(f"MD {label} present & correct", f"{tv:.3e}", hit_m,
        "MATCH" if hit_m else "MISS")

# --- B4. text wording per 8.docx --------------------------------------------
rec("HTML wording: '真实的空间聚集' removed (conservative wording applied)",
    "absent", contains_any(Hn, "真实的空间聚集"),
    "MATCH" if "真实的空间聚集" not in Hn else "MISMATCH")
rec("HTML wording: 'reproducible positive spatial autocorrelation' present",
    "present", contains_any(Hn, "reproducible positive spatial autocorrelation"),
    "MATCH" if "reproducible positive spatial autocorrelation" in Hn else "MISS")

# --- B5. headline spot-checks ------------------------------------------------
checks = [
    ("12/12 I>0", "12/12", lambda t: "12/12" in t),
    ("core delta range -0.465..-0.098", "-0.465", lambda t: "-0.465" in t),
    ("I range 0.057-0.230", "0.057", lambda t: "0.057" in t and "0.230" in t),
    ("meta delta -0.275", "-0.275", lambda t: "-0.275" in t),
    ("pathology meta -0.133", "-0.133", lambda t: "-0.133" in t),
    ("meanZ agreement 10.4%", "10.4", lambda t: "10.4" in t),
    ("Spearman mean -0.244", "-0.244", lambda t: "-0.244" in t),
    ("CASP5 0.49%", "0.49", lambda t: "0.49" in t),
    ("GLS2 0.60%", "0.60", lambda t: "0.60" in t),
]
for name, exp, fn in checks:
    rec(f"HTML spot-check: {name}", exp, fn(html_txt), "MATCH" if fn(html_txt) else "MISS")

# ============================================================================
# PART C -- WRITE AUDIT REPORT
# ============================================================================
def verdict():
    bad = [r for r in R if r["status"] in ("MISMATCH", "MISS")]
    return "FAIL" if bad else "PASS"
lines = []
lines.append("=" * 78)
lines.append("GSE208253 -- NUMERIC CONSISTENCY AUDIT (read-only / summary-only)")
lines.append("=" * 78)
lines.append("")
lines.append("[A] RECOMPUTED FROM RAW step1-step7 TSV  (source of truth)")
lines.append("-" * 78)
for r in R:
    if r["in_document"] is None:
        v = r["recomputed"]
        if isinstance(v, float): v = f"{v:.6g}"
        lines.append(f"  {r['item']:<52} {v}")
lines.append("")
lines.append("[B] DOCUMENT-vs-TSV COMPARISON")
lines.append("-" * 78)
for r in R:
    if r["in_document"] is not None:
        lines.append(f"  [{r['status']:<8}] {r['item']}")
        lines.append(f"             recomputed : {r['recomputed']}")
        lines.append(f"             in doc     : {r['in_document']}")
        if r["note"]: lines.append(f"             note       : {r['note']}")
lines.append("")
lines.append("=" * 78)
lines.append(f"VERDICT: {verdict()}")
lines.append("=" * 78)
lines.append("")
lines.append("KEY FINDING -- the contested depth-adjusted core-vs-nc p:")
lines.append(f"  statistic definition : one-sample t-test on the 12 per-sample")
lines.append(f"                         depth-adjusted Cliff's deltas, df = 11")
lines.append(f"  t                    : {core_d['t']:.4f}")
lines.append(f"  CORRECT p (t dist)   : {core_d['p_t']:.6f}")
lines.append(f"  WRONG   p (normal)   : {core_d['p_norm']:.6f}")
lines.append("")
lines.append("  => 0.371 == correct t-distribution p")
lines.append("  => 0.351 == normal-approximation p (INCORRECT; must be retired)")
lines.append("")

out_txt = os.path.join(OUT_DIR, "numeric_consistency_audit.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

# machine-readable
with open(os.path.join(OUT_DIR, "numeric_consistency_audit.json"), "w", encoding="utf-8") as f:
    json.dump(dict(verdict=verdict(), rows=R,
                   depth_stats=depth_stats,
                   html_pvalues=sorted(html_p), md_pvalues=sorted(md_p)),
              f, ensure_ascii=False, indent=2, default=str)

print("\n".join(lines))
print("\nwrote:", out_txt)
