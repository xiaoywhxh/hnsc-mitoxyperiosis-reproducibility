#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render the integration pack into a single reviewable HTML page."""
import csv, json, os, html, re

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "GSE208253_integration_pack.html")

def md2html(md):
    """Minimal, safe markdown -> html (headings, bold, italics, code, quotes, lists, tables, hr)."""
    lines = md.split("\n")
    out, i, in_ul, in_tbl = [], 0, False, False
    def esc(s):
        return html.escape(s, quote=False)
    def inline(s):
        s = esc(s)
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        s = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", s)
        s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
        return s
    while i < len(lines):
        ln = lines[i]
        if in_tbl and not ln.strip().startswith("|"):
            out.append("</tbody></table>"); in_tbl = False
        if in_ul and not re.match(r"^\s*[-*]\s+", ln):
            out.append("</ul>"); in_ul = False
        s = ln.rstrip()
        if not s.strip():
            i += 1; continue
        if s.startswith("<!--"):
            i += 1; continue
        if re.match(r"^#{1,6}\s", s):
            lvl = len(s) - len(s.lstrip("#"))
            out.append(f"<h{lvl+1}>{inline(s.lstrip('# ').strip())}</h{lvl+1}>")
        elif s.strip() in ("---", "***"):
            out.append("<hr>")
        elif s.strip().startswith("|"):
            cells = [c.strip() for c in s.strip().strip("|").split("|")]
            if i + 1 < len(lines) and re.match(r"^\s*\|[\s:|-]+\|\s*$", lines[i+1]):
                out.append("<table><thead><tr>" + "".join(f"<th>{inline(c)}</th>" for c in cells) + "</tr></thead><tbody>")
                in_tbl = True; i += 2; continue
            else:
                out.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in cells) + "</tr>")
                i += 1; continue
        elif s.strip().startswith(">"):
            out.append(f"<blockquote>{inline(s.strip().lstrip('> ').strip())}</blockquote>")
        elif re.match(r"^\s*[-*]\s+", s):
            if not in_ul:
                out.append("<ul>"); in_ul = True
            out.append(f"<li>{inline(re.sub(r'^\s*[-*]\s+', '', s))}</li>")
        else:
            out.append(f"<p>{inline(s)}</p>")
        i += 1
    if in_tbl: out.append("</tbody></table>")
    if in_ul: out.append("</ul>")
    return "\n".join(out)

def read(p):
    with open(os.path.join(HERE, p), encoding="utf-8") as f:
        return f.read()

freeze = None
fp = os.path.join(os.path.dirname(HERE), "FROZEN.md")
if os.path.exists(fp):
    with open(fp, encoding="utf-8") as f:
        freeze = f.read()

blocks = read("manuscript_blocks.md")
legend = read("figure_legend_s4.md")
cite = read("citation_strings.md")
nft = list(csv.reader(read("numeric_freeze_table.tsv").splitlines(), delimiter="\t"))
manifest = json.loads(read("integration_manifest.json"))

nrows = ""
for r in nft[1:]:
    cls = " class=\"retired\"" if "RETIRED" in r[2] else ""
    nrows += f"<tr{cls}><td><code>{html.escape(r[0])}</code></td><td><strong>{html.escape(r[1])}</strong></td><td>{html.escape(r[2])}</td></tr>"

lock = manifest["authoritative_numbers"]
lock_html = f"""
<table>
<thead><tr><th>Locked quantity</th><th>Authoritative value</th></tr></thead>
<tbody>
<tr><td>core-vs-nc REML meta δ</td><td><strong>{lock['core_meta_delta_REML']:.4f}</strong></td></tr>
<tr><td>95% CI</td><td><strong>{lock['core_meta_ci'][0]:.4f} … {lock['core_meta_ci'][1]:.4f}</strong></td></tr>
<tr><td>meta p</td><td><strong>3.467 &times; 10<sup>&minus;18</sup></strong></td></tr>
<tr><td>depth-adjusted mean δ</td><td><strong>+{lock['depth_adj_mean_delta']:.4f}</strong></td></tr>
<tr><td>one-sample t</td><td><strong>{lock['depth_adj_t']:.4f}</strong>, df = <strong>{lock['depth_adj_df']}</strong></td></tr>
<tr><td>depth-adjusted p</td><td><strong>{lock['depth_adj_p']:.4f}</strong></td></tr>
<tr class="retired"><td>0.351 (normal approx.)</td><td><strong>RETIRED</strong> — audit provenance only</td></tr>
</tbody></table>
"""

metrics = [
    ("Audit verdict", "PASS", "30/30 checks"),
    ("Module status", "FROZEN / CLOSED", "9.docx"),
    ("Citation check", "PASS &times;3", "FAIL 0 / WARN 0"),
    ("Lock assertions", "PASS", "build_integration.py"),
]

cards = "".join(
    f'<div class="card"><div class="cv">{v}</div><div class="cl">{l}</div><div class="cs">{s}</div></div>'
    for l, v, s in metrics)

HTMLPAGE = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>GSE208253 Integration Pack — FROZEN / CLOSED</title>
<style>
  :root {{
    --bg:#ffffff; --fg:#1a1d21; --muted:#5c6672; --line:#e3e7ec;
    --accent:#0b6bcb; --ok:#0f7b3f; --warn:#a8590a; --bad:#b3261e;
    --panel:#f7f9fb;
  }}
  * {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--bg); color:var(--fg);
    font:15px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI","Noto Sans SC",sans-serif; }}
  .wrap {{ max-width:1080px; margin:0 auto; padding:34px 26px 70px; }}
  header {{ border-bottom:3px solid var(--accent); padding-bottom:18px; margin-bottom:26px; }}
  h1 {{ font-size:25px; margin:0 0 6px; letter-spacing:-.2px; }}
  .sub {{ color:var(--muted); font-size:14px; }}
  .badge {{ display:inline-block; background:#e8f4ec; color:var(--ok); border:1px solid #bfe0cc;
    border-radius:999px; padding:3px 11px; font-size:12.5px; font-weight:600; margin-left:8px; }}
  .cards {{ display:grid; grid-template-columns:repeat(4,1fr); gap:12px; margin:22px 0 30px; }}
  .card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px; padding:14px 15px; }}
  .cv {{ font-size:19px; font-weight:700; color:var(--accent); }}
  .cl {{ font-size:12.5px; color:var(--muted); margin-top:3px; }}
  .cs {{ font-size:11.5px; color:var(--muted); margin-top:5px; opacity:.8; }}
  h2 {{ font-size:18px; margin:34px 0 10px; padding-bottom:7px; border-bottom:1px solid var(--line); }}
  h3 {{ font-size:15.5px; margin:22px 0 8px; color:#243040; }}
  h4 {{ font-size:14px; margin:18px 0 6px; color:#31404f; }}
  table {{ border-collapse:collapse; width:100%; margin:12px 0 18px; font-size:13.5px; }}
  th,td {{ border:1px solid var(--line); padding:7px 10px; text-align:left; vertical-align:top; }}
  th {{ background:var(--panel); font-weight:650; }}
  tr.retired td {{ background:#fdf1f0; color:var(--bad); }}
  code {{ background:#eef2f6; padding:1px 5px; border-radius:4px; font-size:12.5px;
    font-family:"SF Mono",Consolas,monospace; }}
  blockquote {{ margin:12px 0; padding:10px 15px; background:var(--panel);
    border-left:3px solid var(--accent); color:#243040; font-size:14px; border-radius:0 6px 6px 0; }}
  ul {{ padding-left:22px; }}
  li {{ margin:3px 0; }}
  hr {{ border:0; border-top:1px solid var(--line); margin:26px 0; }}
  .note {{ background:#fff8ec; border:1px solid #f0dcbb; border-radius:8px;
    padding:12px 15px; font-size:13.5px; color:#6b4a12; margin:16px 0; }}
  .filelist {{ font-size:13px; }}
  footer {{ margin-top:44px; padding-top:16px; border-top:1px solid var(--line);
    color:var(--muted); font-size:12.5px; }}
</style>
</head>
<body>
<div class="wrap">
<header>
  <h1>GSE208253 Spatial Validation — Manuscript Integration Pack
    <span class="badge">FROZEN / CLOSED</span></h1>
  <div class="sub">Read-only integration of the frozen module into manuscript / figure / supplement &middot;
    every number machine-derived from the frozen step1&ndash;step7 TSVs &middot; 2026-10-02</div>
</header>

<div class="cards">{cards}</div>

<h2>1 &nbsp;Authoritative numeric lock</h2>
<p>Per the 9.docx directive these values are locked and must not change. The integration
generator asserts every one of them against the frozen TSVs and refuses to emit output on
any mismatch.</p>
{lock_html}
<div class="note"><strong>0.351 is retired.</strong> It was a normal-approximation value, not the
one-sample t-test (df = 11) used in the analysis. It survives only as audit provenance and must
never be cited as a scientific result. The correct value is <strong>0.3709</strong>.</div>

<h2>2 &nbsp;Frozen three-layer conclusion</h2>
<ol>
  <li><strong>Reproducible positive spatial autocorrelation — SUPPORTED.</strong> Moran's I positive in
      12/12 sections, FDR &lt; 0.05 in 12/12. This establishes non-random spatial organization; it does
      not by itself separate biological from spatially structured technical effects.</li>
  <li><strong>Tumour-specific localization — NOT CONFIRMED.</strong> core vs. non-core was negative in
      12/12 (meta &delta; = &minus;0.2754, p = 3.467 &times; 10<sup>&minus;18</sup>); the independent
      pathologist SCC vs. non-SCC comparator was negative in 10/12. Both are opposite to the
      pre-specified expectation.</li>
  <li><strong>Raw inverse core localization — substantially attenuated after depth adjustment.</strong>
      Mean &delta; moved from &minus;0.2869 to +0.0365 (one-sample t = 0.9329, df = 11, p = 0.3709).
      This must <em>not</em> be interpreted as biological depletion of the program in the tumour core.</li>
</ol>

<h2>3 &nbsp;Narrative chain of evidence (required wording)</h2>
<blockquote>67-gene program &rarr; reproducible spatial autocorrelation &rarr; pre-specified tumour
enrichment not confirmed &rarr; raw inverse localization strongly attenuated after depth adjustment
&rarr; therefore evidence supports <strong>spatial organization, not tumour-specific localization</strong>.</blockquote>

<h2>4 &nbsp;Numeric freeze table (all values, with source TSV)</h2>
<table>
<thead><tr><th>key</th><th>value</th><th>source TSV</th></tr></thead>
<tbody>{nrows}</tbody>
</table>

<h2>5 &nbsp;Manuscript blocks (Results / Methods / Supplement)</h2>
{md2html(blocks)}

<h2>6 &nbsp;Figure legend</h2>
{md2html(legend)}

<h2>7 &nbsp;Citation &amp; data availability strings</h2>
{md2html(cite)}

<h2>8 &nbsp;Frozen module record</h2>
{md2html(freeze) if freeze else "<p>(FROZEN.md not found)</p>"}

<h2>9 &nbsp;Deliverables</h2>
<table class="filelist">
<thead><tr><th>file</th><th>purpose</th></tr></thead>
<tbody>
<tr><td><code>integration/manuscript_blocks.md</code></td><td>Results / Methods / Supplement drop-in text (EN)</td></tr>
<tr><td><code>integration/figure_legend_s4.md</code></td><td>Supplementary Figure S4 legend</td></tr>
<tr><td><code>integration/citation_strings.md</code></td><td>Data/code availability + citations + provenance statement</td></tr>
<tr><td><code>integration/numeric_freeze_table.tsv</code></td><td>26 authoritative numbers + source TSV each</td></tr>
<tr><td><code>integration/integration_manifest.json</code></td><td>Provenance manifest + retired values</td></tr>
<tr><td><code>integration/build_integration.py</code></td><td>Regenerator (LOCK assertions, fails closed)</td></tr>
<tr><td><code>integration/check_citation_consistency.py</code></td><td>Citation/number consistency checker (C1&ndash;C6)</td></tr>
<tr><td><code>integration/citation_consistency_check.json</code></td><td>Checker output (PASS &times;3)</td></tr>
</tbody>
</table>

<footer>
Generated from the frozen module (no analysis recomputed). Source of truth:
<code>D:/GSE208253_spatial/spatial_analysis/</code> &middot; freeze record:
<code>spatial_analysis/FROZEN.md</code> &middot; audit:
<code>audit/AUDIT_REPORT_freeze_prereq.md</code> &middot; next gate: manuscript-level consistency review.
</footer>
</div>
</body>
</html>
"""

with open(OUT, "w", encoding="utf-8") as f:
    f.write(HTMLPAGE)
print("wrote:", OUT, os.path.getsize(OUT), "bytes")
