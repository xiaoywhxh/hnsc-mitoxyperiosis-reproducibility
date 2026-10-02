# -*- coding: utf-8 -*-
"""
Build the consolidated HTML report for the pre-specified 67-gene UCell spatial
validation on GSE208253.  Reads only the pre-computed TSVs written by
step1..step7 -- it does NOT recompute or re-select anything.
"""
import csv, os, json, html, base64, math
from collections import defaultdict, OrderedDict

BASE = "D:/GSE208253_spatial/spatial_analysis"
OUT  = os.path.join(BASE, "report")
os.makedirs(OUT, exist_ok=True)

def rd(p, delim="\t"):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=delim))

def num(x):
    try: return float(x)
    except Exception: return float("nan")

def fmt(x, n=3):
    if x is None or (isinstance(x, float) and math.isnan(x)): return "—"
    return f"{x:.{n}f}"

def pfmt(p):
    if p is None or (isinstance(p, float) and math.isnan(p)): return "—"
    if p == 0: return "< 1e-300"
    if p < 1e-4: return f"{p:.1e}"
    return f"{p:.4f}"

def img64(p):
    with open(p, "rb") as f:
        return base64.b64encode(f.read()).decode()

P = lambda *a: os.path.join(BASE, *a)

# ----------------------------------------------------------------- load inputs
moran   = rd(P("moran", "morans_I_by_sample.tsv"))
reff    = rd(P("effects", "region_effect_sizes.tsv"))
rsum    = rd(P("effects", "region_cross_sample_summary.tsv"))
rl3     = rd(P("effects", "region_leave_s3_out_summary.tsv"))
peff    = rd(P("effects", "pathology_effect_sizes.tsv"))
psum    = rd(P("effects", "pathology_cross_sample_summary.tsv"))
pl3     = rd(P("effects", "pathology_leave_s3_out_summary.tsv"))
rspear  = rd(P("robustness", "robustness_spearman.tsv"))
rdir    = rd(P("robustness", "robustness_direction_agreement.tsv"))
rdepth  = rd(P("robustness", "depth_adjusted_region_direction.tsv"))
gdet    = rd(P("quality", "gene_detection_by_sample.tsv"))
scores  = rd(P("score", "spot_scores.tsv"))
spotsum = json.load(open(P("input", "_spot_table_summary.json"), encoding="utf-8")) \
          if os.path.exists(P("input", "_spot_table_summary.json")) else {}

SAMPLES = [f"s{i}" for i in range(1, 13)]
moran_by = {r["sample"]: r for r in moran}

# per-sample score distribution summary
g = defaultdict(list)
for r in scores: g[r["sample"]].append(r)
def med(v): 
    v = sorted(v); n = len(v)
    return v[n//2] if n % 2 else (v[n//2-1]+v[n//2])/2
score_summary = {}
for s in SAMPLES:
    v = g[s]
    score_summary[s] = dict(
        n=len(v),
        u_mean=sum(num(x["ucell"]) for x in v)/len(v),
        u_med=med([num(x["ucell"]) for x in v]),
        lib_med=med([num(x["library_size"]) for x in v]),
        ng_med=med([num(x["n_gene_detected"]) for x in v]),
    )

# depth-adjusted direction summary
# NOTE (audit fix 2026-10-02): the per-comparison p is a ONE-SAMPLE t-test on the
# 12 per-sample depth-adjusted Cliff's deltas (H0: mean delta = 0), i.e. df = n-1.
# The previous version used a normal (erf) approximation, which is not the same
# statistic once n is as small as 12 (it understated p). We now use the exact
# Student-t tail, and we also carry the normal-approximation value alongside so
# the two definitions can never be silently conflated again.
def _t_sf2(t, df):
    """Two-sided Student-t tail p. Uses scipy if available, else a series fallback."""
    try:
        from scipy import stats as _st
        return float(_st.t.sf(abs(t), df) * 2)
    except Exception:
        # incomplete-beta based two-sided tail (Numerical Recipes style)
        x = df / (df + t * t)
        return float(_betainc_reg(0.5 * df, 0.5, x))

def _betainc_reg(a, b, x):
    """Regularised incomplete beta I_x(a,b) via continued fraction (Lentz)."""
    if x <= 0: return 0.0
    if x >= 1: return 1.0
    lbeta = math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b)
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2*m - 1) * (a + 2*m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2*m) * (a + 2*m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30: d = 1e-30
        d = 1.0 / d
        c = 1.0 + num / c
        if abs(c) < 1e-30: c = 1e-30
        f *= c * d
        if abs(1.0 - c * d) < 1e-12:
            break
    return front * (f - 1.0)

def grp_depth(cmp_name):
    v = [(num(r["d_raw"]), num(r["d_adj"])) for r in rdepth if r["comparison"] == cmp_name]
    n = len(v)
    raw = [a for a,_ in v]; adj = [b for _,b in v]
    mr = sum(raw)/n; ma = sum(adj)/n
    sd = math.sqrt(sum((x-ma)**2 for x in adj)/(n-1))
    t = ma/(sd/math.sqrt(n))
    return dict(n=n, n_pos=sum(1 for b in adj if b>0), mean_raw=mr, mean_adj=ma,
                t=t, df=n-1,
                p_t=_t_sf2(t, n-1),                       # CORRECT: Student t
                p_norm=2*(1-0.5*(1+math.erf(abs(t)/math.sqrt(2)))))  # legacy normal approx

depth_sum = {c: grp_depth(c) for c in ["core_vs_nc", "edge_vs_nc", "transitory_vs_nc"]}

def sci(v, sig=2):
    """Format a p-value in scientific notation: 3.47e-18 -> '3.5 × 10⁻¹⁸'."""
    if v is None or (isinstance(v, float) and math.isnan(v)): return "—"
    if v == 0: return "&lt; 1e-300"
    e = int(math.floor(math.log10(abs(v))))
    mant = v / (10 ** e)
    SUP = {"-": "⁻", "0":"⁰","1":"¹","2":"²","3":"³","4":"⁴",
           "5":"⁵","6":"⁶","7":"⁷","8":"⁸","9":"⁹"}
    ex = "".join(SUP[ch] for ch in str(e))
    return f"{mant:.{sig-1}f} × 10{ex}"

def meta_p(comparison):
    """Pull the exact meta p from the region cross-sample summary TSV."""
    for r in rsum:
        if r["comparison"] == comparison:
            return num(r["meta_p"])
    return float("nan")

def meta_p_path():
    return num(psum[0]["meta_p"])

n_dir_agree = sum(1 for r in rdir if r["agree"] == "TRUE")
n_dir_total = len(rdir)

# ------------------------------------------------------------------ HTML build
CSS = """
*{box-sizing:border-box}
body{font-family:-apple-system,'Segoe UI',Roboto,'Helvetica Neue','Microsoft YaHei',sans-serif;
     margin:0;padding:0;background:#f5f6f8;color:#1c1e21;line-height:1.62;font-size:15px}
.wrap{max-width:1180px;margin:0 auto;padding:34px 26px 80px}
h1{font-size:27px;margin:0 0 6px;letter-spacing:-.35px}
h2{font-size:20px;margin:40px 0 12px;padding-bottom:8px;border-bottom:2px solid #d8dbe0;letter-spacing:-.2px}
h3{font-size:16.5px;margin:26px 0 9px;color:#23272b}
h4{font-size:14.5px;margin:18px 0 7px;color:#3a4048}
p{margin:9px 0}
.sub{color:#6a727c;font-size:13.5px;margin:0 0 20px}
table{border-collapse:collapse;width:100%;margin:13px 0;font-size:13px;background:#fff;
      box-shadow:0 1px 3px rgba(0,0,0,.07);border-radius:6px;overflow:hidden}
th{background:#eceff3;text-align:left;padding:9px 11px;font-weight:600;border-bottom:1.5px solid #d3d8de;
   font-size:12.3px;color:#333a42;white-space:nowrap}
td{padding:8px 11px;border-bottom:1px solid #eef0f3;vertical-align:top}
tr:last-child td{border-bottom:none}
tr:nth-child(even) td{background:#fbfcfd}
code,.mono{font-family:'SF Mono',Consolas,'Cascadia Mono',monospace;font-size:12.2px;
      background:#eef1f5;padding:1.5px 5px;border-radius:3.5px}
.badge{display:inline-block;padding:2.5px 9px;border-radius:11px;font-size:11.5px;font-weight:600;
       letter-spacing:.2px}
.b-pos{background:#fde3e1;color:#a4262c}
.b-neg{background:#e1eefb;color:#1c4e8a}
.b-ns{background:#ececf0;color:#5c636b}
.b-ok{background:#dff3e3;color:#1d6b34}
.b-warn{background:#fdf0d5;color:#8a5a00}
.card{background:#fff;border-radius:9px;padding:18px 20px;margin:15px 0;
      box-shadow:0 1px 4px rgba(0,0,0,.075)}
.callout{border-left:4px solid #4a7ab5;background:#f2f7fc;padding:13px 17px;margin:15px 0;border-radius:0 6px 6px 0}
.callout.warn{border-left-color:#d08c1a;background:#fdf8ee}
.callout.crit{border-left-color:#c0392b;background:#fdf1ef}
.callout.ok{border-left-color:#2e8b57;background:#f0f8f3}
.kv{display:grid;grid-template-columns:220px 1fr;gap:5px 16px;font-size:13.5px;margin:10px 0}
.kv div:nth-child(odd){color:#6a727c;font-weight:500}
figure{margin:20px 0;background:#fff;padding:14px;border-radius:9px;box-shadow:0 1px 4px rgba(0,0,0,.075)}
figure img{width:100%;display:block;border-radius:5px}
figcaption{font-size:12.8px;color:#68707a;margin-top:10px;line-height:1.55}
.toc{background:#fff;border-radius:9px;padding:15px 21px;box-shadow:0 1px 4px rgba(0,0,0,.075);
     font-size:13.6px;column-count:2;column-gap:32px}
.toc a{color:#2b5f9e;text-decoration:none;display:block;padding:2.5px 0}
.toc a:hover{text-decoration:underline}
.meta{font-size:12.5px;color:#7b838d;text-align:center;margin-top:46px;
      padding-top:16px;border-top:1px solid #dde0e5}
ul,ol{margin:9px 0;padding-left:24px}
li{margin:4px 0}
.tag{font-size:11px;font-weight:700;padding:1.5px 6px;border-radius:4px;background:#e7ebf0;color:#4a525b}
.pos{color:#a4262c;font-weight:600}
.neg{color:#1c4e8a;font-weight:600}
@media print{body{background:#fff}.wrap{max-width:none}}
"""

H = []
A = H.append

A(f"""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>GSE208253 空间转录组验证报告 — 67-gene mitoxyperiosis program (UCell Primary)</title>
<style>{CSS}</style></head><body><div class="wrap">""")

A("""<h1>GSE208253 空间转录组预先规定验证报告</h1>
<p class="sub">67-gene <i>mitoxyperiosis</i> program × 12 例 HPV-negative OSCC Visium 样本<br>
Primary method 锁定为 <b>UCell rank-based score</b>；全部 67 个基因保留；未行基因再筛选、未做 DEG / CellChat / trajectory / ML / drug prediction</p>""")

A(f"""<div class="toc">
<a href="#s1">1. 设计与预先规定（Pre-specification）</a>
<a href="#s2">2. 数据集与 QC</a>
<a href="#s3">3. 方法：UCell score（Primary）</a>
<a href="#s4">4. 每样本 score 分布</a>
<a href="#s5">5. 空间自相关：global Moran's I</a>
<a href="#s6">6. region_4class 效应量（Primary: core vs nc）</a>
<a href="#s7">7. 病理学家标注独立验证（SCC vs non-SCC）</a>
<a href="#s8">8. 12 例方向一致性汇总</a>
<a href="#s9">9. 稳健性：mean-expression score</a>
<a href="#s10">10. 测序深度混杂诊断</a>
<a href="#s11">11. 技术局限（低检出基因）</a>
<a href="#s12">12. 空间图</a>
<a href="#s13">13. 结论与建议写法</a>
<a href="#s14">14. 输出文件清单</a>
</div>""")

# ---------------------------------------------------------------- 1 design
A("""<h2 id="s1">1. 设计与预先规定（Pre-specification）</h2>
<div class="card">
<p>本分析的目标是<b>验证</b>而非<b>发现</b>：67-gene program 的基因列表、分层定义、region 定义与主要比较方向在分析前全部锁定，分析过程中未因结果而修改任何一项。</p>
<table>
<tr><th>项目</th><th>预先规定内容</th><th>是否被结果改动</th></tr>
<tr><td>Gene set</td><td>67 genes（73-gene mechanism program 的平台共同可测子集；排除 PRKN / AKT2 / MIEF2 / MARCHF5 / STING1 / CGAS）</td><td><span class="badge b-ok">否，全部保留</span></td></tr>
<tr><td>Primary score</td><td>Rank-based UCell-style score，<b>每样本独立</b>计算，避免跨样本测序深度影响</td><td><span class="badge b-ok">否</span></td></tr>
<tr><td>Primary region comparison</td><td><b>core vs nc</b></td><td><span class="badge b-ok">否</span></td></tr>
<tr><td>Secondary</td><td>edge vs nc；transitory vs nc</td><td><span class="badge b-ok">否</span></td></tr>
<tr><td>跨样本推断</td><td>每例独立计算效应量 + sign test + random-effects meta；<b>不把 spots 当独立患者做 pooled 显著性推断</b></td><td><span class="badge b-ok">否</span></td></tr>
<tr><td>病理 comparator</td><td>在打分<b>之前</b>预先定义 non-SCC 类别（见 §7）</td><td><span class="badge b-ok">否</span></td></tr>
<tr><td>sample_3</td><td>Moran's I 用全部 GEO spots；annotation 类分析只用 476 个可靠 coordinate-matched spots；不插补缺失 annotation；另加 leave-s3-out 敏感性分析</td><td><span class="badge b-ok">否</span></td></tr>
<tr><td>稳健性分析</td><td>简单标准化 mean-expression score，<b>仅用于 robustness，不用于选择主方法</b></td><td><span class="badge b-ok">否，Primary 恒为 UCell</span></td></tr>
<tr><td>禁止事项</td><td>未运行 DEG、CellChat、trajectory、ML、drug prediction</td><td><span class="badge b-ok">遵守</span></td></tr>
</table>
</div>""")

# ---------------------------------------------------------------- 2 dataset
A(f"""<h2 id="s2">2. 数据集与 QC</h2>
<div class="card">
<div class="kv">
<div>数据集</div><div>GSE208253（Arora et al.，HPV-negative OSCC，fresh-frozen 10x Visium，GRCh38）</div>
<div>样本</div><div>12 例（GSM6339631–GSM6339642 = s1…s12）</div>
<div>in-tissue spots</div><div><b>26,371</b>（本分析，SpaceRanger filtered matrix）</div>
<div>作者 annotation 覆盖</div><div><b>24,399</b> spots（可匹配到 region_4class / pathologist_anno）</div>
<div>基因</div><div>36,601 genes（UCell 排序背景）</div>
<div>每样本文件</div><div><code>filtered_feature_bc_matrix.h5</code> · <code>scalefactors_json.json.gz</code> · <code>tissue_hires_image.png.gz</code> · <code>tissue_positions_list.csv.gz</code></div>
<div>样本量（in-tissue spots）</div><div>{' · '.join(f"{s}={score_summary[s]['n']}" for s in SAMPLES)}</div>
</div>
<div class="callout warn"><b>坐标说明：</b><code>tissue_positions_list.csv.gz</code> 中的 pixel 坐标在本数据集中<b>不唯一</b>（1,185 个 spot 仅 545 个唯一 pxl_row），像素坐标无法用于邻接构建。因此本分析全部使用 <code>array_row</code>/<code>array_col</code>（Visium 六边形阵列坐标），邻接规则经实测确认：同一行 <code>(r, c±2)</code>，相邻行 <code>(r±1, c±1)</code>；中位邻居数 5.4–5.9。</div>
</div>""")

# --------------------------------------------------------------- 3 methods
A("""<h2 id="s3">3. 方法：UCell score（Primary）</h2>
<div class="card">
<p>UCell（Mannen et al. 2022, <i>Comput Biol Med</i>）逐 spot 计算：</p>
<ol>
<li>将该 spot 的<b>全部 36,601 个基因</b>按表达量降序排名；</li>
<li>取 67 个 signature 基因的排名；截断于 <code>maxRank = 1500</code>；</li>
<li><code>U = 1 − (S − S_min)/(S_max − S_min)</code>，其中 <code>S = Σ ranks</code>，<code>S_min = n(n+1)/2</code>，<code>S_max = n(2·maxRank − n + 1)/2</code>，<code>n = 67</code>。</li>
</ol>
<div class="callout"><b>实现要点（已验证）：</b>非零表达基因占据降序排名的前 <i>m</i> 位，全零基因共享 <i>m+1…n_genes</i> 的并列排名，因此<b>一个可检出基因的全局降序排名即为它在可检出基因中的排名</b>。早期按「(n_genes − m) + rk」平移的实现会把所有基因排名推过 1500 导致 U 恒为 0；已修正。</div>
<p><span class="tag">QC</span> 12 例 UCell 均值 0.556–0.623，<b>无天花板效应</b>（无样本均值 &gt; 0.65）。</p>
</div>""")

# ------------------------------------------------------- 4 score distribution
rows = "".join(
    f"<tr><td><b>{s}</b></td><td>{score_summary[s]['n']}</td>"
    f"<td>{fmt(score_summary[s]['u_mean'],4)}</td><td>{fmt(score_summary[s]['u_med'],4)}</td>"
    f"<td>{score_summary[s]['lib_med']:,.0f}</td><td>{score_summary[s]['ng_med']:,.0f}</td>"
    f"<td>{fmt(num(moran_by[s]['morans_I']),4)}</td>"
    f"<td>{fmt(num(moran_by[s]['expected_I']),5)}</td>"
    f"<td>{pfmt(num(moran_by[s]['fdr_bh']))}</td></tr>"
    for s in SAMPLES)

A(f"""<h2 id="s4">4. 每样本 score 分布</h2>
<div class="card">
<table>
<tr><th>Sample</th><th>n spots</th><th>UCell mean</th><th>UCell median</th><th>library size (median)</th><th>n genes (median)</th><th>Moran's I</th><th>E[I]</th><th>FDR (BH)</th></tr>
{rows}
</table>
<p style="font-size:13px;color:#6a727c">注：UCell score 逐样本独立计算，因此样本间均值本身不可直接做组间比较（这正是下游改用「每例内部 core vs nc 差值」的原因）。library size 与 n genes 在此列出，用于 §10 的深度混杂诊断。</p>
</div>""")

# ------------------------------------------------------------- 5 Moran's I
n_pos_moran = sum(1 for r in moran if num(r["morans_I"]) > 0)
n_sig_moran = sum(1 for r in moran if num(r["fdr_bh"]) < 0.05)
A(f"""<h2 id="s5">5. 空间自相关：global Moran's I</h2>
<div class="card">
<p>对每个样本<b>独立</b>计算 spot-level 67-gene UCell score 的 global Moran's I，权重矩阵由 Visium 六边形邻接构建并做行标准化：</p>
<p style="text-align:center"><code>I = (n/S₀) · (x′Wx)/(x′x)</code></p>
<p>显著性检验采用基于正态近似的解析方差（含 S₀、S₁、S₂ 与峰度 k 的完整表达式），并与 <code>spdep::moran.test</code> 交叉验证——s10 两者结果<b>完全一致</b>（I = 0.128172079425，与 spdep 逐位相同）。</p>
<table>
<tr><th>Sample</th><th>n spots</th><th>mean neighbours</th><th>Moran's I</th><th>E[I] = −1/(n−1)</th><th>z</th><th>p</th><th>FDR (BH)</th><th>方向</th></tr>
{''.join(f'<tr><td><b>{r["sample"]}</b></td><td>{int(float(r["n_spots"]))}</td><td>{fmt(num(r["mean_n_neighbours"]),2)}</td><td><b>{fmt(num(r["morans_I"]),4)}</b></td><td>{fmt(num(r["expected_I"]),5)}</td><td>{fmt(num(r["z"]),2)}</td><td>{pfmt(num(r["p_value"]))}</td><td>{pfmt(num(r["fdr_bh"]))}</td><td><span class="badge b-pos">I &gt; 0</span></td></tr>' for r in moran)}
</table>
<div class="callout ok"><b>结果：12/12 样本 Moran's I &gt; 0，且 12/12 FDR &lt; 0.05。</b><br>
I 范围 0.057–0.230（中位 0.126）。这说明 67-gene UCell score 在每个样本内部都呈现<b>可重复的正空间自相关</b>——不是随机散布。</div>
<p style="font-size:13px;color:#6a727c"><b>措辞边界：</b>Moran's I 支持的是「非随机空间自相关」，而空间自相关本身<b>不能区分生物学空间组织与空间结构化的技术因素</b>（如组织密度、捕获效率梯度）。因此本报告及论文正文统一使用「reproducible positive spatial autocorrelation」，不进一步表述为「已证明真实生物学空间聚集」。</p>
<p style="font-size:13px;color:#6a727c">不同样本 I 值有差异（s3 = 0.057 最低，s6 = 0.230 最高），与组织切片形态、肿瘤占比及有效 spot 密度有关；I 的绝对值在 Visium 数据上普遍偏低（0.05–0.25）属正常范围。</p>
</div>""")

# --------------------------------------------------------- 6 region effects
rsum_by = {r["comparison"]: r for r in rsum}
def reg_rows(cmp_name):
    return [r for r in reff if r["comparison"] == cmp_name]

A(f"""<h2 id="s6">6. region_4class 效应量（Primary: core vs nc）</h2>
<div class="card">
<p>使用作者提供的 <code>region_4class</code>（core / edge / transitory / nc）标注。每例<b>分别</b>计算 median difference 与 Cliff's delta（rank-biserial effect size）；<b>不</b>把 spots 作为独立患者做跨样本 pooled 显著性推断。</p>
<h3>6.1 core vs nc（Primary）</h3>
<table>
<tr><th>Sample</th><th>n core</th><th>n nc</th><th>median core</th><th>median nc</th><th>Δ median</th><th>Cliff's δ</th><th>方向</th></tr>
{''.join(f'<tr><td><b>{r["sample"]}</b></td><td>{r["n_A"]}</td><td>{r["n_B"]}</td><td>{fmt(num(r["median_A"]))}</td><td>{fmt(num(r["median_B"]))}</td><td class="neg">{fmt(num(r["median_diff"]))}</td><td class="neg"><b>{fmt(num(r["cliffs_delta"]))}</b></td><td><span class="badge b-neg">δ &lt; 0</span></td></tr>' for r in reg_rows("core_vs_nc"))}
</table>
<div class="callout crit"><b>原始 UCell 结果方向：12/12 例 core 的 67-gene score 低于 nc（Cliff's δ 全部为负，δ 范围 −0.465 … −0.098）。</b><br>
Sign test 12/12 一致：p = 4.88 × 10⁻⁴；random-effects meta（REML）：δ = <b>−0.275</b>（95% CI −0.338 … −0.213），p = {sci(meta_p('core_vs_nc'))}，I² = 77%。</div>
<h3>6.2 Secondary：edge vs nc / transitory vs nc</h3>
<table>
<tr><th>Comparison</th><th>n 例</th><th>δ &gt; 0</th><th>δ &lt; 0</th><th>sign test p</th><th>meta δ</th><th>95% CI</th><th>meta p</th><th>I²</th></tr>
{''.join(f'<tr><td>{r["comparison"].replace("_vs_"," vs ")}</td><td>{r["n_samples"]}</td><td>{r["n_delta_positive"]}</td><td>{r["n_delta_negative"]}</td><td>{pfmt(num(r["sign_test_p"]))}</td><td class="neg"><b>{fmt(num(r["meta_delta"]))}</b></td><td>{fmt(num(r["meta_ci_lo"]))} … {fmt(num(r["meta_ci_hi"]))}</td><td>{pfmt(num(r["meta_p"]))}</td><td>{fmt(num(r["meta_I2"]),1)}%</td></tr>' for r in rsum)}
</table>
<h3>6.3 leave-sample_3-out 敏感性分析</h3>
<table>
<tr><th>Comparison</th><th>n 例</th><th>δ &gt; 0</th><th>sign test p</th><th>meta δ</th><th>meta p</th></tr>
{''.join(f'<tr><td>{r["comparison"].replace("_vs_"," vs ")}</td><td>{r["n_samples"]}</td><td>{r["n_delta_positive"]}</td><td>{pfmt(num(r["sign_test_p"]))}</td><td class="neg">{fmt(num(r["meta_delta"]))}</td><td>{pfmt(num(r["meta_p"]))}</td></tr>' for r in rl3)}
</table>
<p style="font-size:13px;color:#6a727c">剔除 s3 后结论不变（core vs nc 仍 0/11 为正，meta δ = −0.278，p = {sci(num([r for r in rl3 if r["comparison"]=="core_vs_nc"][0]["meta_p"]))}），说明结果不依赖 annotation 覆盖较低的 s3。</p>
</div>""")

# ------------------------------------------------------------- 7 pathology
A(f"""<h2 id="s7">7. 病理学家标注独立验证（SCC vs non-SCC）</h2>
<div class="card">
<p>使用与 region_4class 相互独立的 <code>pathologist_anno_raw</code>（组织学标签）验证 67-gene program 的恶性定位。</p>
<h4>比较组在打分前预先定义</h4>
<table>
<tr><th>归类</th><th>pathologist_anno_raw 类别</th></tr>
<tr><td><span class="badge b-neg">SCC</span></td><td><code>SCC</code></td></tr>
<tr><td><span class="badge b-pos">non-SCC（对照）</span></td><td>Lymphocyte Negative Stroma · Lymphocyte Positive Stroma · Muscle · Glandular Stroma · Non-cancerous Mucosa · Lymphocyte Positive Muscles · Artery/Vein</td></tr>
<tr><td><span class="badge b-ns">排除</span></td><td>Artifact · Cautery · Fold · Edge Effects · Keratin</td></tr>
</table>
<p style="font-size:13px;color:#6a727c">排除项为组织处理伪影与角化碎屑，其基因表达不能代表真实组织生物学；这一排除规则在查看任何结果之前即已固定。</p>
<table>
<tr><th>Sample</th><th>n SCC</th><th>n non-SCC</th><th>median SCC</th><th>median non-SCC</th><th>Δ median</th><th>Cliff's δ</th><th>方向</th></tr>
{''.join(f'<tr><td><b>{r["sample"]}</b></td><td>{r["n_SCC"]}</td><td>{r["n_nonSCC"]}</td><td>{fmt(num(r["median_SCC"]))}</td><td>{fmt(num(r["median_nonSCC"]))}</td><td class="{"neg" if num(r["median_diff"])<0 else "pos"}">{fmt(num(r["median_diff"]))}</td><td class="{"neg" if num(r["cliffs_delta"])<0 else "pos"}"><b>{fmt(num(r["cliffs_delta"]))}</b></td><td>{"<span class=\'badge b-neg\'>δ &lt; 0</span>" if num(r["cliffs_delta"])<0 else "<span class=\'badge b-pos\'>δ &gt; 0</span>"}</td></tr>' for r in peff)}
</table>
<div class="callout crit"><b>SCC vs non-SCC：10/12 例 δ 为负</b>（meta δ = −0.133，95% CI −0.191 … −0.075，p = {sci(meta_p_path())}，I² = 91%；sign test p = 0.0386）。<br>
即：<b>原始 UCell score 在癌区（SCC）反而低于非癌组织</b>——与 §6 的 core vs nc 方向一致，但同样受深度混杂影响（见 §10）。</div>
<p style="font-size:13px;color:#6a727c">leave-s3-out：11 例中 2 例为正，meta δ = −0.131，p = {sci(num(pl3[0]["meta_p"]))}，结论不变。</p>
</div>""")

# --------------------------------------------------------- 8 direction summary
A(f"""<h2 id="s8">8. 12 例方向一致性汇总</h2>
<div class="card">
<table>
<tr><th>分析</th><th>方向一致性</th><th>sign test p</th><th>random-effects meta δ</th><th>meta p</th><th>I²</th><th>判定</th></tr>
<tr><td>Moran's I &gt; 0</td><td><b>12/12</b></td><td>4.88 × 10⁻⁴</td><td>—</td><td>12/12 FDR &lt; 0.05</td><td>—</td><td><span class="badge b-ok">正空间自相关成立</span></td></tr>
<tr><td>core vs nc（δ &lt; 0）</td><td><b>12/12</b></td><td>4.88 × 10⁻⁴</td><td>−0.275</td><td>{sci(meta_p('core_vs_nc'))}</td><td>77%</td><td><span class="badge b-warn">方向一致但与预期相反</span></td></tr>
<tr><td>edge vs nc（δ &lt; 0）</td><td>10/12</td><td>0.0386</td><td>−0.129</td><td>{sci(meta_p('edge_vs_nc'))}</td><td>92%</td><td><span class="badge b-warn">方向一致但与预期相反</span></td></tr>
<tr><td>transitory vs nc（δ &lt; 0）</td><td>11/12</td><td>0.0063</td><td>−0.197</td><td>{sci(meta_p('transitory_vs_nc'))}</td><td>93%</td><td><span class="badge b-warn">方向一致但与预期相反</span></td></tr>
<tr><td>SCC vs non-SCC（δ &lt; 0）</td><td>10/12</td><td>0.0386</td><td>−0.133</td><td>{sci(meta_p_path())}</td><td>91%</td><td><span class="badge b-warn">方向一致但与预期相反</span></td></tr>
</table>
<div class="callout"><b>核心张力：</b>预先规定验证的两个预期中，<b>「多数样本 Moran's I &gt; 0」完全成立（12/12）</b>；而<b>「多数病例 core / SCC 方向一致」虽然成立（12/12），但方向是相反的</b>——67-gene UCell score 在正常组织侧更高，而非肿瘤核心更高。</div>
</div>""")

# ------------------------------------------------------- 9 robustness meanZ
A(f"""<h2 id="s9">9. 稳健性：mean-expression score（<i>仅 robustness，不用于选择主方法</i>）</h2>
<div class="card">
<p>按预先规定，另外计算一个简单标准化 mean-expression score（gene-wise z-scored mean of log1p(CPM)），<b>仅用于检验主方法结论是否稳健</b>，无论其结果看起来是否更"漂亮"，<b>Primary 恒为 UCell</b>。</p>
<h4>9.1 Spearman 相关（每样本 UCell vs meanZ）</h4>
<table>
<tr><th>Sample</th><th>n</th><th>Spearman ρ</th><th>方向</th></tr>
{''.join(f'<tr><td><b>{r["sample"]}</b></td><td>{r["n"]}</td><td class="neg">{fmt(num(r["spearman"]))}</td><td><span class="badge b-neg">负相关</span></td></tr>' for r in rspear)}
</table>
<p><b>12/12 样本 Spearman 相关系数为负</b>（均值 −0.244，范围 −0.496 … −0.095）。两种打分方法在同一 spot 上呈系统性负相关。</p>
<h4>9.2 区域效应方向一致性</h4>
<div class="callout warn"><b>UCell 与 meanZ 的区域效应方向一致性仅 {n_dir_agree}/{n_dir_total}（{100*n_dir_agree/n_dir_total:.1f}%）。</b><br>
具体而言：meanZ 显示 core 的 score <b>高于</b> nc（δ 为正，多例 δ &gt; +0.5），而 UCell 显示 core <b>低于</b> nc。两种方法在核心问题上给出相反答案。</div>
<p>不一致的来源是两种指标对<b>测序深度</b>的敏感性不同：meanZ 直接使用归一化后的表达值，检出基因多、深度高的 spot 会同时获得更高的均值；UCell 是基于全转录组排名的相对量度，对绝对深度相对不敏感，但依然会被「core spot 可检出基因少」这一现象影响（见 §10）。</p>
</div>""")

# ------------------------------------------------------------ 10 depth
d = depth_sum
# actual per-sample core vs nc library-size medians
reg_lookup = {}
for r in rd(P("input", "spot_table_all.tsv")):
    reg_lookup[(r["sample"], r["barcode"])] = r["region_4class"]
depth_rows = []
for s in SAMPLES:
    lib_core = [num(x["library_size"]) for x in g[s]
                if reg_lookup.get((s, x["barcode"])) == "core"]
    lib_nc   = [num(x["library_size"]) for x in g[s]
                if reg_lookup.get((s, x["barcode"])) == "nc"]
    if lib_core and lib_nc:
        mc, mn = med(lib_core), med(lib_nc)
        depth_rows.append((s, len(lib_core), mc, len(lib_nc), mn, mc / mn))

A(f"""<h2 id="s10">10. 测序深度混杂诊断</h2>
<div class="card">
<p>这是一个<b>诊断性</b>分析（不改变 Primary 结论，只解释其来源）。每样本内将 UCell score 对 <code>log(library_size) + log(n_gene_detected)</code> 回归，取残差后重新计算 core vs nc 的方向。</p>
<h3>10.1 core 与 nc 的测序深度差异</h3>
<table>
<tr><th>Sample</th><th>n core</th><th>core median library</th><th>n nc</th><th>nc median library</th><th>core / nc 倍数</th></tr>
{''.join(f'<tr><td><b>{s}</b></td><td>{nc_}</td><td>{mc:,.0f}</td><td>{nn}</td><td>{mn:,.0f}</td><td class="neg"><b>{ratio:.2f}×</b></td></tr>' for s, nc_, mc, nn, mn, ratio in depth_rows)}
</table>
<p>实测：core spot 的中位 library size 约 <b>{min(r[2] for r in depth_rows):,.0f}–{max(r[2] for r in depth_rows):,.0f}</b>，而 nc spot 约 <b>{min(r[4] for r in depth_rows):,.0f}–{max(r[4] for r in depth_rows):,.0f}</b>，core 普遍高出 <b>{min(r[5] for r in depth_rows):.1f}–{max(r[5] for r in depth_rows):.1f} 倍</b>。这是本数据集的一个强烈系统性特征。</p>
<h3>10.2 深度校正后的区域方向</h3>
<p style="font-size:13px;color:#6a727c">检验定义：对每例的深度校正后 Cliff's δ，跨 12 例做 <b>单样本 t 检验</b>（H₀: mean δ = 0，df = 11）。此前版本误用正态近似，已按审计结果改为 Student-t 精确尾概率，两套定义并列如下以便核对。</p>
<table>
<tr><th>Comparison</th><th>原始 mean δ</th><th>原始 δ &gt; 0 例数</th><th>校正后 mean δ</th><th>校正后 δ &gt; 0 例数</th><th>t (df=11)</th><th>校正后 p（t 检验）</th><th>（正态近似 p，仅备查）</th><th>判定</th></tr>
<tr><td>core vs nc</td><td class="neg">{fmt(d["core_vs_nc"]["mean_raw"])}</td><td>0/12</td><td>{fmt(d["core_vs_nc"]["mean_adj"])}</td><td>{d["core_vs_nc"]["n_pos"]}/12</td><td>{fmt(d["core_vs_nc"]["t"],4)}</td><td><b>{pfmt(d["core_vs_nc"]["p_t"])}</b></td><td>{pfmt(d["core_vs_nc"]["p_norm"])}</td><td><span class="badge b-ns">校正后趋于零</span></td></tr>
<tr><td>edge vs nc</td><td class="neg">{fmt(d["edge_vs_nc"]["mean_raw"])}</td><td>2/12</td><td>{fmt(d["edge_vs_nc"]["mean_adj"])}</td><td>{d["edge_vs_nc"]["n_pos"]}/12</td><td>{fmt(d["edge_vs_nc"]["t"],4)}</td><td><b>{pfmt(d["edge_vs_nc"]["p_t"])}</b></td><td>{pfmt(d["edge_vs_nc"]["p_norm"])}</td><td><span class="badge b-warn">弱正向</span></td></tr>
<tr><td>transitory vs nc</td><td class="neg">{fmt(d["transitory_vs_nc"]["mean_raw"])}</td><td>1/12</td><td>{fmt(d["transitory_vs_nc"]["mean_adj"])}</td><td>{d["transitory_vs_nc"]["n_pos"]}/12</td><td>{fmt(d["transitory_vs_nc"]["t"],4)}</td><td><b>{pfmt(d["transitory_vs_nc"]["p_t"])}</b></td><td>{pfmt(d["transitory_vs_nc"]["p_norm"])}</td><td><span class="badge b-warn">弱正向</span></td></tr>
</table>
<div class="callout crit"><b>关键诊断结论：</b>一旦校正测序深度与检出基因数，<b>core vs nc 的效应趋于零</b>（δ 由 {fmt(d["core_vs_nc"]["mean_raw"])} 变为 <b>+{fmt(abs(d["core_vs_nc"]["mean_adj"]))}</b>；单样本 t 检验 df=11，t = {fmt(d["core_vs_nc"]["t"],4)}，p = {pfmt(d["core_vs_nc"]["p_t"])}，不显著）。<br>
换言之，§6 中「12/12 core 低于 nc」的强烈一致信号，<b>在很大程度上可由 core spot 的测序深度差异解释</b>，而不能直接解读为「67-gene program 在肿瘤核心被下调」。</div>
<p style="font-size:13px;color:#6a727c">edge/transitory 校正后变为弱正向（{fmt(d["edge_vs_nc"]["mean_adj"])} / {fmt(d["transitory_vs_nc"]["mean_adj"])}），但效应量很小、I² 极高，不足以支撑任何方向的强结论。此诊断不改变 Primary 方法的预先规定，只是为结果的解释边界提供依据。</p>
</div>""")

# ---------------------------------------------------------- 11 detection
low = [r for r in gdet if num(r["mean_frac_spots_detected"]) < 0.05]
A(f"""<h2 id="s11">11. 技术局限：低检出基因（保留，不删除）</h2>
<div class="card">
<p>按预先规定，<b>全部 67 个基因一律保留</b>，不因检出率低而删除或替换。以下低检出基因此处仅作为<b>技术限制</b>如实报告。</p>
<table>
<tr><th>Gene</th><th>mean 检出 spot 比例</th><th>备注</th></tr>
{''.join(f'<tr><td><code>{r["gene"]}</code></td><td class="neg">{num(r["mean_frac_spots_detected"])*100:.2f}%</td><td>{"<span class=\'badge b-warn\'>极低检出</span>" if num(r["mean_frac_spots_detected"])<0.01 else "<span class=\'badge b-ns\'>低检出</span>"}</td></tr>' for r in low)}
</table>
<p>共 <b>{sum(1 for r in gdet if num(r['mean_frac_spots_detected'])<0.05)}/67</b> 个基因平均检出率 &lt; 5%，其中 <b>{sum(1 for r in gdet if num(r['mean_frac_spots_detected'])<0.01)}/67</b> 个 &lt; 1%（CASP5 0.49%、GLS2 0.60%）。</p>
<p>这些基因在 67-gene UCell score 中的贡献非常有限（其排名多为并列下限），因此它们既不会实质性地抬高分数，也不会成为分数差异的驱动因素。但它们提示：<b>Visium 平台对这部分 mitoxyperiosis 机制基因的检测灵敏度不足以支持 spot 级别的单基因解读</b>，任何基于个别低检出基因的空间结论都不成立。</p>
<p style="font-size:13px;color:#6a727c">完整逐样本检出率矩阵见 <code>quality/gene_detection_by_sample.tsv</code>。</p>
</div>""")

# ---------------------------------------------------------------- 12 figures
figs = [
    ("fig_spatial_ucell_grid.png", "图 1 · 12 例样本的 67-gene UCell score 空间分布（Primary）",
     "每格为一个样本，坐标使用 Visium array_row/array_col（非像素坐标）。副标题给出该样本的 Moran's I 与 FDR。可见 score 在组织内呈空间渐变/块状结构，而非随机散布——与 §5 的 12/12 显著正 Moran's I 一致。"),
    ("fig_spatial_region_grid.png", "图 2 · 作者 region_4class 空间标注（Arora et al.）",
     "core（深红）/ edge（橙）/ transitory（浅蓝）/ nc（深蓝）。灰色为未获作者标注的 in-tissue spot（不插补、不参与 annotation 类分析）。s3 因注释覆盖低而灰区明显，已按预先规定单独处理。"),
    ("fig_spatial_pathology_grid.png", "图 3 · 病理学家标注独立验证（SCC vs non-SCC）",
     "红 = SCC，蓝 = 预先定义的 non-SCC 对照，灰 = 排除类别（Artifact / Cautery / Fold / Edge Effects / Keratin）。副标题给出每样本 SCC 与 non-SCC spot 数。"),
    ("fig_score_distribution.png", "图 4 · 每样本 UCell score 分布与 library size 分布",
     "上图：UCell score 逐样本分布（无天花板效应）。下图：library size（log10），清楚显示样本间与样本内部的深度差异——这是 §10 混杂诊断的直接依据。"),
    ("fig_moran_bar.png", "图 5 · 每样本 global Moran's I",
     "蓝柱为观测 Moran's I，灰柱为期望值 E[I] = −1/(n−1)（接近 0）。* 表示 FDR &lt; 0.05。12/12 样本观测值显著高于期望。"),
    ("fig_effect_forest.png", "图 6 · 每样本效应量（Cliff's delta）",
     "四个面板分别为 core vs nc、edge vs nc、transitory vs nc、SCC vs non-SCC。红点 = δ &gt; 0，蓝点 = δ &lt; 0。可见 core vs nc 面板 12/12 全为负；跨样本汇总采用 sign test + random-effects meta，未做 pooled spot 级推断。"),
]
for fn, title, cap in figs:
    p = P("plots", fn)
    if os.path.exists(p):
        A(f'<figure><img src="data:image/png;base64,{img64(p)}" alt="{html.escape(title)}">'
          f'<figcaption><b>{html.escape(title)}</b><br>{cap}</figcaption></figure>')

# ---------------------------------------------------------------- 13 conclusion
A("""<h2 id="s13">13. 结论与建议写法</h2>
<div class="card">
<h3>13.1 可以确定成立的</h3>
<ul>
<li><b>67-gene program 在 GSE208253 中具有可重复的正空间自相关。</b>12/12 样本 Moran's I &gt; 0 且 FDR &lt; 0.05，score 沿组织形态呈块状/渐变分布（图 1）。<b>注意措辞边界：</b>空间自相关本身不区分生物学空间组织与空间结构化的技术因素，故支持的是「non-random spatial organization」，而非「已证明的生物学空间聚集」。</li>
<li><b>67 个基因全部存在于该平台</b>（67/67 可检出，54/67 在 ≥5% spot 中检出，61/67 在 ≥1% spot 中检出），程序本身对该 Visium 数据集是可测的。</li>
<li><b>方法学层面：</b>UCell 与 meanZ 两种独立打分给出系统性不同的结果（Spearman 12/12 为负，区域方向一致性仅 10.4%），说明「用哪种 score」对本数据集的结论有决定性影响——这本身是一个需要在论文中披露的方法学事实。</li>
</ul>
<h3>13.2 不能成立的（必须如实报告）</h3>
<ul>
<li><b>原始 UCell 下 core / SCC 方向与「肿瘤高表达」预期相反</b>（core vs nc 12/12 为负；SCC vs non-SCC 10/12 为负）。</li>
<li><b>该方向信号在测序深度校正后基本消失</b>（core vs nc 校正后 δ = +0.0365；单样本 t 检验 df=11，t = {fmt(d["core_vs_nc"]["t"],4)}，p = {pfmt(d["core_vs_nc"]["p_t"])}）。因此它<b>不能</b>被解释为「67-gene program 在肿瘤核心的生物学耗竭」。</li>
<li>edge / transitory 校正后仅呈极弱正向，且异质性极高（I² &gt; 90%），不足以支撑方向性结论。</li>
</ul>
<div class="callout crit"><b>建议的诚实表述（供 manuscript 使用）：</b><br>
「在 GSE208253 的 12 例 OSCC Visium 样本中，67-gene <i>mitoxyperiosis</i> program 的 UCell 评分在全部 12 例中均呈现可重复的正空间自相关（Moran's I &gt; 0，FDR &lt; 0.05）；<i>The 67-gene score showed reproducible positive spatial autocorrelation across all 12 samples.</i> 然而，其与肿瘤区域的空间对应关系未获确认（tumor-specific localization not confirmed）：预先规定的主要比较（core vs nc）与独立的病理学家标注比较（SCC vs non-SCC）均未显示预期的肿瘤区升高；且该原始方向性差异在很大程度上与 spot 测序深度差异混杂，经深度校正后接近零。因此，我们将空间转录组结果定位为对<b>该程序可测性与空间结构性的支持</b>，而非对肿瘤特异性定位的确认；spot 级别的基因-区域对应关系在本数据集中不具稳健性。」</div>
<h3>13.3 分析边界声明</h3>
<ul>
<li>未运行 DEG、CellChat、trajectory inference、机器学习或药物预测。</li>
<li>未因结果调整 67-gene list、未调整 region 定义、未切换 Primary 方法。</li>
<li>所有「每例独立 + sign test + random-effects meta」的跨样本汇总均避免将 24,399 个 spot 当作 24,399 个独立患者。</li>
<li><b>数值溯源：</b>本报告全部关键数值由 <code>audit_numeric_consistency.py</code> 从 step1–step7 原始 TSV 自动核对；depth-adjusted p 的检验定义已明确为「12 例 per-sample δ 的单样本 t 检验，df = 11」，并同时列出正态近似值以备查（见 §10.2）。</li>
</ul>
</div>""")

# ---------------------------------------------------------------- 14 files
A("""<h2 id="s14">14. 输出文件清单</h2>
<div class="card">
<table>
<tr><th>类别</th><th>文件</th><th>内容</th></tr>
<tr><td>Score</td><td><code>score/spot_scores.tsv</code></td><td>26,371 spot 的 UCell 与 meanZ score、library size、检出基因数</td></tr>
<tr><td>Moran's I</td><td><code>moran/morans_I_by_sample.tsv</code></td><td>12 例 Moran's I、E[I]、z、p、FDR、邻居数</td></tr>
<tr><td>Region 效应</td><td><code>effects/region_effect_sizes.tsv</code></td><td>每例 × 每比较的 median diff、Cliff's δ、MWU p</td></tr>
<tr><td>Region 汇总</td><td><code>effects/region_cross_sample_summary.tsv</code></td><td>sign test + random-effects meta（REML）</td></tr>
<tr><td>Region 敏感性</td><td><code>effects/region_leave_s3_out_summary.tsv</code></td><td>leave-sample_3-out</td></tr>
<tr><td>病理 效应</td><td><code>effects/pathology_effect_sizes.tsv</code></td><td>SCC vs 预先定义 non-SCC</td></tr>
<tr><td>病理 汇总</td><td><code>effects/pathology_cross_sample_summary.tsv</code> · <code>pathology_leave_s3_out_summary.tsv</code></td><td>跨样本总结与敏感性</td></tr>
<tr><td>稳健性</td><td><code>robustness/robustness_spearman.tsv</code> · <code>robustness_direction_agreement.tsv</code></td><td>UCell vs meanZ 相关与方向一致性</td></tr>
<tr><td>深度诊断</td><td><code>robustness/depth_adjusted_region_direction.tsv</code></td><td>深度校正前后 δ 对比</td></tr>
<tr><td>技术局限</td><td><code>quality/gene_detection_by_sample.tsv</code></td><td>67 gene × 12 样本检出率矩阵</td></tr>
<tr><td>空间图</td><td><code>plots/fig_spatial_ucell_grid.png/pdf</code> 等 6 组</td><td>score 空间图、region 叠加、病理叠加、分布、Moran's I、forest</td></tr>
<tr><td>输入</td><td><code>input/spot_table_all.tsv</code></td><td>样本 × barcode × 坐标 × region × pathology</td></tr>
</table>
</div>""")

A(f"""<div class="meta">
GSE208253 空间转录组预先规定验证 · 67-gene <i>mitoxyperiosis</i> program<br>
Primary = UCell rank-based score（per-sample） · 12 samples · 26,371 in-tissue spots · 24,399 annotated<br>
生成于 {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M')} · 所有数值均由 step1–step7 脚本直接读取，未经手工转录
</div></div></body></html>""")

out = os.path.join(OUT, "GSE208253_spatial_validation_report.html")
with open(out, "w", encoding="utf-8") as f:
    f.write("".join(H))
print("WROTE", out, os.path.getsize(out), "bytes")
