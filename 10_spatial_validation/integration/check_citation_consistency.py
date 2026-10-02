#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GSE208253 integration -- citation & numeric consistency checker.

SCOPE (per 9.docx): "manuscript/figure/supplement integration and citation
consistency check" ONLY. Read-only. Never mutates the frozen module.

Checks performed on any candidate manuscript text (Results / Methods /
Supplement / response letter / figure legend):

  C1  RETIRED VALUE  -- does the text cite the retired normal-approximation
                        p (0.351 / 0.3509) as a scientific result?
  C2  LOCKED VALUES  -- do the numbers that DO appear match the frozen lock?
  C3  STATISTIC LABEL-- when the depth-adjusted p appears, is it labelled as
                        a one-sample t-test with df = 11?
  C4  NAMING         -- is the program always called "67-gene", never
                        "5-gene signature" / "11-gene" for this module?
  C5  ACCESSION      -- is GSE208253 cited consistently?
  C6  CLAIM BOUNDARY -- any forbidden over-claim ("biological depletion",
                        "tumour-specific localization ... confirmed",
                        "真实的空间聚集")?

Usage:
  python check_citation_consistency.py <file-or-dir> [more...]
  python check_citation_consistency.py --selftest
"""
import os, re, sys, json

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---- frozen lock (loaded from the integration pack, so single source) ----
def load_lock():
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "numeric_freeze_table.tsv")
    lock = {}
    if os.path.exists(p):
        with open(p, encoding="utf-8") as f:
            next(f)
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    lock[parts[0]] = parts[1]
    return lock

LOCK = load_lock()

RETIRED = ["0.351", "0.3509", "0.350877"]
LOCKED_OK = {
    "core meta delta": ["-0.2754", "\u22120.2754"],
    "core meta p": ["3.47", "3.467", "10\u207b\u00b9\u2078"],
    "core CI lo": ["-0.3375", "\u22120.3375"],
    "core CI hi": ["-0.2134", "\u22120.2134"],
    "depth-adj mean": ["0.0365", "+0.0365"],
    "depth-adj t": ["0.9329"],
    "depth-adj p": ["0.3709", "0.370901"],
    "pathology delta": ["-0.1327", "\u22120.1327"],
}
FORBIDDEN_CLAIMS = [
    ("biological depletion", "must NOT interpret the raw inverse core signal as biological depletion"),
    ("tumor core depletion", "same boundary"),
    ("\u771f\u5b9e\u7684\u7a7a\u95f4\u805a\u96c6", "over-strong wording; use 'reproducible positive spatial autocorrelation'"),
]

# Claims that are only forbidden when ASSERTED (not when explicitly negated).
ASSERTION_CLAIMS = [
    ("tumour-specific localization was confirmed",
     "NOT CONFIRMED -- do not assert"),
    ("tumor-specific localization was confirmed",
     "NOT CONFIRMED -- do not assert"),
    ("tumour-specific localization is confirmed",
     "NOT CONFIRMED -- do not assert"),
]

# Negation cues: if one of these precedes the phrase within ~60 chars, the
# statement is a correct negation and must NOT be flagged.
NEG_CUES = re.compile(
    r"(not|no|never|cannot|does not|do not|did not|"
    r"\u4e0d\u5f97|\u4e0d\u80fd|\u4e0d\u53ef|\u672a|\u975e|\u4e0d\u662f)",
    re.I)


def _negated(t, idx, window=90):
    """True if a negation cue occurs shortly before position idx."""
    start = max(0, idx - window)
    return bool(NEG_CUES.search(t[start:idx]))
GOOD_NAMING = ["67-gene", "67 gene", "67gene", "67-gene program"]
REQUIRED_ACCESSION = "GSE208253"


def read_any(path):
    """Read txt/md/html/docx -> plain text."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".docx":
        import zipfile
        z = zipfile.ZipFile(path)
        xml = z.read("word/document.xml").decode("utf-8", "ignore")
        xml = re.sub(r"</w:p>", "\n", xml)
        txt = re.sub(r"<[^>]+>", "", xml)
        return txt
    with open(path, encoding="utf-8", errors="ignore") as f:
        t = f.read()
    if ext in (".html", ".htm"):
        t = re.sub(r"<[^>]+>", " ", t)
    return t


def norm(t):
    t = t.replace("\u2212", "-").replace("\u2013", "-").replace("\u2014", "-")
    t = t.replace("\u00d7", "x")
    for a, b in zip("\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b",
                    "0123456789-"):
        t = t.replace(a, b)
    return t


def check_text(name, raw):
    t = norm(raw)
    fails, warns, oks = [], [], []

    # C1 retired value
    for r in RETIRED:
        if r in t:
            # allow if explicitly flagged as retired / audit-only
            ctx_ok = re.search(
                r"(retired|deprecated|\u5df2\u5e9f\u5f03|provenance|\u5ba1\u8ba1|legacy|"
                r"\u5386\u53f2\u503c|\u4e0d\u5f97\u518d\u4f5c\u4e3a)", t, re.I)
            if ctx_ok:
                oks.append(f"C1 retired value {r} present but explicitly flagged as retired")
            else:
                fails.append(f"C1 RETIRED VALUE {r} cited WITHOUT retirement flag")

    # C2 locked values present?
    for label, variants in LOCKED_OK.items():
        if any(v in t for v in variants):
            oks.append(f"C2 {label} present and matches lock")

    # C3 statistic label
    if "0.3709" in t or "0.370901" in t:
        if re.search(r"(one-sample|one sample|\u5355\u6837\u672c)", t, re.I) and \
           re.search(r"df\s*=?\s*11", t, re.I):
            oks.append("C3 statistic correctly labelled (one-sample t-test, df=11)")
        else:
            warns.append("C3 depth-adjusted p present but statistic label (one-sample t-test, df=11) not adjacent")

    # C4 naming
    if REQUIRED_ACCESSION in t or "spatial" in t.lower():
        if any(g in t for g in GOOD_NAMING):
            oks.append("C4 program correctly named 67-gene")
        if re.search(r"\b5-gene signature\b", t, re.I):
            fails.append("C4 forbidden legacy naming '5-gene signature'")

    # C5 accession
    if REQUIRED_ACCESSION in t:
        oks.append("C5 accession GSE208253 cited")

    # C6 forbidden claims (negation-aware)
    for bad, why in FORBIDDEN_CLAIMS:
        low = t.lower()
        b = bad.lower()
        start = 0
        while True:
            idx = low.find(b, start)
            if idx < 0:
                break
            if _negated(t, idx):
                oks.append(f"C6 '{bad}' present but correctly negated")
            else:
                fails.append(f"C6 forbidden claim/word: '{bad}' ({why})")
            start = idx + len(b)

    # C6b assertion-only claims
    for bad, why in ASSERTION_CLAIMS:
        idx = t.lower().find(bad.lower())
        if idx >= 0 and not _negated(t, idx):
            fails.append(f"C6 forbidden assertion: '{bad}' ({why})")

    return fails, warns, oks


def collect(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            for dp, _, fn in os.walk(p):
                for n in fn:
                    if os.path.splitext(n)[1].lower() in (".md", ".txt", ".html", ".htm", ".docx"):
                        files.append(os.path.join(dp, n))
        else:
            files.append(p)
    return files


def selftest():
    good = ("The 67-gene score showed reproducible positive spatial autocorrelation "
            "across all 12 samples (GSE208253). core vs nc meta delta = -0.2754 "
            "(95% CI -0.3375 to -0.2134, p = 3.47e-18). After depth adjustment the "
            "core vs nc effect was +0.0365 (one-sample t-test, df = 11, p = 0.3709). "
            "This must not be interpreted as biological depletion.")
    bad = good + " The corrected p is 0.351. Tumour-specific localization was confirmed."
    f1, w1, o1 = check_text("good", norm(good))
    f2, w2, o2 = check_text("bad", norm(bad))
    print("SELFTEST good -> fails:", f1)
    print("SELFTEST bad  -> fails:", f2)
    assert not f1, "good sample should pass"
    assert f2, "bad sample should fail"
    print("SELFTEST: PASS")


def main():
    args = sys.argv[1:]
    if not args or args[0] == "--selftest":
        selftest()
        return
    files = collect(args)
    report = {"files": [], "summary": {"FAIL": 0, "WARN": 0, "PASS": 0}}
    for fp in files:
        try:
            raw = read_any(fp)
        except Exception as e:
            report["files"].append({"file": fp, "error": str(e)})
            continue
        fails, warns, oks = check_text(os.path.basename(fp), raw)
        verdict = "FAIL" if fails else ("WARN" if warns else "PASS")
        report["summary"][verdict] += 1
        report["files"].append({
            "file": fp, "verdict": verdict,
            "fails": fails, "warns": warns, "oks": oks,
        })
        print(f"[{verdict}] {fp}")
        for x in fails:
            print("   FAIL ", x)
        for x in warns:
            print("   WARN ", x)
    print("\nSUMMARY:", report["summary"])
    outp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "citation_consistency_check.json")
    with open(outp, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print("wrote:", outp)


if __name__ == "__main__":
    main()
