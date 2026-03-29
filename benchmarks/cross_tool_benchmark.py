#!/usr/bin/env python3
"""
Cross-tool benchmark: eskaks vs KaKs_Calculator vs BioPython vs PAML yn00.
Generates accuracy scatter plots, performance bar charts, and saves all results to JSON.
"""

import os
import sys
import json
import math
import time
import subprocess
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
ROOT_DIR = SCRIPT_DIR.parent
ESKAKS_BIN = ROOT_DIR / "target" / "release" / "eskaks"
DATA_DIR = SCRIPT_DIR / "data"
RESULTS_DIR = SCRIPT_DIR / "results"
PLOTS_DIR = SCRIPT_DIR / "plots"
JSON_FILE = SCRIPT_DIR / "cross_tool_results.json"

# ── Benchmark configurations ─────────────────────────────────────────
CONFIGS = {
    "small":  {"n_seqs": 20,  "n_codons": 100,  "label": "20 seqs x 300 bp"},
    "medium": {"n_seqs": 100, "n_codons": 1000, "label": "100 seqs x 3000 bp"},
    "large":  {"n_seqs": 500, "n_codons": 1000, "label": "500 seqs x 3000 bp"},
}

# Tools that are too slow for the large dataset
SKIP_LARGE = {"biopython", "yn00"}


# ── Parsing helpers ───────────────────────────────────────────────────

def parse_eskaks(path):
    """Parse eskaks pairwise TSV -> {(seq1,seq2): (dN, dS)}."""
    results = {}
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                key = (parts[0], parts[1])
                dn = float(parts[2]) if parts[2] not in ('NaN', 'nan') else float('nan')
                ds = float(parts[3]) if parts[3] not in ('NaN', 'nan') else float('nan')
                results[key] = (dn, ds)
    return results


def parse_kakscalc(path):
    """Parse KaKs_Calculator output -> {(seq1,seq2): (Ka, Ks)}."""
    results = {}
    with open(path) as f:
        f.readline()  # header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                pair = parts[0].split('-', 1)
                if len(pair) == 2:
                    try:
                        ka = float(parts[2]) if parts[2] not in ('NA', 'nan') else float('nan')
                        ks = float(parts[3]) if parts[3] not in ('NA', 'nan') else float('nan')
                        results[(pair[0], pair[1])] = (ka, ks)
                    except ValueError:
                        pass
    return results


def parse_biopython(path):
    """Parse BioPython output TSV -> {(seq1,seq2): (dN, dS)}."""
    results = {}
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                key = (parts[0], parts[1])
                dn = float(parts[2]) if parts[2] != 'nan' else float('nan')
                ds = float(parts[3]) if parts[3] != 'nan' else float('nan')
                results[key] = (dn, ds)
    return results


def parse_yn00(work_dir):
    """Parse PAML yn00 output -> {(seq1,seq2): (dN, dS)}.

    yn00 outputs a lower-triangular matrix in the Yang & Nielsen (2000) section.
    We parse the 'Nei & Gojobori 1986' section for NG-comparable values,
    and the 'Yang & Nielsen (2000)' section for YN values.
    Returns (ng_results, yn_results) dicts.
    """
    out_file = os.path.join(work_dir, "yn00_output.txt")
    if not os.path.exists(out_file):
        return {}, {}

    with open(out_file) as f:
        content = f.read()

    ng_results = {}
    yn_results = {}

    # Parse the "(A) Nei-Gojobori (1986) method" section
    # Format: seq2         dS(±SE)    dN(±SE)    omega
    import re
    lines = content.split('\n')
    seq_names = []
    section = None

    for i, line in enumerate(lines):
        if 'Nei & Gojobori' in line or 'Nei-Gojobori' in line:
            section = 'ng'
            continue
        if 'Yang & Nielsen' in line:
            section = 'yn'
            continue
        if section is None:
            continue

        # Detect sequence names (lines with just a name at start)
        stripped = line.strip()
        if not stripped:
            continue

        # yn00 outputs matrix rows like:
        # seq2        0.1234(0.0100)  0.0567(0.0080)  0.4567
        # Each row has the row-sequence name, then pairs of dS(SE) dN(SE) omega
        parts = stripped.split()
        if not parts:
            continue

        # Try to parse as a matrix row
        name = parts[0]
        floats = []
        for p in parts[1:]:
            # Extract number before '(' or the whole thing
            m = re.match(r'^([0-9.eE+-]+)', p)
            if m:
                try:
                    floats.append(float(m.group(1)))
                except ValueError:
                    pass

        # Each pair entry is: dS, dN, omega (3 values per pair)
        if len(floats) >= 3 and len(floats) % 3 == 0:
            n_pairs = len(floats) // 3
            if name not in seq_names:
                seq_names.append(name)
            row_idx = seq_names.index(name)
            for col in range(n_pairs):
                ds_val = floats[col * 3]
                dn_val = floats[col * 3 + 1]
                if col < row_idx:
                    s1 = seq_names[col]
                    s2 = name
                    key = (s1, s2)
                    rev_key = (s2, s1)
                    target = ng_results if section == 'ng' else yn_results
                    target[key] = (dn_val, ds_val)
                    target[rev_key] = (dn_val, ds_val)

    return ng_results, yn_results


# ── Tool runners ──────────────────────────────────────────────────────

def time_cmd(cmd, **kwargs):
    """Run a command and return (elapsed_ms, exit_code)."""
    start = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, **kwargs)
    elapsed_ms = int((time.perf_counter() - start) * 1000)
    return elapsed_ms, result.returncode


def run_eskaks(fasta, model, out_prefix, workers=1):
    cmd = [str(ESKAKS_BIN), fasta, "--model", model, "-o", out_prefix, "--workers", str(workers)]
    return time_cmd(cmd)


def run_kakscalc(axt, out_path, method="NG"):
    cmd = ["KaKs_Calculator", "-i", axt, "-o", out_path, "-m", method]
    return time_cmd(cmd)


def run_biopython(fasta, out_path):
    script = str(SCRIPT_DIR / "biopython_bench.py")
    cmd = ["python3", script, fasta, out_path]
    return time_cmd(cmd)


def run_yn00(phy, work_dir):
    os.makedirs(work_dir, exist_ok=True)
    out_file = os.path.join(work_dir, "yn00_output.txt")
    ctl = f"      seqfile = {phy}\n      outfile = {out_file}\n      verbose = 0\n      icode = 0\n      commonf3x4 = 0\n"
    with open(os.path.join(work_dir, "yn00.ctl"), 'w') as f:
        f.write(ctl)
    cmd = ["yn00", "yn00.ctl"]
    return time_cmd(cmd, cwd=work_dir)


# ── Accuracy comparison ──────────────────────────────────────────────

def compare_pair_data(data_a, data_b):
    """Compare two {(s1,s2): (dN,dS)} dicts. Returns lists of (ref_val, test_val) for dN and dS."""
    common = set(data_a.keys()) & set(data_b.keys())
    dn_pairs = []
    ds_pairs = []
    for key in common:
        a_dn, a_ds = data_a[key]
        b_dn, b_ds = data_b[key]
        if math.isfinite(a_dn) and math.isfinite(b_dn):
            dn_pairs.append((a_dn, b_dn))
        if math.isfinite(a_ds) and math.isfinite(b_ds):
            ds_pairs.append((a_ds, b_ds))
    return dn_pairs, ds_pairs


def accuracy_stats(pairs):
    """Compute accuracy statistics from [(ref, test)] pairs."""
    if not pairs:
        return {"n": 0, "mean_abs_diff": None, "max_abs_diff": None, "r_squared": None}
    diffs = [abs(a - b) for a, b in pairs]
    refs = [a for a, b in pairs]
    tests = [b for a, b in pairs]
    mean_ref = sum(refs) / len(refs) if refs else 0
    ss_tot = sum((r - mean_ref)**2 for r in refs)
    ss_res = sum((r - t)**2 for r, t in pairs)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return {
        "n": len(pairs),
        "mean_abs_diff": sum(diffs) / len(diffs),
        "max_abs_diff": max(diffs),
        "r_squared": r2,
    }


# ── Main benchmark ───────────────────────────────────────────────────

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    all_results = {"timings": {}, "accuracy": {}}

    # ── Generate data ──
    print(">>> Generating benchmark data...")
    for size, cfg in CONFIGS.items():
        bp = cfg["n_codons"] * 3
        fasta = DATA_DIR / f"bench_{cfg['n_seqs']}x{bp}.fasta"
        if not fasta.exists():
            gen_script = str(SCRIPT_DIR / "generate_seqs.py")
            subprocess.run(["python3", gen_script, str(cfg["n_seqs"]), str(cfg["n_codons"]),
                          str(DATA_DIR), "42"], check=True, capture_output=True)

    # ── Run all tools ──
    for size, cfg in CONFIGS.items():
        bp = cfg["n_codons"] * 3
        fasta = str(DATA_DIR / f"bench_{cfg['n_seqs']}x{bp}.fasta")
        axt = str(DATA_DIR / f"bench_{cfg['n_seqs']}x{bp}.axt")
        phy = str(DATA_DIR / f"bench_{cfg['n_seqs']}x{bp}.phy")

        print(f"\n{'='*60}")
        print(f"  {cfg['label']} ({size})")
        print(f"{'='*60}")

        timings = {}

        # eskaks Nei (1 thread)
        prefix = str(RESULTS_DIR / f"eskaks_nei_{size}")
        ms, rc = run_eskaks(fasta, "nei", prefix, workers=1)
        timings["eskaks_nei_1t"] = ms
        print(f"  eskaks Nei (1t):        {ms:>8} ms")

        # eskaks Nei (4 threads)
        prefix_4t = str(RESULTS_DIR / f"eskaks_nei_{size}_4t")
        ms, rc = run_eskaks(fasta, "nei", prefix_4t, workers=4)
        timings["eskaks_nei_4t"] = ms
        print(f"  eskaks Nei (4t):        {ms:>8} ms")

        # eskaks Li (1 thread)
        prefix = str(RESULTS_DIR / f"eskaks_li_{size}")
        ms, rc = run_eskaks(fasta, "li", prefix, workers=1)
        timings["eskaks_li_1t"] = ms
        print(f"  eskaks Li  (1t):        {ms:>8} ms")

        # eskaks Li (4 threads)
        prefix_4t = str(RESULTS_DIR / f"eskaks_li_{size}_4t")
        ms, rc = run_eskaks(fasta, "li", prefix_4t, workers=4)
        timings["eskaks_li_4t"] = ms
        print(f"  eskaks Li  (4t):        {ms:>8} ms")

        # KaKs_Calculator NG
        out = str(RESULTS_DIR / f"kakscalc_ng_{size}.tsv")
        ms, rc = run_kakscalc(axt, out, "NG")
        timings["kakscalc_ng"] = ms
        print(f"  KaKs_Calculator NG:     {ms:>8} ms")

        # KaKs_Calculator LPB
        out = str(RESULTS_DIR / f"kakscalc_lpb_{size}.tsv")
        ms, rc = run_kakscalc(axt, out, "LPB")
        timings["kakscalc_lpb"] = ms
        print(f"  KaKs_Calculator LPB:    {ms:>8} ms")

        # PAML yn00
        if size not in ("large",):
            work_dir = str(RESULTS_DIR / f"yn00_{size}_workdir")
            ms, rc = run_yn00(phy, work_dir)
            timings["yn00"] = ms
            print(f"  PAML yn00:              {ms:>8} ms")
        else:
            print(f"  PAML yn00:              (skipped - too slow)")

        # BioPython
        if size not in ("large",):
            out = str(RESULTS_DIR / f"biopython_ng_{size}.tsv")
            ms, rc = run_biopython(fasta, out)
            timings["biopython_ng"] = ms
            print(f"  BioPython NG86:         {ms:>8} ms")
        else:
            print(f"  BioPython NG86:         (skipped - too slow)")

        all_results["timings"][size] = timings

    # ── Accuracy analysis ──
    print(f"\n{'='*60}")
    print(f"  Accuracy Analysis (small dataset)")
    print(f"{'='*60}")

    acc = {}
    eskaks_nei = parse_eskaks(str(RESULTS_DIR / "eskaks_nei_small_pairwise_results.tsv"))
    eskaks_li = parse_eskaks(str(RESULTS_DIR / "eskaks_li_small_pairwise_results.tsv"))
    kakscalc_ng = parse_kakscalc(str(RESULTS_DIR / "kakscalc_ng_small.tsv"))
    kakscalc_lpb = parse_kakscalc(str(RESULTS_DIR / "kakscalc_lpb_small.tsv"))

    bp_path = RESULTS_DIR / "biopython_ng_small.tsv"
    biopython_ng = parse_biopython(str(bp_path)) if bp_path.exists() else {}

    # Parse yn00 results (NG and YN methods)
    yn00_dir = RESULTS_DIR / "yn00_small_workdir"
    yn00_ng, yn00_yn = parse_yn00(str(yn00_dir)) if yn00_dir.exists() else ({}, {})

    comparisons = {
        "eskaks_nei_vs_kakscalc_ng": (eskaks_nei, kakscalc_ng, "eskaks Nei", "KaKs_Calc NG"),
        "eskaks_li_vs_kakscalc_lpb": (eskaks_li, kakscalc_lpb, "eskaks Li", "KaKs_Calc LPB"),
    }
    if biopython_ng:
        comparisons["eskaks_nei_vs_biopython"] = (eskaks_nei, biopython_ng, "eskaks Nei", "BioPython NG86")
        comparisons["kakscalc_ng_vs_biopython"] = (kakscalc_ng, biopython_ng, "KaKs_Calc NG", "BioPython NG86")
    if yn00_ng:
        comparisons["eskaks_nei_vs_yn00_ng"] = (eskaks_nei, yn00_ng, "eskaks Nei", "PAML yn00 (NG)")
    if yn00_yn:
        comparisons["eskaks_nei_vs_yn00_yn"] = (eskaks_nei, yn00_yn, "eskaks Nei", "PAML yn00 (YN)")

    scatter_data = {}
    for comp_name, (data_a, data_b, label_a, label_b) in comparisons.items():
        dn_pairs, ds_pairs = compare_pair_data(data_a, data_b)
        dn_stats = accuracy_stats(dn_pairs)
        ds_stats = accuracy_stats(ds_pairs)
        acc[comp_name] = {"dN": dn_stats, "dS": ds_stats, "label_a": label_a, "label_b": label_b}
        scatter_data[comp_name] = {"dn": dn_pairs, "ds": ds_pairs, "label_a": label_a, "label_b": label_b}

        print(f"\n  {label_a} vs {label_b}:")
        print(f"    dN: n={dn_stats['n']}, mean_diff={dn_stats.get('mean_abs_diff', 0):.6f}, "
              f"max_diff={dn_stats.get('max_abs_diff', 0):.6f}, R²={dn_stats.get('r_squared', 0):.6f}")
        print(f"    dS: n={ds_stats['n']}, mean_diff={ds_stats.get('mean_abs_diff', 0):.6f}, "
              f"max_diff={ds_stats.get('max_abs_diff', 0):.6f}, R²={ds_stats.get('r_squared', 0):.6f}")

    all_results["accuracy"] = acc

    # ── Save JSON ──
    with open(JSON_FILE, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to: {JSON_FILE}")

    # ── Generate plots ──
    generate_accuracy_plots(scatter_data)
    generate_performance_plots(all_results["timings"])

    print(f"\nAll plots saved to: {PLOTS_DIR}/")


def generate_accuracy_plots(scatter_data):
    """Generate accuracy scatter plots: eskaks vs each reference tool."""
    n_comps = len(scatter_data)
    fig, axes = plt.subplots(n_comps, 2, figsize=(12, 5 * n_comps), squeeze=False)
    fig.suptitle("Numerical Accuracy: eskaks vs Reference Tools\n(20 seqs x 300 bp, small dataset)",
                 fontsize=14, fontweight="bold", y=1.01)

    for row, (comp_name, data) in enumerate(scatter_data.items()):
        label_a = data["label_a"]
        label_b = data["label_b"]

        for col, (metric, pairs) in enumerate([("dN", data["dn"]), ("dS", data["ds"])]):
            ax = axes[row][col]
            if not pairs:
                ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{metric}: {label_a} vs {label_b}")
                continue

            xs = [p[0] for p in pairs]
            ys = [p[1] for p in pairs]
            vmax = max(max(xs), max(ys)) * 1.1 if xs else 1

            ax.scatter(xs, ys, alpha=0.4, s=15, c="#2196F3", edgecolors="none")
            ax.plot([0, vmax], [0, vmax], "--", color="red", alpha=0.7, linewidth=1, label="y = x")
            ax.set_xlabel(label_a)
            ax.set_ylabel(label_b)
            ax.set_title(f"{metric}: {label_a} vs {label_b}")
            ax.set_xlim(0, vmax)
            ax.set_ylim(0, vmax)
            ax.set_aspect("equal")
            ax.legend(loc="upper left", fontsize=8)

            stats = accuracy_stats(pairs)
            ax.text(0.97, 0.03,
                    f"n={stats['n']}\nR²={stats['r_squared']:.6f}\nmax|diff|={stats['max_abs_diff']:.4f}",
                    transform=ax.transAxes, fontsize=8, va="bottom", ha="right",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

    plt.tight_layout()
    path = PLOTS_DIR / "accuracy_scatter.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Accuracy plot: {path}")


def generate_performance_plots(timings):
    """Generate performance comparison bar charts."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Performance: eskaks vs Reference Tools", fontsize=14, fontweight="bold")

    colors = {
        "eskaks_nei_1t": "#1565C0",
        "eskaks_nei_4t": "#42A5F5",
        "eskaks_li_1t":  "#C62828",
        "eskaks_li_4t":  "#EF5350",
        "kakscalc_ng":   "#2E7D32",
        "kakscalc_lpb":  "#66BB6A",
        "yn00":          "#F57F17",
        "biopython_ng":  "#7B1FA2",
    }
    labels = {
        "eskaks_nei_1t": "eskaks Nei (1t)",
        "eskaks_nei_4t": "eskaks Nei (4t)",
        "eskaks_li_1t":  "eskaks Li (1t)",
        "eskaks_li_4t":  "eskaks Li (4t)",
        "kakscalc_ng":   "KaKs_Calc NG",
        "kakscalc_lpb":  "KaKs_Calc LPB",
        "yn00":          "PAML yn00",
        "biopython_ng":  "BioPython NG86",
    }

    for idx, (size, size_label) in enumerate([
        ("small", "Small (20x300bp)"),
        ("medium", "Medium (100x3000bp)"),
        ("large", "Large (500x3000bp)")
    ]):
        ax = axes[idx]
        if size not in timings:
            ax.set_title(size_label)
            continue

        data = timings[size]
        # Sort tools: eskaks first, then external
        tool_order = ["eskaks_nei_1t", "eskaks_nei_4t", "eskaks_li_1t", "eskaks_li_4t",
                      "kakscalc_ng", "kakscalc_lpb", "yn00", "biopython_ng"]
        tools_present = [t for t in tool_order if t in data]
        vals = [data[t] for t in tools_present]
        tool_labels = [labels[t] for t in tools_present]
        tool_colors = [colors[t] for t in tools_present]

        bars = ax.barh(range(len(tools_present)), vals, color=tool_colors, edgecolor="white", height=0.7)
        ax.set_yticks(range(len(tools_present)))
        ax.set_yticklabels(tool_labels, fontsize=9)
        ax.set_xlabel("Time (ms)")
        ax.set_title(size_label)
        ax.invert_yaxis()

        # Add value labels on bars
        for bar, val in zip(bars, vals):
            ax.text(bar.get_width() + max(vals) * 0.01, bar.get_y() + bar.get_height() / 2,
                    f"{val:,} ms", va="center", fontsize=8)

        # Log scale if range is > 100x
        if max(vals) > 100 * min(v for v in vals if v > 0):
            ax.set_xscale("log")
            ax.set_xlabel("Time (ms, log scale)")

    plt.tight_layout()
    path = PLOTS_DIR / "performance_bars.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Performance plot: {path}")

    # ── Speedup chart: eskaks vs each tool ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Speedup: eskaks (4 threads) vs Reference Tools", fontsize=14, fontweight="bold")

    for col, (eskaks_key, ext_tools, title) in enumerate([
        ("eskaks_nei_4t", [("kakscalc_ng", "KaKs_Calc NG"), ("biopython_ng", "BioPython NG86"), ("yn00", "PAML yn00")],
         "Nei-Gojobori Model"),
        ("eskaks_li_4t", [("kakscalc_lpb", "KaKs_Calc LPB")],
         "Li (1993) Model"),
    ]):
        ax = axes[col]
        sizes = ["small", "medium"]
        x_pos = np.arange(len(sizes))
        width = 0.25
        offset = 0

        for ext_key, ext_label in ext_tools:
            speedups = []
            for s in sizes:
                if s in timings and eskaks_key in timings[s] and ext_key in timings[s]:
                    eskaks_ms = timings[s][eskaks_key]
                    ext_ms = timings[s][ext_key]
                    speedups.append(ext_ms / eskaks_ms if eskaks_ms > 0 else 0)
                else:
                    speedups.append(0)

            bars = ax.bar(x_pos + offset * width, speedups, width, label=f"vs {ext_label}",
                         edgecolor="white")
            for bar, sp in zip(bars, speedups):
                if sp > 0:
                    ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                            f"{sp:.0f}x", ha="center", va="bottom", fontsize=9, fontweight="bold")
            offset += 1

        ax.set_xticks(x_pos + width * (len(ext_tools) - 1) / 2)
        ax.set_xticklabels([CONFIGS[s]["label"] for s in sizes], fontsize=9)
        ax.set_ylabel("Speedup (x times faster)")
        ax.set_title(title)
        ax.legend(fontsize=9)
        ax.axhline(y=1, color="gray", linestyle="--", alpha=0.5)

    plt.tight_layout()
    path = PLOTS_DIR / "speedup_chart.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Speedup plot: {path}")


if __name__ == "__main__":
    main()
