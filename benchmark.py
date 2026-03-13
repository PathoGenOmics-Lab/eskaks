#!/usr/bin/env python3
"""
Benchmark suite for eskaks — dN/dS calculator.
Measures execution time and peak memory across different:
  - Sequence counts (N)
  - Sequence lengths (L)
  - Models (Nei vs Li)
  - Thread counts
Generates plots saved to benchmark_results/.
"""

import os
import random
import subprocess
import time
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ── Configuration ──────────────────────────────────────────────────────
ESKAKS = "./target/release/eskaks"
OUT_DIR = Path("benchmark_results")
FASTA_DIR = OUT_DIR / "fasta"
RESULTS_FILE = OUT_DIR / "results.json"

# Benchmark dimensions
SEQ_COUNTS = [50, 200, 500, 1000]
SEQ_LENGTHS = [3000, 30000, 120000, 300000]  # nucleotides (must be divisible by 3)
MODELS = ["nei", "li"]
THREAD_COUNTS = [1, 2, 4, 8]

CODONS = [
    "ATG", "GCT", "GCC", "GCA", "GCG",  # Met, Ala x4
    "CTT", "CTC", "CTA", "CTG",          # Leu x4
    "GTT", "GTC", "GTA", "GTG",          # Val x4
    "TCT", "TCC", "TCA", "TCG",          # Ser x4
    "CCT", "CCC", "CCA", "CCG",          # Pro x4
    "ACT", "ACC", "ACA", "ACG",          # Thr x4
    "GAT", "GAC",                         # Asp x2
    "GAA", "GAG",                         # Glu x2
    "AAA", "AAG",                         # Lys x2
    "TTC", "TTT",                         # Phe x2
    "TAT", "TAC",                         # Tyr x2
    "TGT", "TGC",                         # Cys x2
    "CAT", "CAC",                         # His x2
    "CAA", "CAG",                         # Gln x2
    "AAT", "AAC",                         # Asn x2
    "ATT", "ATC", "ATA",                  # Ile x3
    "TGG",                                # Trp
    "AGT", "AGC",                         # Ser x2
    "AGA", "AGG",                         # Arg x2
    "CGT", "CGC", "CGA", "CGG",          # Arg x4
    "GGT", "GGC", "GGA", "GGG",          # Gly x4
]


def generate_fasta(n_seqs: int, seq_len: int, filepath: Path):
    """Generate a FASTA file with n_seqs sequences of seq_len nucleotides."""
    random.seed(42)
    n_codons = seq_len // 3

    # Generate a reference sequence
    ref_codons = [random.choice(CODONS) for _ in range(n_codons)]

    with open(filepath, "w") as f:
        for i in range(n_seqs):
            # Mutate ~5% of codons from reference (mix of syn and nonsyn)
            codons = list(ref_codons)
            n_mutations = max(1, n_codons // 20)
            positions = random.sample(range(n_codons), min(n_mutations, n_codons))
            for pos in positions:
                codons[pos] = random.choice(CODONS)
            seq = "".join(codons)
            f.write(f">seq{i:04d}_group{i % 3}\n{seq}\n")


def get_peak_rss_mb(pid: int) -> float:
    """Read peak RSS from /proc/<pid>/status (VmHWM) in MB."""
    try:
        with open(f"/proc/{pid}/status") as f:
            for line in f:
                if line.startswith("VmHWM:"):
                    return int(line.split()[1]) / 1024.0  # kB -> MB
    except (FileNotFoundError, ProcessLookupError):
        pass
    return 0.0


def run_benchmark(fasta_path: str, model: str, workers: int) -> dict:
    """Run eskaks and measure wall time + peak RSS via /proc polling."""
    out_prefix = str(OUT_DIR / "tmp_output")

    cmd = [
        ESKAKS, fasta_path,
        "--model", model,
        "--workers", str(workers),
        "-o", out_prefix,
    ]

    start = time.perf_counter()
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    peak_mb = 0.0
    while proc.poll() is None:
        mb = get_peak_rss_mb(proc.pid)
        if mb > peak_mb:
            peak_mb = mb
        time.sleep(0.05)
    elapsed = time.perf_counter() - start

    # Final read after process exits (VmHWM is the lifetime peak)
    mb = get_peak_rss_mb(proc.pid)
    if mb > peak_mb:
        peak_mb = mb

    # Clean up output files
    for f in OUT_DIR.glob("tmp_output_*"):
        f.unlink()

    return {
        "wall_time_s": round(elapsed, 3),
        "peak_rss_mb": round(peak_mb, 1),
        "exit_code": proc.returncode,
    }


def run_all_benchmarks():
    """Run the full benchmark matrix."""
    OUT_DIR.mkdir(exist_ok=True)
    FASTA_DIR.mkdir(exist_ok=True)

    results = []

    # 1. Scaling by sequence count (fixed length = 120000 nt, 4 threads)
    print("=== Benchmark 1: Scaling by sequence count ===")
    fixed_len = 120000
    for n in SEQ_COUNTS:
        fasta = FASTA_DIR / f"n{n}_l{fixed_len}.fasta"
        if not fasta.exists():
            print(f"  Generating {fasta.name}...")
            generate_fasta(n, fixed_len, fasta)
        for model in MODELS:
            print(f"  Running: N={n}, L={fixed_len}, model={model}, workers=4")
            r = run_benchmark(str(fasta), model, 4)
            r.update({"n_seqs": n, "seq_len": fixed_len, "model": model, "workers": 4, "bench": "seq_count"})
            results.append(r)
            print(f"    -> {r['wall_time_s']}s, {r['peak_rss_mb']}MB")

    # 2. Scaling by sequence length (fixed count = 50, 4 threads)
    print("\n=== Benchmark 2: Scaling by sequence length ===")
    fixed_n = 50
    for l in SEQ_LENGTHS:
        fasta = FASTA_DIR / f"n{fixed_n}_l{l}.fasta"
        if not fasta.exists():
            print(f"  Generating {fasta.name}...")
            generate_fasta(fixed_n, l, fasta)
        for model in MODELS:
            print(f"  Running: N={fixed_n}, L={l}, model={model}, workers=4")
            r = run_benchmark(str(fasta), model, 4)
            r.update({"n_seqs": fixed_n, "seq_len": l, "model": model, "workers": 4, "bench": "seq_len"})
            results.append(r)
            print(f"    -> {r['wall_time_s']}s, {r['peak_rss_mb']}MB")

    # 3. Thread scaling (fixed: 500 seqs, 120000 nt)
    print("\n=== Benchmark 3: Thread scaling ===")
    fixed_n2, fixed_l2 = 500, 120000
    fasta = FASTA_DIR / f"n{fixed_n2}_l{fixed_l2}.fasta"
    if not fasta.exists():
        generate_fasta(fixed_n2, fixed_l2, fasta)
    for w in THREAD_COUNTS:
        for model in MODELS:
            print(f"  Running: N={fixed_n2}, L={fixed_l2}, model={model}, workers={w}")
            r = run_benchmark(str(fasta), model, w)
            r.update({"n_seqs": fixed_n2, "seq_len": fixed_l2, "model": model, "workers": w, "bench": "threads"})
            results.append(r)
            print(f"    -> {r['wall_time_s']}s, {r['peak_rss_mb']}MB")

    # Save raw results
    with open(RESULTS_FILE, "w") as f:
        json.dump(results, f, indent=2)

    return results


def plot_results(results: list):
    """Generate benchmark plots."""
    plt.style.use("seaborn-v0_8-whitegrid")
    colors = {"nei": "#2196F3", "li": "#FF5722"}
    markers = {"nei": "o", "li": "s"}

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("eskaks Benchmark Results", fontsize=16, fontweight="bold")

    # ── Plot 1: Time vs Sequence Count ──
    ax = axes[0][0]
    for model in MODELS:
        data = [r for r in results if r["bench"] == "seq_count" and r["model"] == model]
        data.sort(key=lambda r: r["n_seqs"])
        xs = [r["n_seqs"] for r in data]
        ys = [r["wall_time_s"] for r in data]
        ax.plot(xs, ys, f"-{markers[model]}", color=colors[model], label=model.capitalize(), linewidth=2, markersize=8)
    ax.set_xlabel("Number of Sequences")
    ax.set_ylabel("Wall Time (seconds)")
    ax.set_title("Time vs Sequence Count\n(L=120k nt, 4 threads)")
    ax.legend()

    # ── Plot 2: Time vs Sequence Length ──
    ax = axes[0][1]
    for model in MODELS:
        data = [r for r in results if r["bench"] == "seq_len" and r["model"] == model]
        data.sort(key=lambda r: r["seq_len"])
        xs = [r["seq_len"] // 1000 for r in data]
        ys = [r["wall_time_s"] for r in data]
        ax.plot(xs, ys, f"-{markers[model]}", color=colors[model], label=model.capitalize(), linewidth=2, markersize=8)
    ax.set_xlabel("Sequence Length (k nucleotides)")
    ax.set_ylabel("Wall Time (seconds)")
    ax.set_title("Time vs Sequence Length\n(N=50 seqs, 4 threads)")
    ax.legend()

    # ── Plot 3: Time vs Threads ──
    ax = axes[0][2]
    for model in MODELS:
        data = [r for r in results if r["bench"] == "threads" and r["model"] == model]
        data.sort(key=lambda r: r["workers"])
        xs = [r["workers"] for r in data]
        ys = [r["wall_time_s"] for r in data]
        ax.plot(xs, ys, f"-{markers[model]}", color=colors[model], label=model.capitalize(), linewidth=2, markersize=8)
    # Add ideal scaling line
    data_nei = [r for r in results if r["bench"] == "threads" and r["model"] == "nei"]
    if data_nei:
        data_nei.sort(key=lambda r: r["workers"])
        t1 = data_nei[0]["wall_time_s"]
        ideal_xs = [r["workers"] for r in data_nei]
        ideal_ys = [t1 / w for w in ideal_xs]
        ax.plot(ideal_xs, ideal_ys, "--", color="gray", alpha=0.5, label="Ideal scaling")
    ax.set_xlabel("Number of Threads")
    ax.set_ylabel("Wall Time (seconds)")
    ax.set_title("Thread Scaling\n(N=500, L=120k nt)")
    ax.set_xticks(THREAD_COUNTS)
    ax.legend()

    # ── Plot 4: Memory vs Sequence Count ──
    ax = axes[1][0]
    for model in MODELS:
        data = [r for r in results if r["bench"] == "seq_count" and r["model"] == model]
        data.sort(key=lambda r: r["n_seqs"])
        xs = [r["n_seqs"] for r in data]
        ys = [r["peak_rss_mb"] for r in data]
        ax.plot(xs, ys, f"-{markers[model]}", color=colors[model], label=model.capitalize(), linewidth=2, markersize=8)
    ax.set_xlabel("Number of Sequences")
    ax.set_ylabel("Peak RSS (MB)")
    ax.set_title("Memory vs Sequence Count\n(L=120k nt, 4 threads)")
    ax.legend()

    # ── Plot 5: Memory vs Sequence Length ──
    ax = axes[1][1]
    for model in MODELS:
        data = [r for r in results if r["bench"] == "seq_len" and r["model"] == model]
        data.sort(key=lambda r: r["seq_len"])
        xs = [r["seq_len"] // 1000 for r in data]
        ys = [r["peak_rss_mb"] for r in data]
        ax.plot(xs, ys, f"-{markers[model]}", color=colors[model], label=model.capitalize(), linewidth=2, markersize=8)
    ax.set_xlabel("Sequence Length (k nucleotides)")
    ax.set_ylabel("Peak RSS (MB)")
    ax.set_title("Memory vs Sequence Length\n(N=50 seqs, 4 threads)")
    ax.legend()

    # ── Plot 6: Throughput (pairs/sec) vs Threads ──
    ax = axes[1][2]
    for model in MODELS:
        data = [r for r in results if r["bench"] == "threads" and r["model"] == model]
        data.sort(key=lambda r: r["workers"])
        n = data[0]["n_seqs"] if data else 100
        total_pairs = n * (n - 1) // 2
        xs = [r["workers"] for r in data]
        ys = [total_pairs / r["wall_time_s"] if r["wall_time_s"] > 0 else 0 for r in data]
        ax.plot(xs, ys, f"-{markers[model]}", color=colors[model], label=model.capitalize(), linewidth=2, markersize=8)
    ax.set_xlabel("Number of Threads")
    ax.set_ylabel("Pairs / second")
    ax.set_title("Throughput vs Threads\n(N=500, L=120k nt)")
    ax.set_xticks(THREAD_COUNTS)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
    ax.legend()

    plt.tight_layout()
    plot_path = OUT_DIR / "benchmark_plots.png"
    fig.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nPlots saved to: {plot_path}")

    # ── Summary table ──
    print("\n" + "=" * 80)
    print(f"{'Bench':<12} {'Model':<6} {'N':>5} {'L(nt)':>8} {'Thr':>4} {'Time(s)':>9} {'RSS(MB)':>9}")
    print("-" * 80)
    for r in results:
        print(f"{r['bench']:<12} {r['model']:<6} {r['n_seqs']:>5} {r['seq_len']:>8} {r['workers']:>4} {r['wall_time_s']:>9.3f} {r['peak_rss_mb']:>9.1f}")
    print("=" * 80)


if __name__ == "__main__":
    print(f"eskaks benchmark suite")
    print(f"Binary: {ESKAKS}")
    print()

    results = run_all_benchmarks()
    plot_results(results)
