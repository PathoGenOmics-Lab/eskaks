#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# eskaks Benchmark Suite
# Compares eskaks against KaKs_Calculator, PAML yn00, and BioPython
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ESKAKS_BIN="${SCRIPT_DIR}/../target/release/eskaks"
BENCH_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results"

mkdir -p "$BENCH_DIR" "$RESULTS_DIR"

# --- Configuration ---
# Small benchmark for accuracy comparison
SMALL_N=20
SMALL_L=100   # codons (300 bp)

# Medium benchmark for performance
MED_N=100
MED_L=1000    # codons (3000 bp)

# Large benchmark (eskaks + KaKs_Calculator only; yn00/biopython too slow)
LARGE_N=500
LARGE_L=1000  # codons (3000 bp)

echo "============================================"
echo " eskaks Benchmark Suite"
echo "============================================"
echo ""

# --- Generate test data ---
echo ">>> Generating benchmark sequences..."
python3 "${SCRIPT_DIR}/generate_seqs.py" $SMALL_N $SMALL_L "$BENCH_DIR" 42
python3 "${SCRIPT_DIR}/generate_seqs.py" $MED_N $MED_L "$BENCH_DIR" 42
python3 "${SCRIPT_DIR}/generate_seqs.py" $LARGE_N $LARGE_L "$BENCH_DIR" 42
echo ""

# --- Helper: time a command and print result ---
time_cmd() {
    local label="$1"
    shift
    local start end elapsed
    start=$(date +%s%N)
    "$@" > /dev/null 2>&1
    end=$(date +%s%N)
    elapsed=$(( (end - start) / 1000000 ))
    printf "  %-45s %8d ms\n" "$label" "$elapsed"
    echo "$label,$elapsed" >> "${RESULTS_DIR}/timings.csv"
}

# ============================================================================
# PART 1: Numerical Accuracy (small dataset)
# ============================================================================
echo "============================================"
echo " PART 1: Numerical Accuracy (${SMALL_N} seqs × ${SMALL_L} codons)"
echo "============================================"
SMALL_FASTA="${BENCH_DIR}/bench_${SMALL_N}x$((SMALL_L*3)).fasta"
SMALL_AXT="${BENCH_DIR}/bench_${SMALL_N}x$((SMALL_L*3)).axt"
SMALL_PHY="${BENCH_DIR}/bench_${SMALL_N}x$((SMALL_L*3)).phy"

echo ""
echo "--- eskaks (Nei-Gojobori) ---"
"$ESKAKS_BIN" "$SMALL_FASTA" --model nei -o "${RESULTS_DIR}/eskaks_nei_small" --workers 1 2>/dev/null
echo "  Output: ${RESULTS_DIR}/eskaks_nei_small_pairwise_results.tsv"

echo ""
echo "--- eskaks (Li 1993) ---"
"$ESKAKS_BIN" "$SMALL_FASTA" --model li -o "${RESULTS_DIR}/eskaks_li_small" --workers 1 2>/dev/null
echo "  Output: ${RESULTS_DIR}/eskaks_li_small_pairwise_results.tsv"

echo ""
echo "--- KaKs_Calculator (NG) ---"
KaKs_Calculator -i "$SMALL_AXT" -o "${RESULTS_DIR}/kakscalc_ng_small.tsv" -m NG 2>/dev/null
echo "  Output: ${RESULTS_DIR}/kakscalc_ng_small.tsv"

echo ""
echo "--- KaKs_Calculator (LPB / Li 1993) ---"
KaKs_Calculator -i "$SMALL_AXT" -o "${RESULTS_DIR}/kakscalc_lpb_small.tsv" -m LPB 2>/dev/null
echo "  Output: ${RESULTS_DIR}/kakscalc_lpb_small.tsv"

echo ""
echo "--- PAML yn00 ---"
# yn00 requires a control file and runs in the current directory
YN00_DIR="${RESULTS_DIR}/yn00_workdir"
mkdir -p "$YN00_DIR"
cat > "${YN00_DIR}/yn00.ctl" <<CTLEOF
      seqfile = ${SMALL_PHY}
      outfile = ${YN00_DIR}/yn00_output.txt
      verbose = 0
      icode = 0
      commonf3x4 = 0
CTLEOF
(cd "$YN00_DIR" && yn00 yn00.ctl > /dev/null 2>&1) || echo "  (yn00 returned non-zero, output may be partial)"
echo "  Output: ${YN00_DIR}/yn00_output.txt"

echo ""
echo "--- BioPython (NG86) ---"
python3 "${SCRIPT_DIR}/biopython_bench.py" "$SMALL_FASTA" "${RESULTS_DIR}/biopython_ng_small.tsv" 2>/dev/null || echo "  (BioPython failed)"

# ============================================================================
# PART 2: Performance Benchmark
# ============================================================================
echo ""
echo "============================================"
echo " PART 2: Performance Benchmark"
echo "============================================"

# Clear timings file
echo "tool,time_ms" > "${RESULTS_DIR}/timings.csv"

# --- Medium dataset (100 seqs × 3000 bp) ---
echo ""
echo "--- Medium: ${MED_N} seqs × $((MED_L*3)) bp ---"
MED_FASTA="${BENCH_DIR}/bench_${MED_N}x$((MED_L*3)).fasta"
MED_AXT="${BENCH_DIR}/bench_${MED_N}x$((MED_L*3)).axt"
MED_PHY="${BENCH_DIR}/bench_${MED_N}x$((MED_L*3)).phy"

time_cmd "eskaks Nei (1 thread)" "$ESKAKS_BIN" "$MED_FASTA" --model nei -o "${RESULTS_DIR}/eskaks_nei_med" --workers 1
time_cmd "eskaks Nei (4 threads)" "$ESKAKS_BIN" "$MED_FASTA" --model nei -o "${RESULTS_DIR}/eskaks_nei_med_4t" --workers 4
time_cmd "eskaks Li (1 thread)" "$ESKAKS_BIN" "$MED_FASTA" --model li -o "${RESULTS_DIR}/eskaks_li_med" --workers 1
time_cmd "eskaks Li (4 threads)" "$ESKAKS_BIN" "$MED_FASTA" --model li -o "${RESULTS_DIR}/eskaks_li_med_4t" --workers 4
time_cmd "KaKs_Calculator NG" KaKs_Calculator -i "$MED_AXT" -o "${RESULTS_DIR}/kakscalc_ng_med.tsv" -m NG
time_cmd "KaKs_Calculator LPB (Li 1993)" KaKs_Calculator -i "$MED_AXT" -o "${RESULTS_DIR}/kakscalc_lpb_med.tsv" -m LPB

# yn00 on medium
YN00_MED_DIR="${RESULTS_DIR}/yn00_med_workdir"
mkdir -p "$YN00_MED_DIR"
cat > "${YN00_MED_DIR}/yn00.ctl" <<CTLEOF
      seqfile = ${MED_PHY}
      outfile = ${YN00_MED_DIR}/yn00_output.txt
      verbose = 0
      icode = 0
      commonf3x4 = 0
CTLEOF
time_cmd "PAML yn00" bash -c "cd '$YN00_MED_DIR' && yn00 yn00.ctl"

time_cmd "BioPython NG86" python3 "${SCRIPT_DIR}/biopython_bench.py" "$MED_FASTA" "${RESULTS_DIR}/biopython_ng_med.tsv"

# --- Large dataset (500 seqs × 3000 bp) ---
echo ""
echo "--- Large: ${LARGE_N} seqs × $((LARGE_L*3)) bp ---"
LARGE_FASTA="${BENCH_DIR}/bench_${LARGE_N}x$((LARGE_L*3)).fasta"
LARGE_AXT="${BENCH_DIR}/bench_${LARGE_N}x$((LARGE_L*3)).axt"

time_cmd "eskaks Nei (1 thread)" "$ESKAKS_BIN" "$LARGE_FASTA" --model nei -o "${RESULTS_DIR}/eskaks_nei_large" --workers 1
time_cmd "eskaks Nei (4 threads)" "$ESKAKS_BIN" "$LARGE_FASTA" --model nei -o "${RESULTS_DIR}/eskaks_nei_large_4t" --workers 4
time_cmd "eskaks Li (1 thread)" "$ESKAKS_BIN" "$LARGE_FASTA" --model li -o "${RESULTS_DIR}/eskaks_li_large" --workers 1
time_cmd "eskaks Li (4 threads)" "$ESKAKS_BIN" "$LARGE_FASTA" --model li -o "${RESULTS_DIR}/eskaks_li_large_4t" --workers 4
time_cmd "KaKs_Calculator NG" KaKs_Calculator -i "$LARGE_AXT" -o "${RESULTS_DIR}/kakscalc_ng_large.tsv" -m NG
time_cmd "KaKs_Calculator LPB (Li 1993)" KaKs_Calculator -i "$LARGE_AXT" -o "${RESULTS_DIR}/kakscalc_lpb_large.tsv" -m LPB

echo ""
echo "(Skipping yn00 and BioPython for large dataset — too slow for ${LARGE_N} seqs)"

# ============================================================================
# PART 3: Summary
# ============================================================================
echo ""
echo "============================================"
echo " Results Summary"
echo "============================================"
echo ""
cat "${RESULTS_DIR}/timings.csv" | column -t -s','
echo ""
echo "All outputs saved to: ${RESULTS_DIR}/"
