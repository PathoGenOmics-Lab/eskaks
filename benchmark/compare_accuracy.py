#!/usr/bin/env python3
"""Compare numerical accuracy between eskaks and KaKs_Calculator outputs.
Usage: compare_accuracy.py <results_dir>
"""
import sys
import os
import re
import math


def parse_eskaks(path):
    """Parse eskaks pairwise TSV output."""
    results = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                key = (parts[0], parts[1])
                dn = float(parts[2]) if parts[2] != 'NaN' else float('nan')
                ds = float(parts[3]) if parts[3] != 'NaN' else float('nan')
                results[key] = (dn, ds)
    return results


def parse_kakscalc(path):
    """Parse KaKs_Calculator output TSV."""
    results = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                # Sequence pair name is "seq1-seq2"
                pair_name = parts[0]
                match = pair_name.split('-', 1)
                if len(match) == 2:
                    key = (match[0], match[1])
                    try:
                        ka = float(parts[2]) if parts[2] != 'NA' and parts[2] != 'nan' else float('nan')
                        ks = float(parts[3]) if parts[3] != 'NA' and parts[3] != 'nan' else float('nan')
                        results[key] = (ka, ks)
                    except (ValueError, IndexError):
                        pass
    return results


def parse_biopython(path):
    """Parse BioPython output TSV."""
    results = {}
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                key = (parts[0], parts[1])
                dn = float(parts[2]) if parts[2] != 'nan' else float('nan')
                ds = float(parts[3]) if parts[3] != 'nan' else float('nan')
                results[key] = (dn, ds)
    return results


def compare(ref_name, ref_data, test_name, test_data):
    """Compare two sets of results and compute statistics."""
    common_keys = set(ref_data.keys()) & set(test_data.keys())
    if not common_keys:
        print(f"  No common pairs between {ref_name} and {test_name}")
        return

    dn_diffs = []
    ds_diffs = []
    dn_rel_diffs = []
    ds_rel_diffs = []
    mismatches = 0

    for key in sorted(common_keys):
        ref_dn, ref_ds = ref_data[key]
        test_dn, test_ds = test_data[key]

        # Skip if either is NaN
        if math.isnan(ref_dn) or math.isnan(test_dn):
            continue
        if math.isnan(ref_ds) or math.isnan(test_ds):
            continue

        dn_diff = abs(ref_dn - test_dn)
        ds_diff = abs(ref_ds - test_ds)
        dn_diffs.append(dn_diff)
        ds_diffs.append(ds_diff)

        if ref_dn > 0.001:
            dn_rel_diffs.append(dn_diff / ref_dn)
        if ref_ds > 0.001:
            ds_rel_diffs.append(ds_diff / ref_ds)

        if dn_diff > 0.01 or ds_diff > 0.01:
            mismatches += 1

    n = len(dn_diffs)
    if n == 0:
        print(f"  No finite pairs to compare between {ref_name} and {test_name}")
        return

    print(f"\n  {ref_name} vs {test_name}  ({n} finite pairs compared)")
    print(f"  {'':>20} {'Mean AbsDiff':>14} {'Max AbsDiff':>14} {'Mean RelDiff':>14} {'Max RelDiff':>14}")

    dn_mean = sum(dn_diffs) / n
    dn_max = max(dn_diffs)
    ds_mean = sum(ds_diffs) / n
    ds_max = max(ds_diffs)

    dn_rel_mean = sum(dn_rel_diffs) / len(dn_rel_diffs) if dn_rel_diffs else 0
    dn_rel_max = max(dn_rel_diffs) if dn_rel_diffs else 0
    ds_rel_mean = sum(ds_rel_diffs) / len(ds_rel_diffs) if ds_rel_diffs else 0
    ds_rel_max = max(ds_rel_diffs) if ds_rel_diffs else 0

    print(f"  {'dN':>20} {dn_mean:>14.6f} {dn_max:>14.6f} {dn_rel_mean:>14.4%} {dn_rel_max:>14.4%}")
    print(f"  {'dS':>20} {ds_mean:>14.6f} {ds_max:>14.6f} {ds_rel_mean:>14.4%} {ds_rel_max:>14.4%}")

    if mismatches > 0:
        print(f"  Pairs with |diff| > 0.01: {mismatches}/{n}")
    else:
        print(f"  All pairs agree within 0.01 tolerance")


def main():
    results_dir = sys.argv[1] if len(sys.argv) > 1 else "benchmark/results"

    print("=" * 70)
    print(" Numerical Accuracy Comparison")
    print("=" * 70)

    # Load eskaks results
    eskaks_nei = parse_eskaks(os.path.join(results_dir, "eskaks_nei_small_pairwise_results.tsv"))
    eskaks_li = parse_eskaks(os.path.join(results_dir, "eskaks_li_small_pairwise_results.tsv"))

    # Load KaKs_Calculator results
    kakscalc_ng = parse_kakscalc(os.path.join(results_dir, "kakscalc_ng_small.tsv"))
    kakscalc_lpb = parse_kakscalc(os.path.join(results_dir, "kakscalc_lpb_small.tsv"))

    print(f"\n  eskaks Nei: {len(eskaks_nei)} pairs")
    print(f"  eskaks Li:  {len(eskaks_li)} pairs")
    print(f"  KaKs_Calc NG:  {len(kakscalc_ng)} pairs")
    print(f"  KaKs_Calc LPB: {len(kakscalc_lpb)} pairs")

    # Nei-Gojobori comparisons
    print("\n" + "-" * 70)
    print(" Nei-Gojobori Model Comparisons")
    print("-" * 70)
    compare("eskaks Nei", eskaks_nei, "KaKs_Calculator NG", kakscalc_ng)

    # Check for BioPython output
    bp_path = os.path.join(results_dir, "biopython_ng_small.tsv")
    if os.path.exists(bp_path):
        biopython_ng = parse_biopython(bp_path)
        print(f"\n  BioPython NG: {len(biopython_ng)} pairs")
        compare("eskaks Nei", eskaks_nei, "BioPython NG86", biopython_ng)
        compare("KaKs_Calculator NG", kakscalc_ng, "BioPython NG86", biopython_ng)

    # Li model comparisons
    print("\n" + "-" * 70)
    print(" Li (1993) Model Comparisons")
    print("-" * 70)
    compare("eskaks Li", eskaks_li, "KaKs_Calculator LPB", kakscalc_lpb)

    print("\n" + "=" * 70)


if __name__ == '__main__':
    main()
