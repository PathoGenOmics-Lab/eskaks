# Performance & Accuracy

eskaks implements the classical substitution models with precomputed lookup tables,
which makes it dramatically faster than the established tools while staying
numerically accurate.

## Speed

| Dataset | eskaks (4t) | KaKs_Calculator | PAML yn00 | BioPython | Speedup |
|---|---|---|---|---|---|
| 20 seq × 300 bp | 2 ms | 34 ms | 8 ms | 610 ms | 17× |
| 100 seq × 3 kb | 6 ms | 7,703 ms | 697 ms | 111,619 ms | **1,280×** |
| 500 seq × 3 kb | 74 ms | 195,456 ms | - | - | **2,641×** |

Output is deterministic regardless of the number of `--workers` threads.

## Accuracy

The Li model achieves **R² = 1.0** against KaKs_Calculator's LPB implementation.
Full accuracy data and the benchmarking methodology are in
[benchmarks/](https://github.com/PathoGenOmics-Lab/eskaks/tree/main/benchmarks).

## Feature comparison

| | eskaks | KaKs_Calculator | BioPython | PAML yn00 |
|---|---|---|---|---|
| Nei-Gojobori model | ✅ | ✅ | ✅ | ✅ |
| Li/LPB93 model | ✅ | ✅ | ❌ | ❌ |
| Per-gene pN/pS from VCF | ✅ | ❌ | ❌ | ❌ |
| Neutrality test + FDR | ✅ | ❌ | ❌ | ❌ |
| Interactive HTML report | ✅ | ❌ | ❌ | ❌ |
| Custom genetic codes | ✅ (20 tables) | ❌ | ❌ | Limited |
| JSON output / stdin pipe | ✅ | ❌ | ❌ | ❌ |
| Parallel | ✅ | ❌ | ❌ | ❌ |
| Speed (100 seq) | **6 ms** | 7,703 ms | 111,619 ms | 697 ms |
