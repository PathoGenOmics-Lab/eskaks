# eskaks — build, test, and benchmark targets
CARGO := cargo
RUSTFLAGS_NATIVE := RUSTFLAGS="-C target-cpu=native"

.PHONY: all build release test clippy clean benchmark bench-generate bench-run bench-plot docs docs-serve check

# ─── Build ────────────────────────────────────────────────────────────────────

all: release

build:
	$(CARGO) build

release:
	$(RUSTFLAGS_NATIVE) $(CARGO) build --release

# ─── Quality ──────────────────────────────────────────────────────────────────

test:
	$(CARGO) test

clippy:
	$(CARGO) clippy -- -D warnings

check: clippy test

# ─── Benchmarks ───────────────────────────────────────────────────────────────
# Requirements: python3, KaKs_Calculator (optional), BioPython (optional)

BENCH_DIR := benchmark

## Generate synthetic datasets for benchmarking
bench-generate:
	python3 $(BENCH_DIR)/generate_seqs.py

## Run cross-tool benchmarks (requires KaKs_Calculator and BioPython)
bench-run: release bench-generate
	python3 $(BENCH_DIR)/cross_tool_benchmark.py

## Generate accuracy/performance plots from results
bench-plot:
	python3 $(BENCH_DIR)/compare_accuracy.py

## Full benchmark pipeline: generate → run → plot
benchmark: bench-run bench-plot
	@echo "Benchmark complete. Results in $(BENCH_DIR)/cross_tool_results.json"
	@echo "Plots in $(BENCH_DIR)/plots/"

# ─── Documentation ────────────────────────────────────────────────────────────

docs:
	mdbook build docs

docs-serve:
	mdbook serve docs --open

# ─── Clean ────────────────────────────────────────────────────────────────────

clean:
	$(CARGO) clean
	rm -rf docs/book
