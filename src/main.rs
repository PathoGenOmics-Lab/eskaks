mod cli;
mod codon;
mod compute;
mod genetic_code;
mod gff;
mod input;
mod models;
mod output;
mod plot;
mod report;
mod stats;
mod textfmt;
mod tree;
mod vcf;
mod vcf_analysis;
mod run_fasta;
mod run_vcf;

use clap::{CommandFactory, Parser};
use log::LevelFilter;
use cli::{Args, SubCmd};

/// Bundled example alignment, embedded so `--demo` works from an installed binary
/// without the repository's `examples/` directory on disk.
const DEMO_FASTA: &str = include_str!("../examples/genes.fasta");

/// Bundled toy genome (reference + annotation + a small 20-sample genotyped VCF, plus
/// a divergence table), embedded so `eskaks --demo` can also showcase the per-gene
/// pN/pS path - the neutrality test, `--variants`, `--diversity` and the report -
/// with no files on disk. The VCF is genotyped so π / θ_W / Tajima's D are defined.
const DEMO_VCF_REF: &str = include_str!("../examples/toy_genome/reference.fasta");
const DEMO_VCF_GFF: &str = include_str!("../examples/toy_genome/genes.gff3");
const DEMO_VCF: &str = include_str!("../examples/toy_genome/variants_multisample.vcf");
const DEMO_VCF_DIVERGENCE: &str = include_str!("../examples/toy_genome/divergence.tsv");

/// Initialise rayon's global thread pool exactly once per process. A normal run touches
/// only one subcommand, but `--demo` runs `eskaks fasta` and `eskaks vcf` back to back
/// in a single process, so the second initialisation must be a no-op instead of the
/// hard "global thread pool has already been initialized" error.
fn init_global_pool(workers: usize) {
    static POOL: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL.get_or_init(|| {
        // build_global() fails only when a global pool already exists (in which case
        // rayon is still fully usable), so a failure here is benign and must not abort
        // the run. The OnceLock also means we attempt it at most once per process.
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .stack_size(4 * 1024 * 1024)
            .build_global();
    });
}

fn main() -> anyhow::Result<()> {
    // A biologist's intuitive first try is `eskaks alignment.fasta` (no subcommand). clap
    // reports "unrecognized subcommand '<file>'"; add a hint pointing at the right form
    // when the offending value looks like a file path, instead of just the bare error.
    let args = match Args::try_parse() {
        Ok(a) => a,
        Err(e) => {
            if e.kind() == clap::error::ErrorKind::InvalidSubcommand {
                if let Some(clap::error::ContextValue::String(val)) =
                    e.get(clap::error::ContextKind::InvalidSubcommand)
                {
                    if val.contains('.') || val.contains('/') {
                        eprintln!(
                            "error: '{val}' is not a subcommand — did you forget the mode?\n  \
                             FASTA alignment:  eskaks fasta {val} -o results\n  \
                             variants (VCF):   eskaks vcf --ref ref.fa --gff genes.gff3 --vcf {val}\n\
                             \nRun `eskaks --help` to see the two subcommands."
                        );
                        std::process::exit(2);
                    }
                }
            }
            e.exit();
        }
    };

    // Show data-quality warnings by default (previously env_logger defaulted to
    // "off", hiding every REF-mismatch / skipped-gene / saturation diagnostic
    // unless the user knew to set RUST_LOG). RUST_LOG still overrides this.
    let level = if args.quiet {
        LevelFilter::Error
    } else {
        match args.verbose {
            0 => LevelFilter::Warn,
            1 => LevelFilter::Info,
            _ => LevelFilter::Debug,
        }
    };
    // Clean, cargo-style log lines ("warning: ...", "error: ...") instead of
    // env_logger's default "[<ISO timestamp> LEVEL module::path] ..." — the timestamp
    // and module path are noise for an interactive scientific CLI. Colour is applied
    // only when stderr is a terminal (env_logger's auto write-style). Info-level lines
    // (the -v narration) print unadorned. RUST_LOG still overrides the level.
    env_logger::Builder::new()
        .filter_level(level)
        .format(|buf, record| {
            use std::io::Write;
            let lvl = record.level();
            let label = match lvl {
                log::Level::Error => "error",
                log::Level::Warn => "warning",
                log::Level::Info => "",
                log::Level::Debug => "debug",
                log::Level::Trace => "trace",
            };
            if label.is_empty() {
                writeln!(buf, "{}", record.args())
            } else {
                // env_logger 0.11 returns an anstyle Style rather than the old
                // wrapper whose `.value()` painted a single value. A Style now
                // renders its own escape sequences through Display: `{style}`
                // opens it and `{style:#}` closes it, so the label is wrapped
                // rather than handed to the style.
                let style = buf.default_level_style(lvl);
                writeln!(buf, "{style}{label}{style:#}: {}", record.args())
            }
        })
        .parse_default_env()
        .init();

    // Generate a shell completion script and exit (top-level flag, no subcommand needed).
    if let Some(shell) = args.completions {
        let mut cmd = Args::command();
        clap_complete::generate(shell, &mut cmd, "eskaks", &mut std::io::stdout());
        return Ok(());
    }

    // Handle --list-codes (top-level flag)
    if args.list_codes {
        eprintln!("Available NCBI genetic code tables:");
        for (id, name) in genetic_code::list_tables() {
            eprintln!("  {:>2}  {}", id, name);
        }
        eprintln!("\nApply one with --genetic-code <N>, e.g. `eskaks fasta aln.fasta --genetic-code 11`.");
        return Ok(());
    }

    // Quick demo on the bundled example data (top-level flag, no input files needed).
    if args.demo {
        return run_demo();
    }

    match args.command {
        Some(SubCmd::Fasta(fasta_args)) => run_fasta::run_fasta(fasta_args),
        Some(SubCmd::Vcf(vcf_args)) => run_vcf::run_vcf(vcf_args),
        None => {
            // No subcommand: print help
            Args::command().print_help()?;
            println!();
            Ok(())
        }
    }
}

/// Run both `eskaks fasta` and `eskaks vcf` on the embedded example data so a brand-new
/// user gets two real, successful runs (each with a summary and an HTML report) without
/// supplying any files, then point them at the real commands for their own data. The
/// VCF stage exercises the flagship per-gene pN/pS path - including `--variants` and
/// `--diversity` - which the FASTA-only demo used to leave unreachable on bundled data.
fn run_demo() -> anyhow::Result<()> {
    use anyhow::Context;
    let dir = std::env::temp_dir();

    // ── Demo 1/2: pairwise dN/dS from a codon-aligned FASTA ────────────────────
    let fasta = dir.join("eskaks_demo.fasta");
    let fa_out = dir.join("eskaks_demo");
    std::fs::write(&fasta, DEMO_FASTA)
        .with_context(|| format!("Failed to write demo input to {}", fasta.display()))?;

    eprintln!("Demo 1/2: `eskaks fasta` - pairwise dN/dS on a bundled 6-strain alignment.\n");
    let fasta_args = cli::FastaArgs::parse_from([
        "eskaks",
        fasta.to_str().unwrap(),
        "-o",
        fa_out.to_str().unwrap(),
        "--summary",
        "--report",
    ]);
    run_fasta::run_fasta(fasta_args)?;

    // ── Demo 2/2: per-gene pN/pS + diversity from a VCF + reference + GFF3 ──────
    let ref_p = dir.join("eskaks_demo_ref.fasta");
    let gff_p = dir.join("eskaks_demo_genes.gff3");
    let vcf_p = dir.join("eskaks_demo_variants.vcf");
    let div_p = dir.join("eskaks_demo_divergence.tsv");
    let vcf_out = dir.join("eskaks_demo_vcf");
    std::fs::write(&ref_p, DEMO_VCF_REF)
        .with_context(|| format!("Failed to write demo reference to {}", ref_p.display()))?;
    std::fs::write(&gff_p, DEMO_VCF_GFF)
        .with_context(|| format!("Failed to write demo annotation to {}", gff_p.display()))?;
    std::fs::write(&vcf_p, DEMO_VCF)
        .with_context(|| format!("Failed to write demo VCF to {}", vcf_p.display()))?;
    std::fs::write(&div_p, DEMO_VCF_DIVERGENCE)
        .with_context(|| format!("Failed to write demo divergence to {}", div_p.display()))?;

    eprintln!(
        "\nDemo 2/2: `eskaks vcf` - per-gene pN/pS, neutrality test, variants and diversity\n  \
         on a bundled 20-sample toy genome (bacterial genetic code 11).\n"
    );
    let vcf_args = cli::VcfArgs::parse_from([
        "eskaks",
        "--ref",
        ref_p.to_str().unwrap(),
        "--gff",
        gff_p.to_str().unwrap(),
        "--vcf",
        vcf_p.to_str().unwrap(),
        "--genetic-code",
        "11",
        "-o",
        vcf_out.to_str().unwrap(),
        "--summary",
        "--report",
        "--variants",
        "--diversity",
        "--mk",
        "--divergence",
        div_p.to_str().unwrap(),
    ]);
    run_vcf::run_vcf(vcf_args)?;

    eprintln!(
        "\nThat was a two-part demo on bundled data (open either _report.html in a browser).\n\
         To analyse your OWN data:\n  \
         pairwise dN/dS:  eskaks fasta your_alignment.fasta -o results\n  \
         per-gene pN/pS:  eskaks vcf --ref ref.fa --gff genes.gff3 --vcf calls.vcf -o results\n\
         \nDocs: https://pathogenomics-lab.github.io/eskaks/   ·   `eskaks --help` for all options."
    );
    Ok(())
}

