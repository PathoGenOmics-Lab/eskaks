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
mod vcf;
mod vcf_analysis;
mod run_fasta;
mod run_vcf;

use clap::{CommandFactory, Parser};
use log::LevelFilter;
use cli::{Args, SubCmd};

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
                let style = buf.default_level_style(lvl);
                writeln!(buf, "{}: {}", style.value(label), record.args())
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

