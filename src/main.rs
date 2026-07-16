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
mod vcf;
mod vcf_analysis;
mod run_fasta;
mod run_vcf;

use clap::Parser;
use log::LevelFilter;
use cli::{Args, SubCmd};

fn main() -> anyhow::Result<()> {
    let args = Args::parse();

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
    env_logger::Builder::new()
        .filter_level(level)
        .parse_default_env()
        .init();

    // Handle --list-codes (top-level flag)
    if args.list_codes {
        eprintln!("Available NCBI genetic code tables:");
        for (id, name) in genetic_code::list_tables() {
            eprintln!("  {:>2}  {}", id, name);
        }
        return Ok(());
    }

    match args.command {
        Some(SubCmd::Fasta(fasta_args)) => run_fasta::run_fasta(fasta_args),
        Some(SubCmd::Vcf(vcf_args)) => run_vcf::run_vcf(vcf_args),
        None => {
            // No subcommand: print help
            use clap::CommandFactory;
            Args::command().print_help()?;
            println!();
            Ok(())
        }
    }
}

