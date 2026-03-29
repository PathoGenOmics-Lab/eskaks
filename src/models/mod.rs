pub mod li;
pub mod nei;

use clap::ValueEnum;

/// Z-value for 95% confidence interval (large sample approximation).
pub const Z_95_CONFIDENCE: f64 = 1.96;

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    Nei,
    Li,
}

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    Tsv,
    Csv,
    Json,
}

impl OutputFormat {
    pub fn separator(self) -> char {
        match self {
            OutputFormat::Tsv => '\t',
            OutputFormat::Csv => ',',
            OutputFormat::Json => ',', // unused for JSON, but needed for the trait
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            OutputFormat::Tsv => "tsv",
            OutputFormat::Csv => "csv",
            OutputFormat::Json => "json",
        }
    }

}

#[derive(Clone, Copy, Debug)]
pub struct DsDn {
    pub dn: f64,
    pub ds: f64,
}
