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
}

impl OutputFormat {
    pub fn separator(self) -> char {
        match self {
            OutputFormat::Tsv => '\t',
            OutputFormat::Csv => ',',
        }
    }
    pub fn extension(self) -> &'static str {
        match self {
            OutputFormat::Tsv => "tsv",
            OutputFormat::Csv => "csv",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct DsDn {
    pub dn: f64,
    pub ds: f64,
}
