pub mod li;
pub mod nei;

use clap::ValueEnum;

/// Valor Z para intervalo de confianza del 95% (aproximacion para muestras grandes).
pub const Z_95_CONFIDENCE: f64 = 1.96;

#[derive(ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum Model {
    Nei,
    Li,
}

#[derive(Clone, Copy, Debug)]
pub struct DsDn {
    pub dn: f64,
    pub ds: f64,
}
