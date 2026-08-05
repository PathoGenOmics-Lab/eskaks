use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

mod accum;
mod dist;
mod diversity;
mod rng;
#[cfg(test)]
mod tests;

pub use accum::{FloatAccum, SummaryStats, WindowStats};
pub use diversity::{tajimas_d, theta_pi_varn, theta_watterson};
pub use dist::{
    benjamini_hochberg, binomial_two_sided_neglog10p, binomial_two_sided_p, bonferroni,
    chi2_from_two_sided_neglog10p, fisher_exact_two_sided, normal_two_sided_p, percentile_sorted,
    wilson_interval,
};
pub use rng::SplitMix64;
// Used only by the test suite.
#[cfg(test)]
pub(crate) use dist::{erfc, inv_normal_cdf, ln_gamma};
