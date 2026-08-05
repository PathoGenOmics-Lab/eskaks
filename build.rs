//! Embed the short git commit into the version string, so `eskaks --version` reports
//! exactly which build is running. Falls back to the bare crate version when `.git`
//! is absent (a crates.io tarball, a release archive, a shallow CI checkout).
use std::process::Command;

fn main() {
    let pkg = std::env::var("CARGO_PKG_VERSION").unwrap_or_default();

    let git = |args: &[&str]| {
        Command::new("git")
            .args(args)
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .filter(|s| !s.is_empty())
    };

    let version = match git(&["rev-parse", "--short", "HEAD"]) {
        Some(hash) => {
            let dirty = git(&["status", "--porcelain"]).is_some();
            let mark = if dirty { "-dirty" } else { "" };
            format!("{pkg} ({hash}{mark})")
        }
        None => pkg,
    };

    println!("cargo:rustc-env=ESKAKS_VERSION={version}");
    // Recompute the stamp when the checked-out commit or branch changes.
    println!("cargo:rerun-if-changed=.git/HEAD");
    println!("cargo:rerun-if-changed=.git/logs/HEAD");
}
