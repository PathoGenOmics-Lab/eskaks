#![no_main]
//! Fuzz the GFF3 parser: on ANY bytes it must return a clean Ok or Err, never panic.
use libfuzzer_sys::fuzz_target;
use std::io::Write;

fuzz_target!(|data: &[u8]| {
    let mut f = match tempfile::NamedTempFile::new() {
        Ok(f) => f,
        Err(_) => return,
    };
    if f.write_all(data).is_err() || f.flush().is_err() {
        return;
    }
    let _ = eskaks::gff::parse_gff3(f.path());
});
