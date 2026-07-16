//! Shared text-formatting helpers for output writers.
//!
//! Centralised so every output path (FASTA pairwise/lineage/window, VCF pN/pS and
//! McDonald-Kreitman) escapes and quotes identically — the alternative, per-writer
//! copies, is exactly what let some paths drift out of sync (unquoted ids corrupting
//! CSV columns, un-escaped gene names breaking JSON).

/// Escape a string for embedding inside a JSON double-quoted string, so an id or
/// gene name containing `"`, `\`, or a control char cannot produce invalid JSON.
pub(crate) fn json_escape(s: &str) -> String {
    let mut o = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            '\r' => o.push_str("\\r"),
            '\t' => o.push_str("\\t"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o
}

/// Quote a field for a delimited (CSV/TSV) file when it contains the separator, a
/// double quote, or a newline, so an id or gene name carrying the delimiter cannot
/// shift the columns of the row. RFC-4180 style: wrap in double quotes and double any
/// interior quote. (For TSV this keeps a fixed column count for RFC-4180-aware
/// readers; a bare embedded tab would otherwise add a spurious field.)
pub(crate) fn delim_field(s: &str, sep: char) -> String {
    if s.contains(sep) || s.contains('"') || s.contains('\n') || s.contains('\r') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Normalise IEEE negative zero to positive zero so the delimited writers never emit
/// a nonsensical "-0.000000" and stay byte-consistent with the JSON formatter (which
/// already maps -0.0 → 0.0).
pub(crate) fn norm_zero(v: f64) -> f64 {
    if v == 0.0 {
        0.0
    } else {
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn json_escape_neutralizes_quotes_backslash_controls() {
        assert_eq!(json_escape(r#"a"b\c"#), r#"a\"b\\c"#);
        assert_eq!(json_escape("a\tb\nc"), "a\\tb\\nc");
        assert_eq!(json_escape("x\u{0001}y"), "x\\u0001y");
    }

    #[test]
    fn delim_field_quotes_only_when_needed() {
        // Plain fields pass through unchanged.
        assert_eq!(delim_field("seqA", ','), "seqA");
        assert_eq!(delim_field("seqA", '\t'), "seqA");
        // A comma triggers quoting in CSV but not in TSV.
        assert_eq!(delim_field("NP_001,alpha", ','), "\"NP_001,alpha\"");
        assert_eq!(delim_field("NP_001,alpha", '\t'), "NP_001,alpha");
        // A tab triggers quoting in TSV.
        assert_eq!(delim_field("a\tb", '\t'), "\"a\tb\"");
        // Interior quotes are doubled; newlines force quoting.
        assert_eq!(delim_field("a\"b", ','), "\"a\"\"b\"");
        assert_eq!(delim_field("a\nb", ','), "\"a\nb\"");
    }

    #[test]
    fn norm_zero_collapses_negative_zero() {
        assert_eq!(norm_zero(-0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(norm_zero(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(norm_zero(1.5), 1.5);
        assert!(norm_zero(f64::NAN).is_nan());
    }
}
