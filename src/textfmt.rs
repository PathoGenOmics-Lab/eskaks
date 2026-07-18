//! Shared text-formatting helpers for output writers.
//!
//! Centralised so every output path (FASTA pairwise/lineage/window, VCF pN/pS and
//! McDonald-Kreitman) escapes and quotes identically — the alternative, per-writer
//! copies, is exactly what let some paths drift out of sync (unquoted ids corrupting
//! CSV columns, un-escaped gene names breaking JSON).

use std::borrow::Cow;

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
    // A bare double quote only needs escaping in CSV (RFC 4180). TSV has no quoting
    // convention, so a field like `gene"1` (a quote but no tab) must pass through
    // verbatim — wrapping it in quotes there would corrupt the value for TSV readers.
    let needs_quote =
        s.contains(sep) || s.contains('\n') || s.contains('\r') || (sep == ',' && s.contains('"'));
    if needs_quote {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Neutralise spreadsheet formula injection in a CSV *name* cell. A gene, sequence,
/// or chromosome name derived from user input can begin with `=`, `@`, or a `+`/`-`
/// followed by more than a lone sign - which Excel / LibreOffice / Sheets evaluate as
/// a formula when the CSV is opened. Prefix such a cell with a single quote so it is
/// read as literal text. Mirrors the interactive report's CSV export
/// (`src/report/script.js`): a lone `+`/`-` (e.g. a strand value) is left intact, and
/// numeric columns are never routed through here (a leading `-` on a negative number
/// would otherwise be mis-neutralised).
fn csv_formula_guard(s: &str) -> Cow<'_, str> {
    let b = s.as_bytes();
    let is_formula = matches!(b.first(), Some(b'=') | Some(b'@'))
        || (b.len() > 1 && matches!(b.first(), Some(b'+') | Some(b'-')));
    if is_formula {
        Cow::Owned(format!("'{s}"))
    } else {
        Cow::Borrowed(s)
    }
}

/// Quote/escape a user-derived NAME field for a delimited file, additionally defusing
/// spreadsheet formula injection when the delimiter is a comma (CSV). Use this for
/// gene / sequence / chromosome names; numeric columns format themselves and must not
/// pass through here. TSV is left un-guarded (its readers are programmatic and a
/// leading apostrophe would corrupt the value), matching the report's CSV-only scope.
pub(crate) fn name_field(s: &str, sep: char) -> String {
    if sep == ',' {
        delim_field(&csv_formula_guard(s), sep)
    } else {
        delim_field(s, sep)
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
        // Interior quotes are doubled in CSV; newlines force quoting (both formats).
        assert_eq!(delim_field("a\"b", ','), "\"a\"\"b\"");
        assert_eq!(delim_field("a\nb", ','), "\"a\nb\"");
        assert_eq!(delim_field("a\nb", '\t'), "\"a\nb\"");
        // A bare quote (no tab/newline) must NOT be quoted in TSV — it would corrupt the value.
        assert_eq!(delim_field("gene\"1", '\t'), "gene\"1");
    }

    #[test]
    fn name_field_defuses_csv_formula_injection() {
        // Formula triggers get a leading apostrophe in CSV...
        assert_eq!(name_field("=cmd()", ','), "'=cmd()");
        assert_eq!(name_field("@SUM(A1)", ','), "'@SUM(A1)");
        assert_eq!(name_field("-2+3", ','), "'-2+3");
        assert_eq!(name_field("+A1", ','), "'+A1");
        // ...but a lone +/- (e.g. a strand value) passes through intact,
        assert_eq!(name_field("-", ','), "-");
        assert_eq!(name_field("+", ','), "+");
        // an ordinary name is untouched,
        assert_eq!(name_field("katG", ','), "katG");
        // TSV is never formula-guarded (programmatic readers),
        assert_eq!(name_field("=cmd()", '\t'), "=cmd()");
        // and guarding composes with RFC-4180 quoting: a formula name that also
        // carries the delimiter is both apostrophe-prefixed and quoted.
        assert_eq!(name_field("=a,b", ','), "\"'=a,b\"");
    }

    #[test]
    fn norm_zero_collapses_negative_zero() {
        assert_eq!(norm_zero(-0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(norm_zero(0.0).to_bits(), 0.0f64.to_bits());
        assert_eq!(norm_zero(1.5), 1.5);
        assert!(norm_zero(f64::NAN).is_nan());
    }
}
