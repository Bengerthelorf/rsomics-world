//! Byte-level compatibility check against upstream `seqkit stats --tabular`.
//!
//! Mechanism: run both binaries on the same golden FASTA, parse the
//! tabular output into a numeric per-column comparison. We do NOT diff
//! the raw bytes — the `file` column always differs (display path) and a
//! literal byte-diff would also pin us to seqkit's exact whitespace
//! handling. Numeric/categorical fields are what matter.

use std::path::PathBuf;
use std::process::{Command, Stdio};

const FIXTURE: &str = "tests/golden/tiny.fa";

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(FIXTURE)
}

fn rsomics_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fasta-stats"))
}

fn seqkit_available() -> bool {
    // seqkit uses `version` (subcommand), not `--version` (flag).
    Command::new("seqkit")
        .arg("version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn run_tabular(bin: &std::path::Path, args: &[&str]) -> String {
    let out = Command::new(bin)
        .args(args)
        .output()
        .expect("subprocess spawn");
    assert!(
        out.status.success(),
        "{} {args:?} failed: stderr=\n{}",
        bin.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf-8 stdout")
}

struct Row {
    seq_type: String,
    num_seqs: u64,
    sum_len: u64,
    min_len: u64,
    avg_len: f64,
    max_len: u64,
    extended: Option<ExtRow>,
}

struct ExtRow {
    // seqkit's `--all` may emit fractional quartiles for small inputs
    // (e.g. an even split), so parse as f64 even though our own output
    // rounds to %.0f for the tabular surface.
    q1: f64,
    q2: f64,
    q3: f64,
    sum_gap: u64,
    n50: u64,
    n50_num: u64,
    gc_percent: f64,
    sum_n: u64,
}

fn parse_tabular(out: &str) -> Row {
    // Header-driven parsing: seqkit's `--all` output for a FASTA input
    // includes the FASTQ-only columns Q20(%) / Q30(%) / AvgQual (all
    // zero), giving 19 cells instead of our 16. Index by name rather
    // than position so both layouts decode the same way.
    let mut lines = out.lines();
    let header = lines.next().expect("header line");
    let data = lines.next().expect("data line");
    let headers: Vec<&str> = header.split('\t').collect();
    let cells: Vec<&str> = data.split('\t').collect();
    let col = |name: &str| -> &str {
        let idx = headers
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("missing column {name} in {header:?}"));
        cells[idx]
    };
    let has_extended = headers.contains(&"Q1");
    let extended = if has_extended {
        Some(ExtRow {
            q1: col("Q1").parse().unwrap(),
            q2: col("Q2").parse().unwrap(),
            q3: col("Q3").parse().unwrap(),
            sum_gap: col("sum_gap").parse().unwrap(),
            n50: col("N50").parse().unwrap(),
            n50_num: col("N50_num").parse().unwrap(),
            gc_percent: col("GC(%)").parse().unwrap(),
            sum_n: col("sum_n").parse().unwrap(),
        })
    } else {
        None
    };
    Row {
        seq_type: col("type").to_string(),
        num_seqs: col("num_seqs").parse().unwrap(),
        sum_len: col("sum_len").parse().unwrap(),
        min_len: col("min_len").parse().unwrap(),
        avg_len: col("avg_len").parse().unwrap(),
        max_len: col("max_len").parse().unwrap(),
        extended,
    }
}

#[test]
fn tabular_basic_matches_seqkit() {
    assert!(
        seqkit_available(),
        "compat test requires seqkit on PATH (install via `brew install seqkit` / `apt install seqkit`)"
    );
    let fixture = fixture_path();
    let ours = parse_tabular(&run_tabular(
        &rsomics_bin(),
        &["--tabular", fixture.to_str().unwrap()],
    ));
    let theirs = parse_tabular(&run_tabular(
        std::path::Path::new("seqkit"),
        &["stats", "--tabular", fixture.to_str().unwrap()],
    ));
    assert_eq!(ours.seq_type, theirs.seq_type, "type");
    assert_eq!(ours.num_seqs, theirs.num_seqs, "num_seqs");
    assert_eq!(ours.sum_len, theirs.sum_len, "sum_len");
    assert_eq!(ours.min_len, theirs.min_len, "min_len");
    assert_eq!(ours.max_len, theirs.max_len, "max_len");
    assert!((ours.avg_len - theirs.avg_len).abs() < 0.05, "avg_len");
}

#[test]
fn tabular_all_matches_seqkit() {
    assert!(
        seqkit_available(),
        "compat test requires seqkit on PATH (install via `brew install seqkit` / `apt install seqkit`)"
    );
    let fixture = fixture_path();
    let ours = parse_tabular(&run_tabular(
        &rsomics_bin(),
        &["--tabular", "--all", fixture.to_str().unwrap()],
    ));
    let theirs = parse_tabular(&run_tabular(
        std::path::Path::new("seqkit"),
        &["stats", "--tabular", "--all", fixture.to_str().unwrap()],
    ));
    let ours_e = ours.extended.expect("our --all extended");
    let theirs_e = theirs.extended.expect("seqkit --all extended");
    assert!((ours_e.q1 - theirs_e.q1).abs() < 0.5, "Q1");
    assert!((ours_e.q2 - theirs_e.q2).abs() < 0.5, "Q2");
    assert!((ours_e.q3 - theirs_e.q3).abs() < 0.5, "Q3");
    assert_eq!(ours_e.sum_gap, theirs_e.sum_gap, "sum_gap");
    assert_eq!(ours_e.n50, theirs_e.n50, "N50");
    assert_eq!(ours_e.n50_num, theirs_e.n50_num, "N50_num");
    assert_eq!(ours_e.sum_n, theirs_e.sum_n, "sum_n");
    assert!(
        (ours_e.gc_percent - theirs_e.gc_percent).abs() < 0.02,
        "GC(%): {} vs {}",
        ours_e.gc_percent,
        theirs_e.gc_percent
    );
}
