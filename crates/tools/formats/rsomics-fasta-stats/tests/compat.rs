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
    Command::new("seqkit")
        .arg("--version")
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
    q1: u64,
    q2: u64,
    q3: u64,
    sum_gap: u64,
    n50: u64,
    n50_num: u64,
    gc_percent: f64,
    sum_n: u64,
}

fn parse_tabular(out: &str) -> Row {
    let mut lines = out.lines();
    let _header = lines.next().expect("header line");
    let data = lines.next().expect("data line");
    let cells: Vec<&str> = data.split('\t').collect();
    assert!(cells.len() == 8 || cells.len() == 16, "{cells:?}");
    let extended = if cells.len() == 16 {
        Some(ExtRow {
            q1: cells[8].parse().unwrap(),
            q2: cells[9].parse().unwrap(),
            q3: cells[10].parse().unwrap(),
            sum_gap: cells[11].parse().unwrap(),
            n50: cells[12].parse().unwrap(),
            n50_num: cells[13].parse().unwrap(),
            gc_percent: cells[14].parse().unwrap(),
            sum_n: cells[15].parse().unwrap(),
        })
    } else {
        None
    };
    Row {
        seq_type: cells[2].to_string(),
        num_seqs: cells[3].parse().unwrap(),
        sum_len: cells[4].parse().unwrap(),
        min_len: cells[5].parse().unwrap(),
        avg_len: cells[6].parse().unwrap(),
        max_len: cells[7].parse().unwrap(),
        extended,
    }
}

#[test]
fn tabular_basic_matches_seqkit() {
    if !seqkit_available() {
        eprintln!("SKIP: seqkit not on PATH; install via brew/conda for compat testing");
        return;
    }
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
    if !seqkit_available() {
        eprintln!("SKIP: seqkit not on PATH");
        return;
    }
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
    assert_eq!(ours_e.q1, theirs_e.q1, "Q1");
    assert_eq!(ours_e.q2, theirs_e.q2, "Q2");
    assert_eq!(ours_e.q3, theirs_e.q3, "Q3");
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
