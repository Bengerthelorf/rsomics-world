use std::path::PathBuf;
use std::process::{Command, Stdio};

const FIXTURE: &str = "tests/golden/tiny.fq";

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(FIXTURE)
}

fn rsomics_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-stats"))
}

fn seqkit_available() -> bool {
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
    q1: f64,
    q2: f64,
    q3: f64,
    sum_gap: u64,
    n50: u64,
    n50_num: u64,
    q20: f64,
    q30: f64,
    avg_qual: f64,
    gc_percent: f64,
    sum_n: u64,
}

fn parse_tabular(out: &str) -> Row {
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
            q20: col("Q20(%)").parse().unwrap(),
            q30: col("Q30(%)").parse().unwrap(),
            avg_qual: col("AvgQual").parse().unwrap(),
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
    let o = ours.extended.expect("our --all extended");
    let t = theirs.extended.expect("seqkit --all extended");
    assert!((o.q1 - t.q1).abs() < 0.5, "Q1 {} vs {}", o.q1, t.q1);
    assert!((o.q2 - t.q2).abs() < 0.5, "Q2 {} vs {}", o.q2, t.q2);
    assert!((o.q3 - t.q3).abs() < 0.5, "Q3 {} vs {}", o.q3, t.q3);
    assert_eq!(o.sum_gap, t.sum_gap, "sum_gap");
    assert_eq!(o.n50, t.n50, "N50");
    assert_eq!(o.n50_num, t.n50_num, "N50_num");
    assert_eq!(o.sum_n, t.sum_n, "sum_n");
    assert!(
        (o.q20 - t.q20).abs() < 1.0,
        "Q20(%): {} vs {}",
        o.q20,
        t.q20
    );
    assert!(
        (o.q30 - t.q30).abs() < 1.0,
        "Q30(%): {} vs {}",
        o.q30,
        t.q30
    );
    assert!(
        (o.avg_qual - t.avg_qual).abs() < 0.05,
        "AvgQual: {} vs {}",
        o.avg_qual,
        t.avg_qual
    );
    assert!(
        (o.gc_percent - t.gc_percent).abs() < 0.02,
        "GC(%): {} vs {}",
        o.gc_percent,
        t.gc_percent
    );
}

// Heterogeneous fixture (DNA first record, protein-looking later record).
// seqkit guesses `type` from the first record only — this would fail with a
// cross-record alphabet sample, so it pins the first-record-scope contract
// against the real binary (the homogeneous golden masks it).
#[test]
fn heterogeneous_type_is_first_record_like_seqkit() {
    assert!(
        seqkit_available(),
        "compat test requires seqkit on PATH (install via `brew install seqkit` / `apt install seqkit`)"
    );
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/hetero.fq");
    let ours = parse_tabular(&run_tabular(
        &rsomics_bin(),
        &["--tabular", "--all", fixture.to_str().unwrap()],
    ));
    let theirs = parse_tabular(&run_tabular(
        std::path::Path::new("seqkit"),
        &["stats", "--tabular", "--all", fixture.to_str().unwrap()],
    ));
    assert_eq!(ours.seq_type, theirs.seq_type, "type (first-record scope)");
    assert_eq!(ours.seq_type, "DNA", "first record is DNA");
    let o = ours.extended.expect("our --all extended");
    let t = theirs.extended.expect("seqkit --all extended");
    assert_eq!(o.sum_n, t.sum_n, "sum_n (type-dependent N vs X)");
    assert!(
        (o.gc_percent - t.gc_percent).abs() < 0.02,
        "GC(%): {} vs {}",
        o.gc_percent,
        t.gc_percent
    );
}
