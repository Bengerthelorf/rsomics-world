use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-filter"))
}

fn fastp_available() -> bool {
    Command::new("fastp")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn run_to_path(bin: &Path, args: &[&str]) {
    let out = Command::new(bin)
        .args(args)
        .output()
        .expect("subprocess spawn");
    assert!(
        out.status.success(),
        "{} {:?} failed:\nstdout: {}\nstderr: {}",
        bin.display(),
        args,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
}

fn read_ids(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .expect("read fastq")
        .lines()
        .filter(|l| l.starts_with('@'))
        .map(|l| l[1..].split_whitespace().next().unwrap_or("").to_owned())
        .collect()
}

// fastp flags used to isolate quality+length filter, disabling everything else:
//   -A  disable adapter trimming
//   -G  disable poly-G trimming
//   --disable_correction  (PE only; PE overlap correction off)
//   --qualified_quality_phred N
//   --unqualified_percent_limit U
//   --length_required M
//   --length_limit L (0 = omit / not passed since fastp omits it when 0)
//   -j / -h  required fastp report sinks

#[test]
fn se_quality_filter_matches_fastp_defaults() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_quality.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-L",
        ],
    );
    run_to_path(
        Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-G",
            "--disable_length_filtering",
            "--qualified_quality_phred",
            "15",
            "--unqualified_percent_limit",
            "40",
            "--n_base_limit",
            "5",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );

    let ours_ids = read_ids(&ours_out);
    let theirs_ids = read_ids(&theirs_out);
    assert_eq!(
        ours_ids, theirs_ids,
        "SE quality filter: read set diverges from fastp"
    );
    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE quality filter: byte-level FASTQ output diverges from fastp",
    );
}

#[test]
fn se_quality_filter_matches_fastp_strict() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_quality.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-L",
            "--qualified_quality_phred",
            "20",
            "--unqualified_percent_limit",
            "20",
            "--n_base_limit",
            "2",
        ],
    );
    run_to_path(
        Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-G",
            "--disable_length_filtering",
            "--qualified_quality_phred",
            "20",
            "--unqualified_percent_limit",
            "20",
            "--n_base_limit",
            "2",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );

    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE strict quality filter: byte-level FASTQ output diverges from fastp",
    );
}

#[test]
fn se_length_filter_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_length.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-Q",
            "-l",
            "15",
        ],
    );
    run_to_path(
        Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-G",
            "--disable_quality_filtering",
            "--length_required",
            "15",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );

    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE length filter: byte-level FASTQ output diverges from fastp",
    );
}

fn write_pe_record(w: &mut impl Write, id: &str, seq: &[u8], qual_byte: u8) {
    writeln!(w, "@{id}").unwrap();
    w.write_all(seq).unwrap();
    writeln!(w).unwrap();
    writeln!(w, "+").unwrap();
    for _ in 0..seq.len() {
        w.write_all(&[qual_byte]).unwrap();
    }
    writeln!(w).unwrap();
}

fn make_pe_quality_fixture(in1: &Path, in2: &Path) {
    let high_qual: u8 = b'I'; // ASCII 73 = Q40; above any threshold we test
    let low_qual: u8 = b'!'; // ASCII 33 = Q0; below any threshold we test

    let mut w1 = std::fs::File::create(in1).unwrap();
    let mut w2 = std::fs::File::create(in2).unwrap();

    // pair1: both high-qual → should survive
    write_pe_record(&mut w1, "pair1", b"ACGTACGTACGTACGTACGT", high_qual);
    write_pe_record(&mut w2, "pair1", b"TGCATGCATGCATGCATGCA", high_qual);

    // pair2: R1 has 10 low-qual out of 20 (>40% → fail), R2 high-qual → whole pair dropped
    let mut r1_mixed = vec![low_qual; 10];
    r1_mixed.extend_from_slice(&[high_qual; 10]);
    let mut w1_buf = std::io::BufWriter::new(&mut w1);
    writeln!(w1_buf, "@pair2").unwrap();
    w1_buf.write_all(b"ACGTACGTACGTACGTACGT").unwrap();
    writeln!(w1_buf, "\n+").unwrap();
    w1_buf.write_all(&r1_mixed).unwrap();
    writeln!(w1_buf).unwrap();
    drop(w1_buf);
    write_pe_record(&mut w2, "pair2", b"TGCATGCATGCATGCATGCA", high_qual);

    // pair3: R1 high-qual, R2 has 10 low-qual → whole pair dropped
    write_pe_record(&mut w1, "pair3", b"ACGTACGTACGTACGTACGT", high_qual);
    let mut r2_mixed = vec![low_qual; 10];
    r2_mixed.extend_from_slice(&[high_qual; 10]);
    let mut w2_buf = std::io::BufWriter::new(&mut w2);
    writeln!(w2_buf, "@pair3").unwrap();
    w2_buf.write_all(b"TGCATGCATGCATGCATGCA").unwrap();
    writeln!(w2_buf, "\n+").unwrap();
    w2_buf.write_all(&r2_mixed).unwrap();
    writeln!(w2_buf).unwrap();
    drop(w2_buf);

    // pair4: both mates fail (all low-qual) → whole pair dropped
    write_pe_record(&mut w1, "pair4", b"ACGTACGTACGTACGTACGT", low_qual);
    write_pe_record(&mut w2, "pair4", b"TGCATGCATGCATGCATGCA", low_qual);
}

#[test]
fn pe_quality_filter_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let in1 = tmp.path().join("pe_r1.fq");
    let in2 = tmp.path().join("pe_r2.fq");
    make_pe_quality_fixture(&in1, &in2);

    let ours_out1 = tmp.path().join("ours_r1.fq");
    let ours_out2 = tmp.path().join("ours_r2.fq");
    let theirs_out1 = tmp.path().join("theirs_r1.fq");
    let theirs_out2 = tmp.path().join("theirs_r2.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");

    run_to_path(
        &ours(),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            ours_out1.to_str().unwrap(),
            "-O",
            ours_out2.to_str().unwrap(),
            "-L",
        ],
    );
    run_to_path(
        Path::new("fastp"),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            theirs_out1.to_str().unwrap(),
            "-O",
            theirs_out2.to_str().unwrap(),
            "-A",
            "-G",
            "--disable_length_filtering",
            "--qualified_quality_phred",
            "15",
            "--unqualified_percent_limit",
            "40",
            "--n_base_limit",
            "5",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );

    assert_eq!(
        std::fs::read(&ours_out1).unwrap(),
        std::fs::read(&theirs_out1).unwrap(),
        "PE quality filter R1: byte-level output diverges from fastp",
    );
    assert_eq!(
        std::fs::read(&ours_out2).unwrap(),
        std::fs::read(&theirs_out2).unwrap(),
        "PE quality filter R2: byte-level output diverges from fastp",
    );
}
