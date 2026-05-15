use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-trim"))
}

fn fastp_available() -> bool {
    Command::new("fastp")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

fn run_to_path(bin: &std::path::Path, args: &[&str]) {
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
}

#[test]
fn se_adapter_trim_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_adapter.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-a",
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-a",
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "-Q",
            "-L",
            "-G",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours_bytes = std::fs::read(&ours_out).unwrap();
    let theirs_bytes = std::fs::read(&theirs_out).unwrap();
    assert_eq!(ours_bytes, theirs_bytes, "FASTQ output differs");
}

#[test]
fn se_polyg_trim_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_polyg.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-A",
            "-g",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-g",
            "-Q",
            "-L",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours_bytes = std::fs::read(&ours_out).unwrap();
    let theirs_bytes = std::fs::read(&theirs_out).unwrap();
    assert_eq!(ours_bytes, theirs_bytes, "polyG output differs");
}

#[test]
fn se_fixed_trim_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let fastp_json = tmp.path().join("fastp.json");
    let fastp_html = tmp.path().join("fastp.html");
    let input = fixture("se_adapter.fastq");

    run_to_path(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "-A",
            "-f",
            "4",
            "--trim_tail1",
            "2",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            theirs_out.to_str().unwrap(),
            "-A",
            "-f",
            "4",
            "-t",
            "2",
            "-Q",
            "-L",
            "-G",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours_bytes = std::fs::read(&ours_out).unwrap();
    let theirs_bytes = std::fs::read(&theirs_out).unwrap();
    assert_eq!(ours_bytes, theirs_bytes, "fixed-trim output differs");
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

fn write_pe_record(w: &mut impl Write, id: &str, seq: &[u8]) {
    writeln!(w, "@{id}").unwrap();
    w.write_all(seq).unwrap();
    w.write_all(b"\n+\n").unwrap();
    for _ in 0..seq.len() {
        w.write_all(b"I").unwrap();
    }
    w.write_all(b"\n").unwrap();
}

fn make_pe_overlap_fixture(in1: &Path, in2: &Path) {
    // Inserts ≥30 bp (fastp's overlap_len_require default), non-repetitive,
    // no homopolymer: tandem repeats cause near-matches; trailing G/C runs
    // interact with fastp's adapter heuristics in ways that diverge from ours.
    let inserts: &[&[u8]] = &[
        b"GCATATCAGTGCATATCAGTAATGCATGCAT",
        b"CGCGCATGCATGCATGCATTAGTCAGGACGT",
    ];
    let adapter1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
    let adapter2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";

    let mut w1 = std::fs::File::create(in1).unwrap();
    let mut w2 = std::fs::File::create(in2).unwrap();
    for (i, insert) in inserts.iter().enumerate() {
        let mut r1 = insert.to_vec();
        r1.extend_from_slice(adapter1);
        let rc = rev_comp(insert);
        let mut r2 = rc.clone();
        r2.extend_from_slice(adapter2);
        write_pe_record(&mut w1, &format!("pair_{i}"), &r1);
        write_pe_record(&mut w2, &format!("pair_{i}"), &r2);
    }
}

#[test]
fn pe_overlap_detect_matches_fastp() {
    assert!(
        fastp_available(),
        "compat test requires fastp on PATH (install via `brew install fastp` / `apt install fastp`)"
    );
    let tmp = tempfile::tempdir().unwrap();
    let in1 = tmp.path().join("pe_r1.fq");
    let in2 = tmp.path().join("pe_r2.fq");
    make_pe_overlap_fixture(&in1, &in2);

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
            "-2",
            "-L",
        ],
    );
    run_to_path(
        std::path::Path::new("fastp"),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            theirs_out1.to_str().unwrap(),
            "-O",
            theirs_out2.to_str().unwrap(),
            "--detect_adapter_for_pe",
            "-Q",
            "-L",
            "-G",
            "-j",
            fastp_json.to_str().unwrap(),
            "-h",
            fastp_html.to_str().unwrap(),
        ],
    );
    let ours1 = std::fs::read(&ours_out1).unwrap();
    let theirs1 = std::fs::read(&theirs_out1).unwrap();
    assert_eq!(ours1, theirs1, "PE R1 differs");
    let ours2 = std::fs::read(&ours_out2).unwrap();
    let theirs2 = std::fs::read(&theirs_out2).unwrap();
    assert_eq!(ours2, theirs2, "PE R2 differs");
}
