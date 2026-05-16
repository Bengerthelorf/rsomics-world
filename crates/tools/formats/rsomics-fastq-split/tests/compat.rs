//! Byte-compat of `--split_by_lines` vs the pinned fastp 0.20.1 oracle
//! (deterministic both sides); skipped, not asserted, on a non-0.20 fastp.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-split"))
}

fn fastp_reference() -> Option<bool> {
    let out = Command::new("fastp")
        .arg("--version")
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    let v = String::from_utf8_lossy(&out.stderr) + String::from_utf8_lossy(&out.stdout);
    Some(v.contains("0.20"))
}

fn require_reference_fastp() -> bool {
    match fastp_reference() {
        None => {
            eprintln!("SKIP: fastp not on PATH — compat oracle unavailable");
            false
        }
        Some(false) => {
            eprintln!(
                "SKIP: local fastp is not the 0.20 compat reference (split file \
                 naming/distribution is version-specific); authoritative on 4090/CI"
            );
            false
        }
        Some(true) => true,
    }
}

fn run(bin: &Path, args: &[&str]) {
    let out = Command::new(bin).args(args).output().expect("spawn");
    assert!(
        out.status.success(),
        "{} {:?} failed:\nstdout: {}\nstderr: {}",
        bin.display(),
        args,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
}

/// The sorted `(file_name, bytes)` of every split file in `dir` whose name
/// ends with `.<base>` (fastp split naming is `<digits>.<base>`).
fn split_set(dir: &Path, base: &str) -> Vec<(String, Vec<u8>)> {
    let suffix = format!(".{base}");
    let mut v: Vec<(String, Vec<u8>)> = std::fs::read_dir(dir)
        .expect("read split dir")
        .filter_map(|e| {
            let p = e.ok()?.path();
            let name = p.file_name()?.to_string_lossy().into_owned();
            if name.ends_with(&suffix) && name != base {
                Some((name, std::fs::read(&p).expect("read split file")))
            } else {
                None
            }
        })
        .collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

#[test]
fn se_split_by_lines_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let od = tmp.path().join("ours");
    let td = tmp.path().join("theirs");
    std::fs::create_dir_all(&od).unwrap();
    std::fs::create_dir_all(&td).unwrap();
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let input = fixture("se5.fastq");

    run(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            od.join("o.fq").to_str().unwrap(),
            "--split_by_lines",
            "8",
        ],
    );
    run(
        Path::new("fastp"),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            td.join("o.fq").to_str().unwrap(),
            "--split_by_lines",
            "8",
            "--disable_adapter_trimming",
            "--disable_quality_filtering",
            "--disable_length_filtering",
            "--disable_trim_poly_g",
            "-j",
            j.to_str().unwrap(),
            "-h",
            h.to_str().unwrap(),
        ],
    );

    assert_eq!(
        split_set(&od, "o.fq"),
        split_set(&td, "o.fq"),
        "SE --split_by_lines: split file set (names+bytes) diverges from fastp 0.20.1",
    );
}

#[test]
fn pe_split_by_lines_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let od = tmp.path().join("ours");
    let td = tmp.path().join("theirs");
    std::fs::create_dir_all(&od).unwrap();
    std::fs::create_dir_all(&td).unwrap();
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let in1 = fixture("pe.fastq.r1");
    let in2 = fixture("pe.fastq.r2");

    run(
        &ours(),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            od.join("o1.fq").to_str().unwrap(),
            "-O",
            od.join("o2.fq").to_str().unwrap(),
            "--split_by_lines",
            "8",
        ],
    );
    run(
        Path::new("fastp"),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            td.join("o1.fq").to_str().unwrap(),
            "-O",
            td.join("o2.fq").to_str().unwrap(),
            "--split_by_lines",
            "8",
            "--disable_adapter_trimming",
            "--disable_quality_filtering",
            "--disable_length_filtering",
            "--disable_trim_poly_g",
            "-j",
            j.to_str().unwrap(),
            "-h",
            h.to_str().unwrap(),
        ],
    );

    assert_eq!(
        split_set(&od, "o1.fq"),
        split_set(&td, "o1.fq"),
        "PE R1 --split_by_lines: diverges from fastp 0.20.1",
    );
    assert_eq!(
        split_set(&od, "o2.fq"),
        split_set(&td, "o2.fq"),
        "PE R2 --split_by_lines: diverges from fastp 0.20.1",
    );
}
