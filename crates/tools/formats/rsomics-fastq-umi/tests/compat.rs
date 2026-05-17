use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-umi"))
}

// UMI read-name format and flag semantics differ across major fastp versions; oracle is pinned to 0.20.x (CI/4090 reference).
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
                "SKIP: local fastp is not the 0.20 compat reference (UMI format \
                 is version-specific); authoritative on 4090/CI fastp 0.20.1"
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

// fastp flags isolating the UMI transform (0.20.1 long forms):
//   --disable_adapter_trimming --disable_quality_filtering
//   --disable_length_filtering --disable_trim_poly_g
const FASTP_ISOLATE: &[&str] = &[
    "--disable_adapter_trimming",
    "--disable_quality_filtering",
    "--disable_length_filtering",
    "--disable_trim_poly_g",
];

#[test]
fn se_umi_read1_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let input = fixture("se_umi.fastq");

    run(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "--umi_len",
            "6",
        ],
    );
    let mut fp = vec![
        "-i",
        input.to_str().unwrap(),
        "-o",
        theirs_out.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "read1",
        "--umi_len",
        "6",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE UMI read1: byte-level FASTQ output diverges from fastp 0.20.1",
    );
}

#[test]
fn se_umi_with_skip_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let input = fixture("se_umi.fastq");

    run(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "--umi_len",
            "5",
            "--umi_skip",
            "2",
        ],
    );
    let mut fp = vec![
        "-i",
        input.to_str().unwrap(),
        "-o",
        theirs_out.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "read1",
        "--umi_len",
        "5",
        "--umi_skip",
        "2",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE UMI read1 + skip: byte-level output diverges from fastp 0.20.1",
    );
}

#[test]
fn pe_umi_read1_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    let t1 = tmp.path().join("t1.fq");
    let t2 = tmp.path().join("t2.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let in1 = fixture("pe_umi.fastq.r1");
    let in2 = fixture("pe_umi.fastq.r2");

    run(
        &ours(),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
            "--umi_len",
            "8",
        ],
    );
    let mut fp = vec![
        "-i",
        in1.to_str().unwrap(),
        "-I",
        in2.to_str().unwrap(),
        "-o",
        t1.to_str().unwrap(),
        "-O",
        t2.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "read1",
        "--umi_len",
        "8",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&o1).unwrap(),
        std::fs::read(&t1).unwrap(),
        "PE UMI R1: byte-level output diverges from fastp 0.20.1",
    );
    assert_eq!(
        std::fs::read(&o2).unwrap(),
        std::fs::read(&t2).unwrap(),
        "PE UMI R2: byte-level output diverges from fastp 0.20.1",
    );
}

// Byte-equality proves fastp 0.20.1 Read::trimFront clamp-to-length()-1 (keep ≥1 base) is matched for reads shorter than umi_len+skip.
#[test]
fn se_umi_short_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let input = fixture("se_umi_short.fastq");

    run(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "--umi_len",
            "8",
        ],
    );
    let mut fp = vec![
        "-i",
        input.to_str().unwrap(),
        "-o",
        theirs_out.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "read1",
        "--umi_len",
        "8",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE short-read UMI: trimFront clamp diverges from fastp 0.20.1",
    );
}

#[test]
fn se_umi_index1_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let ours_out = tmp.path().join("ours.fq");
    let theirs_out = tmp.path().join("theirs.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let input = fixture("se_umi.fastq");

    run(
        &ours(),
        &[
            "-i",
            input.to_str().unwrap(),
            "-o",
            ours_out.to_str().unwrap(),
            "--umi_loc",
            "index1",
        ],
    );
    let mut fp = vec![
        "-i",
        input.to_str().unwrap(),
        "-o",
        theirs_out.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "index1",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&ours_out).unwrap(),
        std::fs::read(&theirs_out).unwrap(),
        "SE UMI index1: read-name index extraction diverges from fastp 0.20.1",
    );
}

#[test]
fn pe_umi_index2_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    let t1 = tmp.path().join("t1.fq");
    let t2 = tmp.path().join("t2.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let in1 = fixture("pe_index.fastq.r1");
    let in2 = fixture("pe_index.fastq.r2");

    run(
        &ours(),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
            "--umi_loc",
            "index2",
        ],
    );
    let mut fp = vec![
        "-i",
        in1.to_str().unwrap(),
        "-I",
        in2.to_str().unwrap(),
        "-o",
        t1.to_str().unwrap(),
        "-O",
        t2.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "index2",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&o1).unwrap(),
        std::fs::read(&t1).unwrap(),
        "PE UMI index2 R1: diverges from fastp 0.20.1",
    );
    assert_eq!(
        std::fs::read(&o2).unwrap(),
        std::fs::read(&t2).unwrap(),
        "PE UMI index2 R2: diverges from fastp 0.20.1",
    );
}

#[test]
fn pe_umi_per_index_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    let t1 = tmp.path().join("t1.fq");
    let t2 = tmp.path().join("t2.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let in1 = fixture("pe_index.fastq.r1");
    let in2 = fixture("pe_index.fastq.r2");

    run(
        &ours(),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
            "--umi_loc",
            "per_index",
        ],
    );
    let mut fp = vec![
        "-i",
        in1.to_str().unwrap(),
        "-I",
        in2.to_str().unwrap(),
        "-o",
        t1.to_str().unwrap(),
        "-O",
        t2.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "per_index",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&o1).unwrap(),
        std::fs::read(&t1).unwrap(),
        "PE UMI per_index R1: firstIndex_lastIndex merge diverges from fastp 0.20.1",
    );
    assert_eq!(
        std::fs::read(&o2).unwrap(),
        std::fs::read(&t2).unwrap(),
        "PE UMI per_index R2: firstIndex_lastIndex merge diverges from fastp 0.20.1",
    );
}

#[test]
fn pe_umi_per_read_matches_fastp() {
    if !require_reference_fastp() {
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let o1 = tmp.path().join("o1.fq");
    let o2 = tmp.path().join("o2.fq");
    let t1 = tmp.path().join("t1.fq");
    let t2 = tmp.path().join("t2.fq");
    let j = tmp.path().join("fastp.json");
    let h = tmp.path().join("fastp.html");
    let in1 = fixture("pe_umi.fastq.r1");
    let in2 = fixture("pe_umi.fastq.r2");

    run(
        &ours(),
        &[
            "-i",
            in1.to_str().unwrap(),
            "-I",
            in2.to_str().unwrap(),
            "-o",
            o1.to_str().unwrap(),
            "-O",
            o2.to_str().unwrap(),
            "--umi_loc",
            "per_read",
            "--umi_len",
            "6",
        ],
    );
    let mut fp = vec![
        "-i",
        in1.to_str().unwrap(),
        "-I",
        in2.to_str().unwrap(),
        "-o",
        t1.to_str().unwrap(),
        "-O",
        t2.to_str().unwrap(),
        "--umi",
        "--umi_loc",
        "per_read",
        "--umi_len",
        "6",
    ];
    fp.extend_from_slice(FASTP_ISOLATE);
    fp.extend_from_slice(&["-j", j.to_str().unwrap(), "-h", h.to_str().unwrap()]);
    run(Path::new("fastp"), &fp);

    assert_eq!(
        std::fs::read(&o1).unwrap(),
        std::fs::read(&t1).unwrap(),
        "PE UMI per_read R1: umi1_umi2 merge + trim diverges from fastp 0.20.1",
    );
    assert_eq!(
        std::fs::read(&o2).unwrap(),
        std::fs::read(&t2).unwrap(),
        "PE UMI per_read R2: umi1_umi2 merge + trim diverges from fastp 0.20.1",
    );
}
