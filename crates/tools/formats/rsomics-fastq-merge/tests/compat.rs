use std::collections::BTreeMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const PINNED: &str = "0.20.1";
const N_PAIRS: usize = 500;
const FRAG: usize = 200;
const READ: usize = 150;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-merge"))
}

fn fastp_pinned() -> bool {
    let Ok(out) = Command::new("fastp").arg("--version").output() else {
        return false;
    };
    // fastp prints its version banner to stderr, not stdout.
    String::from_utf8_lossy(&out.stderr).contains(PINNED)
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

fn synth(r1p: &Path, r2p: &Path) {
    let mut w1 = BufWriter::new(fs::File::create(r1p).unwrap());
    let mut w2 = BufWriter::new(fs::File::create(r2p).unwrap());
    let mut rng = 0x1234_5678_9abc_def0u64;
    for i in 0..N_PAIRS {
        let mut frag = Vec::with_capacity(FRAG);
        for _ in 0..FRAG {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            frag.push(b"ACGT"[((rng >> 33) & 3) as usize]);
        }
        let r1 = &frag[..READ];
        let r2 = revcomp(&frag[FRAG - READ..]);
        let q = vec![b'I'; READ];
        writeln!(w1, "@frag{i}/1").unwrap();
        w1.write_all(r1).unwrap();
        w1.write_all(b"\n+\n").unwrap();
        w1.write_all(&q).unwrap();
        w1.write_all(b"\n").unwrap();
        writeln!(w2, "@frag{i}/2").unwrap();
        w2.write_all(&r2).unwrap();
        w2.write_all(b"\n+\n").unwrap();
        w2.write_all(&q).unwrap();
        w2.write_all(b"\n").unwrap();
    }
}

fn merged_seqs(fq: &str) -> BTreeMap<String, usize> {
    let mut m = BTreeMap::new();
    let mut lines = fq.lines();
    while let Some(h) = lines.next() {
        if !h.starts_with('@') {
            continue;
        }
        if let Some(seq) = lines.next() {
            *m.entry(seq.to_string()).or_insert(0) += 1;
        }
        lines.next(); // +
        lines.next(); // qual
    }
    m
}

#[test]
fn merged_set_matches_fastp() {
    if !fastp_pinned() {
        eprintln!(
            "SKIP: fastp v{PINNED} not on PATH — merge oracle unavailable \
             (newer fastp drifts; authoritative on the pinned CI lane)"
        );
        return;
    }
    let tmp = tempfile::tempdir().unwrap();
    let r1 = tmp.path().join("r1.fq");
    let r2 = tmp.path().join("r2.fq");
    synth(&r1, &r2);

    let fp_out = tmp.path().join("fastp_merged.fq");
    let ours_out = tmp.path().join("ours_merged.fq");

    let fp = Command::new("fastp")
        .args(["--merge", "--merged_out"])
        .arg(&fp_out)
        .args(["--in1"])
        .arg(&r1)
        .args(["--in2"])
        .arg(&r2)
        .args([
            "--overlap_len_require",
            "30",
            "--overlap_diff_limit",
            "5",
            "--overlap_diff_percent_limit",
            "20",
            "--disable_adapter_trimming",
            "--disable_quality_filtering",
            "--disable_length_filtering",
            "--disable_trim_poly_g",
            "--json",
            "/dev/null",
            "--html",
            "/dev/null",
        ])
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("spawn fastp");
    assert!(fp.success(), "fastp --merge failed");

    let ours = Command::new(bin())
        .args(["--in1"])
        .arg(&r1)
        .args(["--in2"])
        .arg(&r2)
        .args(["-m"])
        .arg(&ours_out)
        .output()
        .expect("spawn ours");
    assert!(
        ours.status.success(),
        "rsomics-fastq-merge failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let fp_seqs = merged_seqs(&fs::read_to_string(&fp_out).unwrap());
    let our_seqs = merged_seqs(&fs::read_to_string(&ours_out).unwrap());
    assert!(!fp_seqs.is_empty(), "fastp produced no merged reads");
    assert_eq!(
        our_seqs, fp_seqs,
        "merged-sequence multiset differs from fastp v{PINNED}"
    );
}
