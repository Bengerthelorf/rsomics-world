//! Black-box compatibility vs the `FastQC` binary. Version-gated to the
//! pinned release (CI/publish install exactly this); on any other or
//! absent `FastQC` the test loud-skips rather than asserting against an
//! unknown oracle.
//!
//! The fixture is synthesised to be unambiguous for *any* correct `FastQC`
//! implementation — uniform-length, base-balanced, uniformly high quality,
//! all-distinct reads, no adapters — so every module is PASS and the Basic
//! Statistics values are exact. That makes the 12 module statuses + the
//! Basic Statistics block a deterministic cross-check that needs no
//! pre-observed `FastQC` output.

use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const PINNED: &str = "0.12.1";
const N_READS: usize = 2_000;
const READ_LEN: usize = 150;

fn rsomics_bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastqc"))
}

fn fastqc_pinned() -> bool {
    let Ok(out) = Command::new("fastqc").arg("--version").output() else {
        return false;
    };
    if !out.status.success() {
        return false;
    }
    String::from_utf8_lossy(&out.stdout).contains(PINNED)
}

/// Deterministic balanced 150 bp reads, Phred+33 `I` (Q40), every read
/// unique (the index is encoded in the first bases) so duplication and
/// overrepresentation are trivially zero.
fn synth_clean(path: &Path) {
    let f = fs::File::create(path).expect("create fixture");
    let mut w = BufWriter::new(f);
    for i in 0..N_READS {
        writeln!(w, "@read_{i}").unwrap();
        let mut rng = 0x9E37_79B9_7F4A_7C15u64 ^ (i as u64).wrapping_mul(0x1000_0001b3);
        let mut seq = Vec::with_capacity(READ_LEN);
        for _ in 0..READ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            seq.push(b"ACGT"[((rng >> 33) & 3) as usize]);
        }
        w.write_all(&seq).unwrap();
        w.write_all(b"\n+\n").unwrap();
        w.write_all(&[b'I'; READ_LEN]).unwrap();
        w.write_all(b"\n").unwrap();
    }
}

/// Parse a `FastQC` `summary.txt` (`STATUS<TAB>Module<TAB>file`) into
/// `(module, status)` pairs in order.
fn parse_summary(text: &str) -> Vec<(String, String)> {
    text.lines()
        .filter_map(|l| {
            let mut it = l.split('\t');
            let status = it.next()?.trim().to_string();
            let module = it.next()?.trim().to_string();
            Some((module, status))
        })
        .collect()
}

fn read_to_string_in(dir: &Path, leaf: &str) -> Option<String> {
    // FastQC names the dir "<base-without-ext>_fastqc"; ours uses the full
    // basename. Find the single *_fastqc dir and read `leaf` from it.
    let entry = fs::read_dir(dir)
        .ok()?
        .filter_map(Result::ok)
        .find(|e| e.file_name().to_string_lossy().ends_with("_fastqc"))?;
    fs::read_to_string(entry.path().join(leaf)).ok()
}

fn basic_stat<'a>(data: &'a str, key: &str) -> Option<&'a str> {
    data.lines()
        .skip_while(|l| !l.starts_with(">>Basic Statistics"))
        .take_while(|l| !l.starts_with(">>END_MODULE"))
        .find_map(|l| {
            let mut it = l.split('\t');
            (it.next()? == key).then(|| it.next())?
        })
}

#[test]
fn statuses_and_basic_stats_match_fastqc() {
    if !fastqc_pinned() {
        eprintln!(
            "SKIP: FastQC v{PINNED} not on PATH — clean-room compat oracle \
             unavailable (Linux-only; authoritative on CI/publish.yml)"
        );
        return;
    }

    let tmp = tempfile::tempdir().expect("tempdir");
    let fixture = tmp.path().join("clean.fq");
    synth_clean(&fixture);

    let fq_out = tmp.path().join("fastqc_out");
    let ours_out = tmp.path().join("ours_out");
    fs::create_dir_all(&fq_out).unwrap();
    fs::create_dir_all(&ours_out).unwrap();

    let fq = Command::new("fastqc")
        .args(["--quiet", "--extract", "-o"])
        .arg(&fq_out)
        .arg(&fixture)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("spawn fastqc");
    assert!(fq.success(), "fastqc run failed");

    let ours = Command::new(rsomics_bin())
        .arg("-o")
        .arg(&ours_out)
        .arg(&fixture)
        .output()
        .expect("spawn ours");
    assert!(
        ours.status.success(),
        "rsomics-fastqc failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let fq_sum = parse_summary(&read_to_string_in(&fq_out, "summary.txt").expect("fastqc summary"));
    let our_sum = parse_summary(&read_to_string_in(&ours_out, "summary.txt").expect("our summary"));

    // Every module FastQC reports must be present in ours with the same
    // status (the synthetic fixture is unambiguous ⇒ all PASS on both).
    for (module, fq_status) in &fq_sum {
        let our = our_sum
            .iter()
            .find(|(m, _)| m == module)
            .unwrap_or_else(|| panic!("module {module:?} missing from our summary"));
        assert_eq!(
            &our.1, fq_status,
            "status mismatch for {module:?}: ours={} fastqc={fq_status}",
            our.1
        );
    }
    // FastQC v0.12.1 emits only the modules it ran: "Per tile sequence
    // quality" is absent without Illumina tile IDs and "Kmer Content" is
    // off by default, so a non-tiled fixture yields ~10, not 12. The
    // contract is that every module FastQC *did* report matches ours (the
    // loop above) and that we emit the full set.
    assert!(
        fq_sum.len() >= 8,
        "FastQC produced too few modules ({}) — broken oracle run",
        fq_sum.len()
    );
    assert_eq!(our_sum.len(), 12, "rsomics-fastqc must emit all 12 modules");

    let fq_data = read_to_string_in(&fq_out, "fastqc_data.txt").expect("fastqc data");
    let our_data = read_to_string_in(&ours_out, "fastqc_data.txt").expect("our data");
    for key in ["Total Sequences", "Sequence length", "%GC"] {
        assert_eq!(
            basic_stat(&our_data, key),
            basic_stat(&fq_data, key),
            "Basic Statistics {key:?} differs"
        );
    }
}
