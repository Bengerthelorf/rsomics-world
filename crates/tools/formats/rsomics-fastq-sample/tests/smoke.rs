use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf { PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-sample")) }
fn fixture(name: &str) -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(name) }

#[test]
fn sample_50_percent() {
    let out = Command::new(bin())
        .args(["-p", "0.5"])
        .arg(fixture("hundred.fq"))
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let count = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('@'))
        .count();
    assert!(count > 10 && count < 90, "~50% of 100 should be 10-90, got {count}");
}
