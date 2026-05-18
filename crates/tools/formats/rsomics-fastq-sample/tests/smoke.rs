use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastq-sample"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn sample_50_percent_gives_roughly_half() {
    let out = Command::new(bin())
        .args(["-p", "0.5", "--seed", "42"])
        .arg(fixture("hundred.fq"))
        .output()
        .expect("spawn");
    assert!(out.status.success());
    let count = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| l.starts_with('@'))
        .count();
    assert!(
        count > 20 && count < 80,
        "50% of 100 should be ~50, got {count}"
    );
}

#[test]
fn sample_is_deterministic_with_seed() {
    let run = |seed: &str| {
        let out = Command::new(bin())
            .args(["-p", "0.3", "--seed", seed])
            .arg(fixture("hundred.fq"))
            .output()
            .expect("spawn");
        String::from_utf8(out.stdout).unwrap()
    };
    assert_eq!(run("42"), run("42"), "same seed must give same output");
    assert_ne!(run("42"), run("99"), "different seeds should differ");
}
