use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-random"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn deterministic_with_seed() {
    let run = || {
        let out = Command::new(ours())
            .args(["-g"])
            .arg(fixture("genome.txt"))
            .args(["-n", "10", "-l", "100", "--seed", "42"])
            .output()
            .expect("spawn");
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8(out.stdout).unwrap()
    };
    assert_eq!(run(), run(), "same seed = same output");
}
