use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-flank"))
}
fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn flank_both_50() {
    let out = Command::new(bin())
        .args(["-i"])
        .arg(fixture("one.bed"))
        .args(["-g"])
        .arg(fixture("genome.txt"))
        .args(["-b", "50"])
        .output()
        .expect("spawn");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let s = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = s.trim().lines().collect();
    assert_eq!(lines.len(), 2, "left flank + right flank = 2 intervals");
    assert_eq!(lines[0], "chr1\t50\t100", "left flank: 100-50=50 to 100");
    assert_eq!(lines[1], "chr1\t200\t250", "right flank: 200 to 200+50=250");
}
