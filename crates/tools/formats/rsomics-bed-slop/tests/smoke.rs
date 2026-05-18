use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf { PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bed-slop")) }
fn fixture(name: &str) -> PathBuf { PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden").join(name) }

#[test]
fn slop_both_sides() {
    let out = Command::new(bin())
        .args(["-i"]).arg(fixture("two.bed"))
        .args(["-g"]).arg(fixture("genome.txt"))
        .args(["-b", "5"])
        .output().expect("spawn");
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let s = String::from_utf8(out.stdout).unwrap();
    let lines: Vec<&str> = s.trim().lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0], "chr1\t5\t25", "10-5=5, 20+5=25");
    assert_eq!(lines[1], "chr1\t45\t65", "50-5=45, 60+5=65");
}

#[test]
fn slop_clamps_to_zero() {
    let out = Command::new(bin())
        .args(["-i"]).arg(fixture("two.bed"))
        .args(["-g"]).arg(fixture("genome.txt"))
        .args(["-b", "20"])
        .output().expect("spawn");
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let first_line = s.lines().next().unwrap();
    assert!(first_line.starts_with("chr1\t0\t"), "start should clamp to 0: {first_line}");
}
