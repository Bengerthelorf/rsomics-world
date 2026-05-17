use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-fastqc"))
}

fn fixture() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/tiny.fq")
}

fn run_stdout(args: &[&str]) -> String {
    let out = Command::new(bin()).args(args).output().expect("spawn");
    assert!(
        out.status.success(),
        "failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

#[test]
fn emits_parseable_fastqc_data() {
    let f = fixture();
    let data = run_stdout(&["--stdout", f.to_str().unwrap()]);

    assert!(data.starts_with("##FastQC\t"), "version header");
    assert!(data.contains(">>Basic Statistics\tpass\n"));
    assert!(data.contains("Total Sequences\t3\n"));
    // tiny.fq has lengths 10/10/8 → unequal ⇒ Sequence Length WARN.
    assert!(data.contains(">>Sequence Length Distribution\twarn\n"));
    // Every opened module is closed.
    let opens = data.matches(">>").count();
    let ends = data.matches(">>END_MODULE").count();
    assert_eq!(opens, ends * 2, "each >>Name pairs with one >>END_MODULE");
}

#[test]
fn rejects_non_fastq() {
    let out = Command::new(bin())
        .args(["--stdout", "/dev/null"])
        .output()
        .expect("spawn");
    assert!(!out.status.success(), "empty/non-FASTQ must fail loud");
}
