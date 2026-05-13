//! Golden-fixture tests that don't depend on the upstream fastp binary.
//! These exercise the internal pipeline end-to-end against small hand-curated
//! FASTQ files under `tests/golden/`.

use std::fs;
use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("golden")
        .join(name)
}

#[test]
fn copy_se_is_identity() {
    let input = fixture("se_basic.fastq");
    let tmp = tempfile::Builder::new()
        .suffix(".fastq")
        .tempfile()
        .expect("tempfile");

    rsomics_fastp::io::copy_se(&input, tmp.path()).expect("copy_se");

    let expected = fs::read(&input).expect("read input");
    let actual = fs::read(tmp.path()).expect("read output");
    assert_eq!(
        actual, expected,
        "identity copy must preserve the input byte-for-byte"
    );
}
