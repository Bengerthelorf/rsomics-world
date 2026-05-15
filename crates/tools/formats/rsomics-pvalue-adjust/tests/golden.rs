use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-pvalue-adjust"))
}

#[test]
fn bh_single_column_matches_textbook() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.txt");
    let out = tmp.path().join("out.txt");
    let mut f = std::fs::File::create(&inp).unwrap();
    for p in [0.01, 0.02, 0.04, 0.10, 0.50] {
        writeln!(f, "{p}").unwrap();
    }
    let status = Command::new(ours())
        .args(["-i", inp.to_str().unwrap(), "-o", out.to_str().unwrap()])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = got.lines().collect();
    assert!(lines[0].ends_with("0.050000"), "{}", lines[0]);
    assert!(lines[1].ends_with("0.050000"), "{}", lines[1]);
    assert!(lines[2].contains("0.066"), "{}", lines[2]);
    assert!(lines[3].ends_with("0.125000"), "{}", lines[3]);
    assert!(lines[4].ends_with("0.500000"), "{}", lines[4]);
}

#[test]
fn bonferroni_caps_at_one() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.txt");
    let out = tmp.path().join("out.txt");
    let mut f = std::fs::File::create(&inp).unwrap();
    for p in [0.5, 0.4, 0.6] {
        writeln!(f, "{p}").unwrap();
    }
    let status = Command::new(ours())
        .args([
            "-i",
            inp.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-m",
            "bonferroni",
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(&out).unwrap();
    assert!(got.lines().next().unwrap().ends_with("1.000000"));
}

#[test]
fn picks_column_from_multi_column_input() {
    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.tsv");
    let out = tmp.path().join("out.tsv");
    let mut f = std::fs::File::create(&inp).unwrap();
    writeln!(f, "gene_a\t1.5\t0.01").unwrap();
    writeln!(f, "gene_b\t2.0\t0.02").unwrap();
    let status = Command::new(ours())
        .args([
            "-i",
            inp.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "-c",
            "3",
        ])
        .status()
        .unwrap();
    assert!(status.success());
    let got = std::fs::read_to_string(out).unwrap();
    assert!(got.starts_with("gene_a\t1.5\t0.01\t"), "{got}");
}
