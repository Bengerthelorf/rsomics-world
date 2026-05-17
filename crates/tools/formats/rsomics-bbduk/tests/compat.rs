use std::collections::BTreeSet;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

const ADAPTER: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
const N: usize = 400;
const LEN: usize = 60;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-bbduk"))
}

fn bbduk_present() -> Option<String> {
    let out = Command::new("bbduk.sh").arg("--version").output().ok()?;
    let v = format!(
        "{}{}",
        String::from_utf8_lossy(&out.stderr),
        String::from_utf8_lossy(&out.stdout)
    );
    v.contains("BBMap").then_some(v)
}

// all reads 60 bp ≥ minlength so the length filter is inert — the test isolates the k-mer decision
fn synth(path: &Path) {
    let mut w = BufWriter::new(fs::File::create(path).unwrap());
    let mut rng = 0x1234_5678_9abc_def0u64;
    let rand_base = |rng: &mut u64| {
        *rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        b"ACGT"[((*rng >> 33) & 3) as usize]
    };
    for i in 0..N {
        let mut seq = Vec::with_capacity(LEN);
        if i % 2 == 0 {
            while seq.len() < 13 {
                seq.push(rand_base(&mut rng));
            }
            seq.extend_from_slice(ADAPTER.as_bytes());
            while seq.len() < LEN {
                seq.push(rand_base(&mut rng));
            }
        } else {
            while seq.len() < LEN {
                seq.push(rand_base(&mut rng));
            }
        }
        seq.truncate(LEN);
        writeln!(w, "@read{i}").unwrap();
        w.write_all(&seq).unwrap();
        w.write_all(b"\n+\n").unwrap();
        w.write_all(&[b'I'; LEN]).unwrap();
        w.write_all(b"\n").unwrap();
    }
}

fn surviving_ids(fq: &str) -> BTreeSet<String> {
    let mut ids = BTreeSet::new();
    let mut lines = fq.lines();
    while let Some(h) = lines.next() {
        if let Some(name) = h.strip_prefix('@') {
            ids.insert(name.split_whitespace().next().unwrap_or("").to_string());
        }
        lines.next();
        lines.next();
        lines.next();
    }
    ids
}

#[test]
fn kfilter_survivor_set_matches_bbduk() {
    // gate on bbduk.sh presence, not a pinned version: kfilter k-mer-set membership is version-stable (unlike heuristic adapter detectors)
    let Some(ver) = bbduk_present() else {
        eprintln!(
            "SKIP: bbduk.sh not on PATH — kfilter oracle unavailable \
             (authoritative on the CI/publish lane where BBMap is installed)"
        );
        return;
    };
    eprintln!("bbduk oracle: {}", ver.lines().next().unwrap_or("?"));

    let tmp = tempfile::tempdir().unwrap();
    let inp = tmp.path().join("in.fq");
    let bb_out = tmp.path().join("bb.fq");
    let our_out = tmp.path().join("our.fq");
    synth(&inp);

    let bb = Command::new("bbduk.sh")
        .arg("-Xmx1g")
        .arg(format!("in={}", inp.display()))
        .arg(format!("out={}", bb_out.display()))
        .arg(format!("literal={ADAPTER}"))
        .args(["k=27", "threads=1", "overwrite=t"])
        .output()
        .expect("spawn bbduk.sh");
    assert!(
        bb.status.success(),
        "bbduk.sh failed: {}",
        String::from_utf8_lossy(&bb.stderr)
    );

    let ours = Command::new(bin())
        .args(["--literal", ADAPTER, "-i"])
        .arg(&inp)
        .args(["-o"])
        .arg(&our_out)
        .output()
        .expect("spawn ours");
    assert!(
        ours.status.success(),
        "rsomics-bbduk failed: {}",
        String::from_utf8_lossy(&ours.stderr)
    );

    let bb_ids = surviving_ids(&fs::read_to_string(&bb_out).unwrap());
    let our_ids = surviving_ids(&fs::read_to_string(&our_out).unwrap());
    assert!(
        !bb_ids.is_empty(),
        "bbduk kept no reads — fixture/oracle wrong"
    );
    assert!(
        bb_ids.len() < N,
        "bbduk removed nothing — contaminants not detected"
    );
    assert_eq!(
        our_ids,
        bb_ids,
        "kfilter survivor set differs from bbduk.sh ({} vs {} kept)",
        our_ids.len(),
        bb_ids.len()
    );
}
