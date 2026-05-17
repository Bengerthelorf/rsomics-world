use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_PAIRS: usize = 100_000;
const FRAG: usize = 200;
const READ: usize = 150;
const SEED: u64 = 0x00C0_FFEE;

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

fn ensure_fixture() -> (PathBuf, PathBuf) {
    let mut d = std::env::temp_dir();
    d.push(format!("rsomics-fastq-merge-bench-{N_PAIRS}"));
    let r1 = d.with_extension("r1.fq");
    let r2 = d.with_extension("r2.fq");
    if r1.exists() && r2.exists() {
        return (r1, r2);
    }
    let mut w1 = BufWriter::new(File::create(&r1).unwrap());
    let mut w2 = BufWriter::new(File::create(&r2).unwrap());
    let mut rng = SEED;
    for i in 0..N_PAIRS {
        let mut frag = Vec::with_capacity(FRAG);
        for _ in 0..FRAG {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            frag.push(b"ACGT"[((rng >> 33) & 3) as usize]);
        }
        let q = [b'I'; READ];
        writeln!(w1, "@f{i}/1").unwrap();
        w1.write_all(&frag[..READ]).unwrap();
        w1.write_all(b"\n+\n").unwrap();
        w1.write_all(&q).unwrap();
        w1.write_all(b"\n").unwrap();
        writeln!(w2, "@f{i}/2").unwrap();
        w2.write_all(&revcomp(&frag[FRAG - READ..])).unwrap();
        w2.write_all(b"\n+\n").unwrap();
        w2.write_all(&q).unwrap();
        w2.write_all(b"\n").unwrap();
    }
    (r1, r2)
}

fn bench(c: &mut Criterion) {
    let (r1, r2) = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-fastq-merge");
    let mut group = c.benchmark_group(format!("fastq_merge/{N_PAIRS}pairs"));
    group.sample_size(20);
    group.bench_function("rsomics-fastq-merge", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["--in1"])
                .arg(&r1)
                .args(["--in2"])
                .arg(&r2)
                .args(["-m", "/dev/null"])
                .output()
                .expect("run");
            assert!(out.status.success());
        });
    });
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
