use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N: usize = 200_000;
const LEN: usize = 150;
const SEED: u64 = 0x00C0_FFEE;
const ADAPTER: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rsomics-bbduk-bench-{N}.fq"));
    if p.exists() {
        return p;
    }
    let mut w = BufWriter::new(File::create(&p).unwrap());
    let mut rng = SEED;
    let q = [b'I'; LEN];
    for i in 0..N {
        let mut seq = Vec::with_capacity(LEN);
        for _ in 0..LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            seq.push(b"ACGT"[((rng >> 33) & 3) as usize]);
        }
        writeln!(w, "@read{i}").unwrap();
        w.write_all(&seq).unwrap();
        w.write_all(b"\n+\n").unwrap();
        w.write_all(&q).unwrap();
        w.write_all(b"\n").unwrap();
    }
    p
}

fn bench(c: &mut Criterion) {
    let fq = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-bbduk");
    let mut group = c.benchmark_group(format!("bbduk/{N}reads"));
    group.sample_size(20);
    group.bench_function("rsomics-bbduk-kfilter", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["--literal", ADAPTER, "-i"])
                .arg(&fq)
                .args(["-o", "/dev/null"])
                .output()
                .expect("run");
            assert!(out.status.success());
        });
    });
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
