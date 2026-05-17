use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_READS: usize = 50_000;
const READ_LEN: usize = 150;
const SEED: u64 = 0x00C0_FFEE;

fn synth_fastq(path: &PathBuf) {
    let f = File::create(path).expect("create bench fixture");
    let mut w = BufWriter::new(f);
    let mut rng = SEED;
    for i in 0..N_READS {
        writeln!(w, "@read_{i}").unwrap();
        for _ in 0..READ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            w.write_all(&[b"ACGT"[((rng >> 33) & 3) as usize]]).unwrap();
        }
        w.write_all(b"\n+\n").unwrap();
        for _ in 0..READ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            let q = 35 + ((rng >> 40) % 5) as u8;
            w.write_all(&[q]).unwrap();
        }
        w.write_all(b"\n").unwrap();
    }
}

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rsomics-fastqc-bench-{N_READS}x{READ_LEN}.fq"));
    if !p.exists() {
        synth_fastq(&p);
    }
    p
}

fn bench(c: &mut Criterion) {
    let fixture = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-fastqc");
    let mut group = c.benchmark_group(format!("fastqc/{N_READS}x{READ_LEN}"));
    group.sample_size(20);
    group.bench_function("rsomics-fastqc", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["--stdout", fixture.to_str().unwrap()])
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "rsomics-fastqc failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
