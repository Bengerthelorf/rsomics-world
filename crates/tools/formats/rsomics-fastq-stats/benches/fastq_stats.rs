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
            // Phred 2..=40 over the sanger range, deterministic per seed.
            let q = 35 + ((rng >> 40) % 5) as u8;
            w.write_all(&[q]).unwrap();
        }
        w.write_all(b"\n").unwrap();
    }
}

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rsomics-fastq-stats-bench-{N_READS}x{READ_LEN}.fq"));
    if !p.exists() {
        synth_fastq(&p);
    }
    p
}

fn seqkit_available() -> bool {
    Command::new("seqkit")
        .arg("version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn bench(c: &mut Criterion) {
    let fixture = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-fastq-stats");
    let mut group = c.benchmark_group(format!("fastq_stats/{N_READS}x{READ_LEN}"));
    group.sample_size(20);
    group.bench_function("rsomics-fastq-stats", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args(["--all", fixture.to_str().unwrap()])
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "rsomics-fastq-stats failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });
    if seqkit_available() {
        let path = fixture.to_str().unwrap().to_string();
        group.bench_function("seqkit-stats", |b| {
            b.iter(|| {
                let out = Command::new("seqkit")
                    .args(["stats", "-a", &path])
                    .output()
                    .expect("seqkit run");
                assert!(
                    out.status.success(),
                    "seqkit failed: {}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("seqkit not on PATH — skipping upstream comparison");
    }
    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
