//! Criterion bench vs `fastp` on a deterministic synthetic FASTQ.
//!
//! Fixture: 100 000 SE reads × 150 bp with varied quality profiles to exercise
//! both the pass and fail paths of the filter. Both binaries pinned to 1 thread
//! with quality filter ON and length filter ON (fastp defaults).

use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_READS: usize = 100_000;
const READ_LEN: usize = 150;
const SEED: u64 = 0x0000_BEEF;

fn synth_fastq(path: &PathBuf) {
    let f = File::create(path).expect("create bench fixture");
    let mut w = BufWriter::new(f);
    let mut rng = SEED;
    for i in 0..N_READS {
        writeln!(w, "@read_{i}").unwrap();
        for _ in 0..READ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            w.write_all(&[b"ACGTN"[((rng >> 33) % 5) as usize]])
                .unwrap();
        }
        writeln!(w).unwrap();
        writeln!(w, "+").unwrap();
        // Quality profile: ~30% of reads have a run of 10 low-qual bases,
        // exercising the filter fail path at the fastp 40% default threshold.
        let low_run = if (i % 3) == 0 { 10 } else { 0 };
        for j in 0..READ_LEN {
            // Low-qual run at positions 20..30 for affected reads.
            let q = if j >= 20 && j < 20 + low_run {
                b'!'
            } else {
                b'I'
            };
            w.write_all(&[q]).unwrap();
        }
        writeln!(w).unwrap();
    }
}

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!(
        "rsomics-fastq-filter-bench-{N_READS}x{READ_LEN}.fq"
    ));
    if !p.exists() {
        synth_fastq(&p);
    }
    p
}

fn fastp_available() -> bool {
    Command::new("fastp")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn bench(c: &mut Criterion) {
    let fixture = ensure_fixture();
    let ours = env!("CARGO_BIN_EXE_rsomics-fastq-filter");
    let outdir = tempfile::tempdir().expect("bench outdir");
    let out_ours = outdir.path().join("ours.fq");
    let out_fastp = outdir.path().join("fastp.fq");
    let json_fastp = outdir.path().join("fastp.json");
    let html_fastp = outdir.path().join("fastp.html");

    let mut group = c.benchmark_group(format!("fastq_filter/{N_READS}x{READ_LEN}"));
    group.sample_size(10);

    group.bench_function("rsomics-fastq-filter", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args([
                    "-i",
                    fixture.to_str().unwrap(),
                    "-o",
                    out_ours.to_str().unwrap(),
                    "-t",
                    "1",
                ])
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "rsomics-fastq-filter failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    if fastp_available() {
        group.bench_function("fastp", |b| {
            b.iter(|| {
                let out = Command::new("fastp")
                    .args([
                        "-i",
                        fixture.to_str().unwrap(),
                        "-o",
                        out_fastp.to_str().unwrap(),
                        "--thread",
                        "1",
                        "-A",
                        "-G",
                        "--json",
                        json_fastp.to_str().unwrap(),
                        "--html",
                        html_fastp.to_str().unwrap(),
                    ])
                    .output()
                    .expect("fastp run");
                assert!(
                    out.status.success(),
                    "fastp failed: {}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("fastp not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
