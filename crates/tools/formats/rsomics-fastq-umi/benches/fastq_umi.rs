//! Criterion bench vs `fastp --umi` on a deterministic synthetic FASTQ.
//!
//! Fixture: 100 000 SE reads × 150 bp. Both binaries pinned to 1 thread,
//! extracting an 8 bp inline UMI from read1.

use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_READS: usize = 100_000;
const READ_LEN: usize = 150;
const SEED: u64 = 0x0000_BEEF;
const UMI_LEN: usize = 8;

fn synth_fastq(path: &PathBuf) {
    let f = File::create(path).expect("create bench fixture");
    let mut w = BufWriter::new(f);
    let mut rng = SEED;
    for i in 0..N_READS {
        writeln!(w, "@read_{i} 1:N:0:ACGT").unwrap();
        for _ in 0..READ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            w.write_all(&[b"ACGT"[((rng >> 33) % 4) as usize]]).unwrap();
        }
        writeln!(w).unwrap();
        writeln!(w, "+").unwrap();
        for _ in 0..READ_LEN {
            w.write_all(b"I").unwrap();
        }
        writeln!(w).unwrap();
    }
}

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!("rsomics-fastq-umi-bench-{N_READS}x{READ_LEN}.fq"));
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
    let ours = env!("CARGO_BIN_EXE_rsomics-fastq-umi");
    let outdir = tempfile::tempdir().expect("bench outdir");
    let out_ours = outdir.path().join("ours.fq");
    let out_fastp = outdir.path().join("fastp.fq");
    let json_fastp = outdir.path().join("fastp.json");
    let html_fastp = outdir.path().join("fastp.html");

    let mut group = c.benchmark_group(format!("fastq_umi/{N_READS}x{READ_LEN}"));
    group.sample_size(10);

    group.bench_function("rsomics-fastq-umi", |b| {
        b.iter(|| {
            let out = Command::new(ours)
                .args([
                    "-i",
                    fixture.to_str().unwrap(),
                    "-o",
                    out_ours.to_str().unwrap(),
                    "--umi_len",
                    &UMI_LEN.to_string(),
                    "-t",
                    "1",
                ])
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "rsomics-fastq-umi failed: {}",
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
                        "--umi",
                        "--umi_loc",
                        "read1",
                        "--umi_len",
                        &UMI_LEN.to_string(),
                        "--thread",
                        "1",
                        "-A",
                        "-G",
                        "-Q",
                        "-L",
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
