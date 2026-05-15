use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_RECORDS: usize = 5_000;
const SEQ_LEN: usize = 1_000;
const SEED: u64 = 0x00C0_FFEE;

fn synth_fasta(path: &PathBuf) {
    let f = File::create(path).expect("create bench fixture");
    let mut w = BufWriter::new(f);
    let mut rng = SEED;
    for i in 0..N_RECORDS {
        writeln!(w, ">contig_{i}").unwrap();
        for _ in 0..SEQ_LEN {
            rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            w.write_all(&[b"ACGT"[((rng >> 33) & 3) as usize]]).unwrap();
        }
        w.write_all(b"\n").unwrap();
    }
}

fn ensure_fixture() -> PathBuf {
    let mut p = std::env::temp_dir();
    p.push(format!(
        "rsomics-fasta-stats-bench-{N_RECORDS}x{SEQ_LEN}.fa"
    ));
    if !p.exists() {
        synth_fasta(&p);
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
    let ours = env!("CARGO_BIN_EXE_rsomics-fasta-stats");
    let mut group = c.benchmark_group(format!("fasta_stats/{N_RECORDS}x{SEQ_LEN}"));
    group.sample_size(20);
    group.bench_function("rsomics-fasta-stats", |b| {
        b.iter(|| {
            let out = Command::new(ours).arg(&fixture).output().expect("ours run");
            assert!(
                out.status.success(),
                "rsomics-fasta-stats failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });
    if seqkit_available() {
        let path = fixture.to_str().unwrap().to_string();
        group.bench_function("seqkit-stats", |b| {
            b.iter(|| {
                let out = Command::new("seqkit")
                    .args(["stats", &path])
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
