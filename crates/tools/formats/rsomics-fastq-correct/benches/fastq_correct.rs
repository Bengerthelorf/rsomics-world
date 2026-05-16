use std::hint::black_box;
use std::io::Write;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_fastq_correct::{CorrectConfig, Pipeline};

fn write_synthetic_fastq(path: &std::path::Path, n: usize, len: usize) {
    let bases = b"ACGT";
    let mut f = std::fs::File::create(path).unwrap();
    for r in 0..n {
        let seq: Vec<u8> = (0..len).map(|i| bases[(r / 4 + i * 3) % 4]).collect();
        writeln!(f, "@r{r}").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f, "\n+").unwrap();
        writeln!(f, "{}", "I".repeat(len)).unwrap();
    }
}

fn bench_correct(c: &mut Criterion) {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("in.fq");
    let output = dir.path().join("out.fq");
    write_synthetic_fastq(&input, 2000, 100);
    let cfg = CorrectConfig {
        k: 17,
        ..CorrectConfig::default()
    };
    c.bench_function("correct_2k_reads_len100_k17", |b| {
        b.iter(|| {
            let p = Pipeline::new(&cfg, 4);
            black_box(p.run(black_box(&input), black_box(&output)).unwrap())
        });
    });
}

criterion_group!(benches, bench_correct);
criterion_main!(benches);
