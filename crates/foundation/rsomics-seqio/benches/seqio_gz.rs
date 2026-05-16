use std::io::Write;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_seqio::open_fastq;

fn make_gz_fixture(n_records: usize) -> tempfile::NamedTempFile {
    let mut fq = Vec::with_capacity(n_records * 60);
    for i in 0..n_records {
        write!(
            fq,
            "@read{i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
        )
        .unwrap();
    }
    let mut f = tempfile::Builder::new()
        .suffix(".fq.gz")
        .tempfile()
        .unwrap();
    {
        let mut enc =
            flate2::write::GzEncoder::new(f.as_file_mut(), flate2::Compression::default());
        enc.write_all(&fq).unwrap();
        enc.finish().unwrap();
    }
    f.flush().unwrap();
    f
}

fn bench_gz_decode(c: &mut Criterion) {
    let fixture = make_gz_fixture(100_000);
    let path = fixture.path().to_path_buf();

    c.bench_function("gz_decode_100k_records", |b| {
        b.iter(|| {
            let mut count: usize = 0;
            for r in open_fastq(&path).unwrap() {
                r.unwrap();
                count += 1;
            }
            assert_eq!(count, 100_000);
        });
    });
}

criterion_group!(benches, bench_gz_decode);
criterion_main!(benches);
