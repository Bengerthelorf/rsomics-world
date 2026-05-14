//! Criterion benches for `rsomics-bam`.
//!
//! Generates a synthetic 100k-record unaligned BAM at bench startup,
//! then measures `view -c` throughput on it. If upstream `samtools` is
//! on PATH a second sample compares against `samtools view -c`;
//! otherwise that arm is skipped.
//!
//! Note on metric choice: `view -c` is record-iteration bound, so
//! `Throughput::Elements(n)` (records/sec) is the meaningful comparison
//! axis — `Bytes` would conflate BGZF compression ratio with iteration
//! cost. Future rsomics-* benches: pick `Elements` for record-oriented
//! pipelines, `Bytes` only for true streaming throughput.

use std::path::Path;
use std::process::{Command, Stdio};
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rsomics_bam::cmd::view::count_records;
use rsomics_common::test_support::tool_on_path;
use rust_htslib::bam::header::{Header, HeaderRecord};
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{Format, Writer};

const N_RECORDS: u64 = 100_000;

fn write_synth_bam(path: &Path, n: u64) {
    let mut header = Header::new();
    let mut hd = HeaderRecord::new(b"HD");
    hd.push_tag(b"VN", "1.6");
    hd.push_tag(b"SO", "unsorted");
    header.push_record(&hd);

    let mut w = Writer::from_path(path, &header, Format::Bam).expect("create BAM");
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual = [30u8; 100];
    for i in 0..n {
        let mut rec = Record::new();
        let qname = format!("r_{i:08}");
        rec.set(qname.as_bytes(), None, seq, &qual);
        rec.set_unmapped();
        w.write(&rec).expect("write record");
    }
}

fn bench_view_count(c: &mut Criterion) {
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let input = tmpdir.path().join("synth.bam");
    write_synth_bam(&input, N_RECORDS);

    let mut group = c.benchmark_group("view_count");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(8));
    group.throughput(Throughput::Elements(N_RECORDS));

    group.bench_function(BenchmarkId::new("rsomics_bam", "synth_100k"), |b| {
        b.iter(|| {
            let n = count_records(&input).expect("count");
            assert_eq!(n, N_RECORDS);
        });
    });

    if tool_on_path("samtools") {
        group.bench_function(BenchmarkId::new("upstream_samtools", "synth_100k"), |b| {
            b.iter(|| {
                let status = Command::new("samtools")
                    .arg("view")
                    .arg("-c")
                    .arg(&input)
                    .stdout(Stdio::null())
                    .stderr(Stdio::null())
                    .status()
                    .expect("spawn samtools");
                assert!(status.success());
            });
        });
    } else {
        eprintln!("note: samtools not on PATH; comparison bench skipped");
    }

    group.finish();
}

criterion_group!(benches, bench_view_count);
criterion_main!(benches);
