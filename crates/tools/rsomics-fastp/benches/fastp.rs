//! Criterion benches for `rsomics-fastp`.
//!
//! Generates a 100k-read synthetic single-end FASTQ at bench startup, then
//! measures `process_se` throughput on it. If upstream `fastp` is on PATH a
//! second sample compares wall-clock per-read against it; otherwise that
//! arm is skipped.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rsomics_fastp::filter::FilterConfig;
use rsomics_fastp::io::process_se;

const N_READS: usize = 100_000;
const READ_LEN: usize = 100;

/// Build a deterministic synthetic FASTQ in a temp dir. Sequence content is
/// a fixed 4-base cycle so adapter / poly-G logic doesn't accidentally
/// trigger; qualities are all `I` (Phred 40) so the read passes the default
/// filter. 100k reads × 100 bp ≈ 30 MB on disk.
fn synth_fastq(path: &Path) -> std::io::Result<u64> {
    let mut f = std::fs::File::create(path)?;
    let seq: Vec<u8> = (0..READ_LEN).map(|i| b"ACGT"[i % 4]).collect();
    let qual: Vec<u8> = vec![b'I'; READ_LEN];
    for i in 0..N_READS {
        writeln!(f, "@read_{i:08}")?;
        f.write_all(&seq)?;
        f.write_all(b"\n+\n")?;
        f.write_all(&qual)?;
        f.write_all(b"\n")?;
    }
    f.flush()?;
    Ok(std::fs::metadata(path)?.len())
}

fn fastp_on_path() -> bool {
    Command::new("fastp")
        .arg("--version")
        .output()
        .is_ok_and(|out| out.status.success())
}

fn bench_se_filter(c: &mut Criterion) {
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let input = tmpdir.path().join("synth.fastq");
    let bytes = synth_fastq(&input).expect("write fixture");

    let mut group = c.benchmark_group("se_filter");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(8));
    group.throughput(Throughput::Bytes(bytes));

    group.bench_function(BenchmarkId::new("rsomics_fastp", "synth_100k"), |b| {
        b.iter_with_setup(
            || {
                tempfile::Builder::new()
                    .suffix(".fastq")
                    .tempfile_in(tmpdir.path())
                    .expect("tempfile")
            },
            |out| {
                process_se(
                    &input,
                    out.path(),
                    None,
                    FilterConfig::default(),
                    None,
                    None,
                    None,
                )
                .expect("process_se");
            },
        );
    });

    if fastp_on_path() {
        group.bench_function(BenchmarkId::new("upstream_fastp", "synth_100k"), |b| {
            b.iter_with_setup(
                || {
                    tempfile::Builder::new()
                        .suffix(".fastq")
                        .tempfile_in(tmpdir.path())
                        .expect("tempfile")
                },
                |out| {
                    let out_path: PathBuf = out.path().to_path_buf();
                    let status = Command::new("fastp")
                        .arg("-i")
                        .arg(&input)
                        .arg("-o")
                        .arg(&out_path)
                        .arg("--disable_adapter_trimming")
                        .arg("--json")
                        .arg(out_path.with_extension("fastp.json"))
                        .arg("--html")
                        .arg(out_path.with_extension("fastp.html"))
                        .arg("--thread")
                        .arg("1")
                        .stderr(std::process::Stdio::null())
                        .stdout(std::process::Stdio::null())
                        .status()
                        .expect("spawn fastp");
                    assert!(status.success());
                },
            );
        });
    } else {
        eprintln!("note: upstream fastp not on PATH; comparison bench skipped");
    }

    group.finish();
}

criterion_group!(benches, bench_se_filter);
criterion_main!(benches);
