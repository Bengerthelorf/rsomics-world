//! Single-pipeline rsomics-* template — every flag sits on one `Args` struct.
//! Appropriate when the tool does exactly one thing end-to-end.
//!
//! For multi-subcommand tools (`view` / `sort` / `index` / …) see
//! [`rsomics-bam`](https://crates.io/crates/rsomics-bam) instead — that
//! crate's `src/main.rs` and `src/cmd/` layout is the right template to
//! clone when the tool has more than one distinct operating mode.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, StderrLog, run};

use rsomics_fastp::filter::FilterConfig;
use rsomics_fastp::polyg::PolyGConfig;
use rsomics_fastp::trim::AdapterConfig;
use rsomics_fastp::umi::{UmiConfig, UmiLoc};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-fastp",
    version,
    about = "Pure-Rust FASTQ preprocessor — output-compatible with OpenGene/fastp.",
    long_about = None,
)]
struct Args {
    /// Input FASTQ R1. Compression auto-detected by file content.
    #[arg(short = 'i', long = "in1")]
    in1: PathBuf,

    /// Output FASTQ R1. Gzip-encoded iff the path ends in `.gz`.
    #[arg(short = 'o', long = "out1")]
    out1: PathBuf,

    /// Input FASTQ R2 (paired-end). When set, also requires `--out2`.
    #[arg(short = 'I', long = "in2")]
    in2: Option<PathBuf>,

    /// Output FASTQ R2 (paired-end). Required iff `--in2` is set.
    #[arg(short = 'O', long = "out2")]
    out2: Option<PathBuf>,

    /// JSON report output path (fastp-compatible schema subset). Long
    /// form is `--json_report` to leave the workspace-wide `--json` flag
    /// in [`CommonFlags`] free as the standard machine-readable-output
    /// switch.
    #[arg(short = 'j', long = "json_report")]
    json_report: Option<PathBuf>,

    /// Per-base quality threshold (Phred); bases below this count as unqualified.
    #[arg(long = "qualified_quality_phred", default_value_t = 15)]
    qualified_quality_phred: u8,

    /// Reject a read if the fraction of unqualified bases (percent) exceeds this.
    #[arg(long = "unqualified_percent_limit", default_value_t = 40)]
    unqualified_percent_limit: u32,

    /// Reject reads shorter than this many bases.
    #[arg(long = "length_required", default_value_t = 15)]
    length_required: usize,

    /// Reject reads with more than this many N bases.
    #[arg(long = "n_base_limit", default_value_t = 5)]
    n_base_limit: usize,

    /// Adapter sequence to trim from the 3' end. Defaults to Illumina `TruSeq`
    /// (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`). Pass an empty string to disable.
    #[arg(
        long = "adapter_sequence",
        default_value = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    )]
    adapter_sequence: String,

    /// Trim 3' poly-G runs (2-color chemistry artifact on NextSeq/NovaSeq).
    /// Off by default; pass `--trim_poly_g` to enable.
    #[arg(long = "trim_poly_g", default_value_t = false)]
    trim_poly_g: bool,

    /// Minimum poly-G run length to trim.
    #[arg(long = "poly_g_min_len", default_value_t = 10)]
    poly_g_min_len: usize,

    /// Enable UMI extraction. Off by default.
    #[arg(long = "umi", default_value_t = false)]
    umi: bool,

    /// Which mate holds the UMI: `read1` or `read2`. SE only supports `read1`.
    #[arg(long = "umi_loc", default_value = "read1")]
    umi_loc: String,

    /// UMI length (bases to take from the 5' end of the donor mate).
    #[arg(long = "umi_len", default_value_t = 8)]
    umi_len: usize,

    #[command(flatten)]
    common: CommonFlags,
}

fn pipeline(args: Args) -> Result<()> {
    let log = StderrLog::from_flags(&args.common);
    if args.common.threads.is_some_and(|n| n > 1) {
        log.info(format_args!(
            "note: rsomics-fastp's SE/PE pipeline is single-threaded; \
             --threads={} sizes the global rayon pool but does not \
             accelerate this run",
            args.common.threads.unwrap()
        ));
    }
    let cfg = FilterConfig {
        qualified_quality_phred: args.qualified_quality_phred,
        unqualified_percent_limit: args.unqualified_percent_limit,
        length_required: args.length_required,
        n_base_limit: args.n_base_limit,
    };
    let adapter = if args.adapter_sequence.is_empty() {
        None
    } else {
        Some(AdapterConfig {
            sequence: args.adapter_sequence.as_bytes().to_vec(),
            min_match_len: 5,
            max_mismatch_rate: 0.2,
        })
    };
    let polyg = args.trim_poly_g.then_some(PolyGConfig {
        min_len: args.poly_g_min_len,
        max_mismatches: 0,
    });
    let umi = if args.umi {
        let loc = match args.umi_loc.as_str() {
            "read1" => UmiLoc::Read1,
            "read2" => UmiLoc::Read2,
            other => {
                return Err(RsomicsError::ConfigError(format!(
                    "--umi_loc must be read1 or read2, got: {other}"
                )));
            }
        };
        Some(UmiConfig {
            loc,
            len: args.umi_len,
        })
    } else {
        None
    };
    match (args.in2, args.out2) {
        (Some(in2), Some(out2)) => {
            let outcome = rsomics_fastp::io::process_pe(
                &args.in1,
                &in2,
                &args.out1,
                &out2,
                args.json_report.as_deref(),
                cfg,
                adapter.as_ref(),
                polyg,
                umi,
            )?;
            let total = outcome.filtering.passed_filter_reads
                + outcome.filtering.low_quality_reads
                + outcome.filtering.too_many_n_reads
                + outcome.filtering.too_short_reads;
            log.info(format_args!(
                "PE: kept {passed}/{total} pairs",
                passed = outcome.filtering.passed_filter_reads,
            ));
        }
        (None, None) => {
            let outcome = rsomics_fastp::io::process_se(
                &args.in1,
                &args.out1,
                args.json_report.as_deref(),
                cfg,
                adapter.as_ref(),
                polyg,
                umi,
            )?;
            let total = outcome.filtering.passed_filter_reads
                + outcome.filtering.low_quality_reads
                + outcome.filtering.too_many_n_reads
                + outcome.filtering.too_short_reads;
            log.info(format_args!(
                "SE: kept {passed}/{total} reads",
                passed = outcome.filtering.passed_filter_reads,
            ));
        }
        _ => {
            return Err(RsomicsError::ConfigError(
                "--in2 and --out2 must both be set, or both unset".into(),
            ));
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    let args = Args::parse();
    let common = args.common.clone();
    run(&common, || pipeline(args))
}
