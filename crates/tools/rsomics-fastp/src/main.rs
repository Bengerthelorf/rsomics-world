use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

use rsomics_fastp::filter::FilterConfig;
use rsomics_fastp::polyg::PolyGConfig;
use rsomics_fastp::trim::AdapterConfig;

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

    /// JSON report path (fastp-compatible schema subset).
    #[arg(short = 'j', long = "json")]
    json: Option<PathBuf>,

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
}

fn main() -> Result<()> {
    let args = Args::parse();
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
    match (args.in2, args.out2) {
        (Some(in2), Some(out2)) => {
            rsomics_fastp::io::process_pe(
                &args.in1,
                &in2,
                &args.out1,
                &out2,
                args.json.as_deref(),
                cfg,
                adapter.as_ref(),
                polyg,
            )?;
        }
        (None, None) => {
            rsomics_fastp::io::process_se(
                &args.in1,
                &args.out1,
                args.json.as_deref(),
                cfg,
                adapter.as_ref(),
                polyg,
            )?;
        }
        _ => {
            anyhow::bail!("--in2 and --out2 must both be set, or both unset");
        }
    }
    Ok(())
}
