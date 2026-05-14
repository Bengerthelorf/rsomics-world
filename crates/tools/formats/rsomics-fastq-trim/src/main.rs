use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta, run};

use rsomics_fastq_trim::{
    AdapterConfig, FixedTrimConfig, OverlapConfig, Pipeline, PipelineConfig, PolyXConfig,
    TrimReport,
};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// FASTQ adapter / `polyG` / `polyX` / fixed-length trimming.
///
/// Flag names follow fastp's `snake_case` convention so existing scripts
/// can swap `fastp` for `rsomics-fastq-trim` without re-learning the
/// option surface. Quality / N-content / sliding-window filtering, UMI
/// extraction, and per-cycle statistics live in sibling crates
/// (`rsomics-fastq-quality`, `rsomics-fastq-umi`, `rsomics-fastq-stats`).
#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-trim", version, about, long_about = None)]
#[allow(clippy::struct_excessive_bools)]
struct Cli {
    #[arg(short = 'i', long = "in1")]
    in1: PathBuf,

    #[arg(short = 'o', long = "out1")]
    out1: PathBuf,

    #[arg(short = 'I', long = "in2")]
    in2: Option<PathBuf>,

    #[arg(short = 'O', long = "out2")]
    out2: Option<PathBuf>,

    /// R1 adapter. Defaults to Illumina `TruSeq` R1 prefix when adapter
    /// trim is enabled and no sequence is supplied; pass an empty string
    /// to disable adapter trim entirely.
    #[arg(short = 'a', long = "adapter_sequence")]
    adapter_sequence: Option<String>,

    /// R2 adapter. Used only in PE mode; falls back to `TruSeq` R2 prefix
    /// if unspecified.
    #[arg(long = "adapter_sequence_r2")]
    adapter_sequence_r2: Option<String>,

    /// Disable static-sequence adapter trim. Useful when relying purely
    /// on PE overlap detection or running upstream-clean data.
    #[arg(short = 'A', long = "disable_adapter_trimming")]
    disable_adapter_trimming: bool,

    /// Force poly-G trim (default-off; fastp auto-enables on
    /// `NextSeq` / `NovaSeq` instruments — we don't auto-detect yet).
    #[arg(short = 'g', long = "trim_poly_g")]
    trim_poly_g: bool,

    /// Minimum run length for poly-G detection (default 10, matches fastp).
    #[arg(long = "poly_g_min_len", default_value_t = 10)]
    poly_g_min_len: usize,

    /// Force poly-X trim. Detects the dominant tail base via the
    /// generalised poly-X scan.
    #[arg(short = 'x', long = "trim_poly_x")]
    trim_poly_x: bool,

    /// Poly-X minimum run length (default 10).
    #[arg(long = "poly_x_min_len", default_value_t = 10)]
    poly_x_min_len: usize,

    /// Bases trimmed from R1 5'. Per-mate.
    /// Short alias matches fastp's `-f`.
    #[arg(short = 'f', long = "trim_front1", default_value_t = 0)]
    trim_front1: usize,

    /// Bases trimmed from R1 3'. Per-mate. No short alias because fastp's
    /// `-t` clashes with this crate family's reserved `-t/--threads`.
    #[arg(long = "trim_tail1", default_value_t = 0)]
    trim_tail1: usize,

    /// Bases trimmed from R2 5'. PE only.
    /// Short alias matches fastp's `-F`.
    #[arg(short = 'F', long = "trim_front2", default_value_t = 0)]
    trim_front2: usize,

    /// Bases trimmed from R2 3'. PE only.
    /// No short alias (T is too easily confused with t).
    #[arg(long = "trim_tail2", default_value_t = 0)]
    trim_tail2: usize,

    /// Enable PE overlap-based adapter detection. Defaults off — caller
    /// opts in. When on, fires before the static-sequence fallback.
    #[arg(short = '2', long = "detect_adapter_for_pe")]
    detect_adapter_for_pe: bool,

    /// Minimum overlap length for the PE detector (default 30).
    #[arg(long = "overlap_len_require", default_value_t = 30)]
    overlap_len_require: usize,

    /// Hard cap on mismatches in the overlap (default 5).
    #[arg(long = "overlap_diff_limit", default_value_t = 5)]
    overlap_diff_limit: usize,

    /// Per-position mismatch cap as a fraction of overlap length
    /// (default 0.20, i.e. 20%).
    #[arg(long = "overlap_diff_percent_limit", default_value_t = 0.20)]
    overlap_diff_percent_limit: f32,

    /// Reads shorter than this after all trim layers are discarded.
    /// Matches fastp's `--length_required` (default 15).
    #[arg(short = 'l', long = "length_required", default_value_t = 15)]
    length_required: usize,

    /// Disable the length filter. Equivalent to `-l 1`: every non-empty
    /// trimmed read is emitted. Mirrors fastp's `-L`.
    #[arg(short = 'L', long = "disable_length_filtering")]
    disable_length_filtering: bool,

    #[command(flatten)]
    common: CommonFlags,
}

fn build_config(cli: &Cli) -> PipelineConfig {
    let adapter1 = if cli.disable_adapter_trimming {
        None
    } else {
        match cli.adapter_sequence.as_deref() {
            Some("") => None,
            Some(s) => Some(AdapterConfig {
                sequence: s.as_bytes().to_vec(),
                min_match_len: 5,
                max_mismatch_rate: 0.2,
            }),
            None => Some(AdapterConfig::illumina_truseq_r1()),
        }
    };
    let adapter2 = if cli.disable_adapter_trimming {
        None
    } else {
        match cli.adapter_sequence_r2.as_deref() {
            Some("") => None,
            Some(s) => Some(AdapterConfig {
                sequence: s.as_bytes().to_vec(),
                min_match_len: 5,
                max_mismatch_rate: 0.2,
            }),
            None => Some(AdapterConfig::illumina_truseq_r2()),
        }
    };

    let poly_g = if cli.trim_poly_g {
        Some(PolyXConfig {
            base: b'G',
            min_len: cli.poly_g_min_len,
            ..PolyXConfig::default()
        })
    } else {
        None
    };
    let poly_x = if cli.trim_poly_x {
        Some(PolyXConfig {
            base: b'G',
            min_len: cli.poly_x_min_len,
            ..PolyXConfig::default()
        })
    } else {
        None
    };

    let overlap = if cli.detect_adapter_for_pe {
        Some(OverlapConfig {
            overlap_require: cli.overlap_len_require,
            diff_limit: cli.overlap_diff_limit,
            diff_percent_limit: cli.overlap_diff_percent_limit,
        })
    } else {
        None
    };

    let min_length_required = if cli.disable_length_filtering {
        1
    } else {
        cli.length_required
    };

    PipelineConfig {
        fixed1: FixedTrimConfig {
            trim_front: cli.trim_front1,
            trim_tail: cli.trim_tail1,
        },
        fixed2: FixedTrimConfig {
            trim_front: cli.trim_front2,
            trim_tail: cli.trim_tail2,
        },
        adapter1,
        adapter2,
        poly_g,
        poly_x,
        overlap,
        min_length_required,
    }
}

fn pipeline(args: &Cli) -> Result<TrimReport> {
    let cfg = build_config(args);
    let p = Pipeline::new(&cfg);

    match (args.in2.as_ref(), args.out2.as_ref()) {
        (Some(in2), Some(out2)) => p.run_pe(&args.in1, in2, &args.out1, out2),
        (None, None) => p.run_se(&args.in1, &args.out1),
        _ => Err(RsomicsError::ConfigError(
            "--in2 and --out2 must be supplied together for PE input".into(),
        )),
    }
}

fn main() -> ExitCode {
    let args = Cli::parse();
    let common = args.common.clone();
    run(&common, META, || pipeline(&args))
}
