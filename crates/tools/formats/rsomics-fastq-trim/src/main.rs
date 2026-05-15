use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta, run};
use rsomics_help::{
    Example, FlagSpec, HelpSpec, Origin, Section, intercept_help, render as render_help,
};

use rsomics_fastq_trim::{
    AdapterConfig, FixedTrimConfig, OverlapConfig, Pipeline, PipelineConfig, PolyXConfig,
    TrimReport,
};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str =
    "FASTQ adapter / polyG / polyX / fixed-length trimming (per-function partition of fastp).";

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-trim", version, about, long_about = None, disable_help_flag = true)]
#[allow(clippy::struct_excessive_bools)]
struct Cli {
    /// R1 input. `.fq` / `.fq.gz` / `.fq.bz2` / `.fq.xz` / `.fq.zst` autodetected.
    #[arg(short = 'i', long = "in1", alias = "in-1")]
    in1: PathBuf,

    /// R1 output. `.gz` suffix triggers parallel libdeflate compression.
    #[arg(short = 'o', long = "out1", alias = "out-1")]
    out1: PathBuf,

    /// R2 input (PE mode).
    #[arg(short = 'I', long = "in2", alias = "in-2")]
    in2: Option<PathBuf>,

    /// R2 output (PE mode).
    #[arg(short = 'O', long = "out2", alias = "out-2")]
    out2: Option<PathBuf>,

    /// R1 adapter sequence. Default = Illumina `TruSeq` R1 prefix
    /// (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`); matches fastp's default
    /// and covers the majority of Illumina library kits. Pass an empty
    /// string to disable adapter trim.
    #[arg(short = 'a', long = "adapter_sequence", alias = "adapter-sequence")]
    adapter_sequence: Option<String>,

    /// R2 adapter sequence. Default = Illumina `TruSeq` R2 prefix
    /// (`AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`). PE only.
    #[arg(long = "adapter_sequence_r2", alias = "adapter-sequence-r2")]
    adapter_sequence_r2: Option<String>,

    /// Minimum bases of compared adapter prefix required for a match.
    /// Default 5 (fastp default). Raise to be stricter, lower to be
    /// more aggressive on short reads.
    #[arg(
        long = "adapter_min_len",
        alias = "adapter-min-len",
        default_value_t = 5
    )]
    adapter_min_len: usize,

    /// Maximum mismatch rate across the compared adapter region.
    /// Default 0.20 (fastp default, allows 1 mismatch per 5 bases).
    #[arg(
        long = "adapter_max_mismatch_rate",
        alias = "adapter-max-mismatch-rate",
        default_value_t = 0.20
    )]
    adapter_max_mismatch_rate: f32,

    /// Disable static-sequence adapter trim. Useful when relying purely
    /// on PE overlap detection or running upstream-clean data.
    #[arg(
        short = 'A',
        long = "disable_adapter_trimming",
        alias = "disable-adapter-trimming"
    )]
    disable_adapter_trimming: bool,

    /// Force poly-G trim. fastp auto-enables on `NextSeq` / `NovaSeq`
    /// 2-color-chemistry instruments where dark cycles read as G; we
    /// don't auto-detect from FASTQ headers — pass explicitly.
    #[arg(short = 'g', long = "trim_poly_g", alias = "trim-poly-g")]
    trim_poly_g: bool,

    /// Base trimmed by `--trim_poly_g`. Default `G`. Override only if
    /// you need a poly-A / poly-T / poly-C single-base scan with the
    /// same parameters; for auto-detect use `--trim_poly_x` instead.
    #[arg(
        long = "poly_g_base",
        alias = "poly-g-base",
        default_value = "G",
        value_parser = parse_base
    )]
    poly_g_base: u8,

    /// Minimum poly-G run length. Default 10 (fastp default — empirical
    /// false-positive rate is acceptable below 1% on real `NovaSeq`
    /// libraries; raise for short-amplicon protocols).
    #[arg(
        long = "poly_g_min_len",
        alias = "poly-g-min-len",
        default_value_t = 10
    )]
    poly_g_min_len: usize,

    /// Trim the 3' poly-X tail by dominant-base detection. Counts
    /// A/C/G/T simultaneously across the tail and trims at the last
    /// occurrence of the most-represented base.
    #[arg(short = 'x', long = "trim_poly_x", alias = "trim-poly-x")]
    trim_poly_x: bool,

    /// Poly-X minimum run length. Default 10 (fastp default).
    #[arg(
        long = "poly_x_min_len",
        alias = "poly-x-min-len",
        default_value_t = 10
    )]
    poly_x_min_len: usize,

    /// Hard cap on mismatches inside a poly-X run regardless of length.
    /// Default 5 (fastp default).
    #[arg(
        long = "polyx_max_mismatches",
        alias = "polyx-max-mismatches",
        default_value_t = 5
    )]
    polyx_max_mismatches: usize,

    /// Rate cap: one allowed mismatch per N scanned bases. Default 8
    /// (fastp default — `floor(scanned / 8)` interspersed non-target
    /// bases tolerated). Must be non-zero.
    #[arg(
        long = "polyx_mismatch_per_bases",
        alias = "polyx-mismatch-per-bases",
        default_value_t = 8
    )]
    polyx_mismatch_per_bases: usize,

    /// Bases trimmed from R1 5'. Short alias matches fastp's `-f`.
    #[arg(
        short = 'f',
        long = "trim_front1",
        alias = "trim-front1",
        default_value_t = 0
    )]
    trim_front1: usize,

    /// Bases trimmed from R1 3'. No short alias because fastp's `-t`
    /// collides with this crate family's reserved `-t/--threads`.
    #[arg(long = "trim_tail1", alias = "trim-tail1", default_value_t = 0)]
    trim_tail1: usize,

    /// Bases trimmed from R2 5'. PE only.
    #[arg(
        short = 'F',
        long = "trim_front2",
        alias = "trim-front2",
        default_value_t = 0
    )]
    trim_front2: usize,

    /// Bases trimmed from R2 3'. PE only.
    #[arg(long = "trim_tail2", alias = "trim-tail2", default_value_t = 0)]
    trim_tail2: usize,

    /// Enable PE overlap-based adapter detection. Off by default — when
    /// on, the geometry of R1 vs reverse-complemented R2 is used to
    /// find the adapter cut-point. fastp recommends turning this on for
    /// "ultra-clean" data; static-seq fallback still fires when overlap
    /// is not detected.
    #[arg(
        short = '2',
        long = "detect_adapter_for_pe",
        alias = "detect-adapter-for-pe"
    )]
    detect_adapter_for_pe: bool,

    /// Minimum overlap length for the PE detector. Default 30 (fastp
    /// default — below this the overlap is statistically likely to be
    /// random chance, not a real adapter).
    #[arg(
        long = "overlap_len_require",
        alias = "overlap-len-require",
        default_value_t = 30
    )]
    overlap_len_require: usize,

    /// Hard cap on mismatches inside the PE overlap. Default 5
    /// (fastp default).
    #[arg(
        long = "overlap_diff_limit",
        alias = "overlap-diff-limit",
        default_value_t = 5
    )]
    overlap_diff_limit: usize,

    /// Per-position mismatch cap as a fraction of overlap length.
    /// Default 0.20 (fastp default = 20%). Clamped to `[0.0, 1.0]`.
    #[arg(
        long = "overlap_diff_percent_limit",
        alias = "overlap-diff-percent-limit",
        default_value_t = 0.20
    )]
    overlap_diff_percent_limit: f32,

    /// Reads shorter than this after all trim layers are discarded.
    /// Default 15 (fastp default — tuned for 150 bp WGS; amplicon /
    /// miRNA protocols need a lower value).
    #[arg(
        short = 'l',
        long = "length_required",
        alias = "length-required",
        default_value_t = 15
    )]
    length_required: usize,

    /// Disable the length filter. Equivalent to `-l 1`: every non-empty
    /// trimmed read is emitted. Mirrors fastp's `-L`.
    #[arg(
        short = 'L',
        long = "disable_length_filtering",
        alias = "disable-length-filtering"
    )]
    disable_length_filtering: bool,

    /// libdeflate gzip compression level for `.gz` output. Default 4
    /// (fastp default — best ratio/speed trade-off). 1 = fastest /
    /// largest, 12 = slowest / smallest.
    #[arg(
        long = "compression",
        alias = "compression-level",
        default_value_t = 4,
        value_parser = clap::value_parser!(i32).range(1..=12),
    )]
    compression: i32,

    #[command(flatten)]
    common: CommonFlags,
}

fn parse_base(s: &str) -> std::result::Result<u8, String> {
    let bytes = s.as_bytes();
    if bytes.len() != 1 {
        return Err(format!("expected a single character, got {s:?}"));
    }
    let b = bytes[0].to_ascii_uppercase();
    // N is rejected — the forced-base poly-X scanner has no meaningful
    // semantics when the target base is the ambiguity placeholder.
    if matches!(b, b'A' | b'C' | b'G' | b'T') {
        Ok(b)
    } else {
        Err(format!("expected one of A C G T, got {s:?}"))
    }
}

fn build_config(cli: &Cli) -> Result<PipelineConfig> {
    let adapter_with = |s: &str| AdapterConfig {
        sequence: s.as_bytes().to_vec(),
        min_match_len: cli.adapter_min_len,
        max_mismatch_rate: cli.adapter_max_mismatch_rate,
    };
    let adapter1 = if cli.disable_adapter_trimming {
        None
    } else {
        match cli.adapter_sequence.as_deref() {
            Some("") => None,
            Some(s) => Some(adapter_with(s)),
            None => Some(AdapterConfig {
                min_match_len: cli.adapter_min_len,
                max_mismatch_rate: cli.adapter_max_mismatch_rate,
                ..AdapterConfig::illumina_truseq_r1()
            }),
        }
    };
    let adapter2 = if cli.disable_adapter_trimming {
        None
    } else {
        match cli.adapter_sequence_r2.as_deref() {
            Some("") => None,
            Some(s) => Some(adapter_with(s)),
            None => Some(AdapterConfig {
                min_match_len: cli.adapter_min_len,
                max_mismatch_rate: cli.adapter_max_mismatch_rate,
                ..AdapterConfig::illumina_truseq_r2()
            }),
        }
    };

    let mismatch_per_bases =
        std::num::NonZeroUsize::new(cli.polyx_mismatch_per_bases).ok_or_else(|| {
            RsomicsError::ConfigError("--polyx_mismatch_per_bases must be > 0".into())
        })?;
    let poly_g = if cli.trim_poly_g {
        Some(PolyXConfig {
            base: cli.poly_g_base,
            min_len: cli.poly_g_min_len,
            max_mismatches: cli.polyx_max_mismatches,
            mismatch_per_bases,
        })
    } else {
        None
    };
    let poly_x = if cli.trim_poly_x {
        Some(PolyXConfig {
            min_len: cli.poly_x_min_len,
            max_mismatches: cli.polyx_max_mismatches,
            mismatch_per_bases,
            ..PolyXConfig::default()
        })
    } else {
        None
    };

    let overlap = if cli.detect_adapter_for_pe {
        Some(OverlapConfig::sanitised(
            cli.overlap_len_require,
            cli.overlap_diff_limit,
            cli.overlap_diff_percent_limit,
        ))
    } else {
        None
    };

    let min_length_required = if cli.disable_length_filtering {
        1
    } else {
        cli.length_required
    };

    Ok(PipelineConfig {
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
        compression: cli.compression,
    })
}

fn pipeline(args: &Cli) -> Result<TrimReport> {
    let cfg = build_config(args)?;
    let p = Pipeline::new(&cfg);

    let report = match (args.in2.as_ref(), args.out2.as_ref()) {
        (Some(in2), Some(out2)) => p.run_pe(&args.in1, in2, &args.out1, out2)?,
        (None, None) => p.run_se(&args.in1, &args.out1)?,
        _ => {
            return Err(RsomicsError::ConfigError(
                "--in2 and --out2 must be supplied together for PE input".into(),
            ));
        }
    };

    if !args.common.json && report.reads_too_short_after_trim > 0 {
        eprintln!(
            "warning: {} reads dropped (too short after trim; -l adjusts the threshold, -L disables)",
            report.reads_too_short_after_trim
        );
    }

    Ok(report)
}

const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: TAGLINE,
    origin: Some(Origin {
        upstream: "fastp",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/bty560"),
    }),
    usage_lines: &[
        "[OPTIONS] --in1 <PATH> --out1 <PATH>",
        "[OPTIONS] --in1 <R1> --in2 <R2> --out1 <O1> --out2 <O2>   (PE)",
    ],
    sections: &[
        Section {
            title: "INPUT / OUTPUT",
            flags: &[
                FlagSpec {
                    short: Some('i'),
                    long: "in1",
                    aliases: &["in-1"],
                    value: Some("<path>"),
                    type_hint: Some("PathBuf"),
                    required: true,
                    default: None,
                    description: "R1 input (gz/bz2/xz/zst autodetect)",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('o'),
                    long: "out1",
                    aliases: &["out-1"],
                    value: Some("<path>"),
                    type_hint: Some("PathBuf"),
                    required: true,
                    default: None,
                    description: "R1 output (.gz uses parallel libdeflate)",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('I'),
                    long: "in2",
                    aliases: &["in-2"],
                    value: Some("<path>"),
                    type_hint: Some("Option<PathBuf>"),
                    required: false,
                    default: None,
                    description: "R2 input (PE mode)",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('O'),
                    long: "out2",
                    aliases: &["out-2"],
                    value: Some("<path>"),
                    type_hint: Some("Option<PathBuf>"),
                    required: false,
                    default: None,
                    description: "R2 output (PE mode)",
                    why_default: None,
                },
            ],
        },
        Section {
            title: "ADAPTER TRIM",
            flags: &[
                FlagSpec {
                    short: Some('a'),
                    long: "adapter_sequence",
                    aliases: &["adapter-sequence"],
                    value: Some("<seq>"),
                    type_hint: Some("Option<String>"),
                    required: false,
                    default: Some("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"),
                    description: "R1 adapter sequence",
                    why_default: Some("Illumina TruSeq R1 prefix; matches fastp's default"),
                },
                FlagSpec {
                    short: None,
                    long: "adapter_sequence_r2",
                    aliases: &["adapter-sequence-r2"],
                    value: Some("<seq>"),
                    type_hint: Some("Option<String>"),
                    required: false,
                    default: None,
                    description: "R2 adapter sequence (PE)",
                    why_default: Some("falls back to TruSeq R2 prefix if unspecified"),
                },
                FlagSpec {
                    short: None,
                    long: "adapter_min_len",
                    aliases: &["adapter-min-len"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("5"),
                    description: "Min adapter prefix match length",
                    why_default: Some("matches fastp's default"),
                },
                FlagSpec {
                    short: None,
                    long: "adapter_max_mismatch_rate",
                    aliases: &["adapter-max-mismatch-rate"],
                    value: Some("<f>"),
                    type_hint: Some("f32"),
                    required: false,
                    default: Some("0.20"),
                    description: "Max mismatch fraction across compared region",
                    why_default: Some("fastp default — 1 mismatch per 5 bases"),
                },
                FlagSpec {
                    short: Some('A'),
                    long: "disable_adapter_trimming",
                    aliases: &["disable-adapter-trimming"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Skip static-seq adapter trim",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('2'),
                    long: "detect_adapter_for_pe",
                    aliases: &["detect-adapter-for-pe"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Enable PE overlap-based adapter detection",
                    why_default: Some("static-seq fallback still fires when overlap not found"),
                },
            ],
        },
        Section {
            title: "POLY-G / POLY-X",
            flags: &[
                FlagSpec {
                    short: Some('g'),
                    long: "trim_poly_g",
                    aliases: &["trim-poly-g"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Force poly-G trim",
                    why_default: Some(
                        "fastp auto-enables on NextSeq/NovaSeq; we don't auto-detect from headers",
                    ),
                },
                FlagSpec {
                    short: Some('x'),
                    long: "trim_poly_x",
                    aliases: &["trim-poly-x"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Dominant-base poly-X scan (A/C/G/T auto-detect)",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "poly_g_min_len",
                    aliases: &["poly-g-min-len"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("10"),
                    description: "Poly-G min run length",
                    why_default: Some("fastp default"),
                },
                FlagSpec {
                    short: None,
                    long: "polyx_max_mismatches",
                    aliases: &["polyx-max-mismatches"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("5"),
                    description: "Hard cap on mismatches in a poly-X run",
                    why_default: Some("fastp default"),
                },
                FlagSpec {
                    short: None,
                    long: "polyx_mismatch_per_bases",
                    aliases: &["polyx-mismatch-per-bases"],
                    value: Some("<n>"),
                    type_hint: Some("NonZeroUsize"),
                    required: false,
                    default: Some("8"),
                    description: "Rate cap: one mismatch per N scanned bases",
                    why_default: Some(
                        "fastp default — floor(scanned/8) interspersed non-target bases tolerated",
                    ),
                },
            ],
        },
        Section {
            title: "FIXED + OUTPUT",
            flags: &[
                FlagSpec {
                    short: Some('f'),
                    long: "trim_front1",
                    aliases: &["trim-front1"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("0"),
                    description: "Bases trimmed from R1 5' end",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "trim_tail1",
                    aliases: &["trim-tail1"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("0"),
                    description: "Bases trimmed from R1 3' end",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('l'),
                    long: "length_required",
                    aliases: &["length-required"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("15"),
                    description: "Discard reads shorter after trim",
                    why_default: Some(
                        "fastp default — tuned for 150 bp WGS; amplicon/miRNA needs lower",
                    ),
                },
                FlagSpec {
                    short: Some('L'),
                    long: "disable_length_filtering",
                    aliases: &["disable-length-filtering"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Disable length filter (= -l 1)",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "compression",
                    aliases: &["compression-level"],
                    value: Some("<lvl>"),
                    type_hint: Some("i32"),
                    required: false,
                    default: Some("4"),
                    description: "libdeflate gz compression level 1-12",
                    why_default: Some("fastp default — best ratio/speed trade-off"),
                },
                FlagSpec {
                    short: None,
                    long: "json",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "AI-friendly JSON envelope on stdout",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('t'),
                    long: "threads",
                    aliases: &[],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: None,
                    description: "Worker threads (default: available cores)",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('h'),
                    long: "help",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: None,
                    description: "Show this help (add --plain or --json for alt modes)",
                    why_default: None,
                },
            ],
        },
    ],
    examples: &[
        Example {
            description: "SE adapter trim, default adapter",
            command: "rsomics-fastq-trim -i in.fq.gz -o out.fq.gz",
        },
        Example {
            description: "PE with overlap-detect, 8 threads, JSON envelope",
            command: "rsomics-fastq-trim -i r1.fq.gz -I r2.fq.gz -o o1.fq.gz -O o2.fq.gz -2 -t 8 --json | jq .result",
        },
        Example {
            description: "Poly-G trim only (Illumina 2-color chemistry)",
            command: "rsomics-fastq-trim -i in.fq.gz -o out.fq.gz -A -g",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-trim/0.4/#json-output-schema"),
};

fn main() -> ExitCode {
    let raw_args: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&raw_args) {
        render_help(&HELP, mode);
        return ExitCode::SUCCESS;
    }
    let args = Cli::parse();
    let common = args.common.clone();
    run(&common, META, || pipeline(&args))
}
