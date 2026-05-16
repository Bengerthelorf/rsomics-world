use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta, run};
use rsomics_fastq_umi::{Pipeline, UmiConfig, UmiLoc, UmiReport};
use rsomics_help::{
    Example, FlagSpec, HelpSpec, Origin, Section, intercept_help, render as render_help,
};

const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str = "FASTQ inline-UMI extract + stamp: move the 5' UMI into the read name (per-function partition of fastp).";

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-umi", version, about, long_about = None, disable_help_flag = true)]
struct Cli {
    /// R1 input. `.fq` / `.fq.gz` autodetected by magic bytes.
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

    /// Which read carries the inline UMI: `read1` or `read2` (fastp
    /// `--umi_loc`). `read2` requires PE. Index-based locations are not
    /// implemented (no consumer needs them yet).
    #[arg(long = "umi_loc", alias = "umi-loc", default_value = "read1")]
    umi_loc: String,

    /// UMI length in bases (fastp `--umi_len`). Clamped to the read length.
    #[arg(long = "umi_len", alias = "umi-len")]
    umi_len: usize,

    /// Extra 5' bases removed after the UMI (fastp `--umi_skip`). Default 0.
    #[arg(long = "umi_skip", alias = "umi-skip", default_value_t = 0)]
    umi_skip: usize,

    /// UMI name prefix (fastp `--umi_prefix`); joined to the UMI with `_`.
    #[arg(long = "umi_prefix", alias = "umi-prefix", default_value = "")]
    umi_prefix: String,

    /// Read-name / UMI separator. fastp default `:`.
    #[arg(long = "umi_delim", alias = "umi-delim", default_value = ":")]
    umi_delim: String,

    /// libdeflate gzip compression level for `.gz` output. Default 4 (fastp default).
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

fn build_config(cli: &Cli) -> Result<UmiConfig> {
    let loc = match cli.umi_loc.as_str() {
        "read1" => UmiLoc::Read1,
        "read2" => UmiLoc::Read2,
        other => {
            return Err(RsomicsError::ConfigError(format!(
                "--umi_loc must be read1 or read2 (index-based locations are not implemented), got {other:?}"
            )));
        }
    };
    if cli.umi_len == 0 {
        return Err(RsomicsError::ConfigError("--umi_len must be > 0".into()));
    }
    let delim = cli.umi_delim.as_bytes();
    if delim.len() != 1 {
        return Err(RsomicsError::ConfigError(
            "--umi_delim must be a single byte".into(),
        ));
    }
    Ok(UmiConfig {
        loc,
        len: cli.umi_len,
        skip: cli.umi_skip,
        prefix: cli.umi_prefix.as_bytes().to_vec(),
        delimiter: delim[0],
    })
}

fn pipeline(cli: &Cli) -> Result<UmiReport> {
    let cfg = build_config(cli)?;
    let p = Pipeline::new(&cfg, cli.compression);
    match (cli.in2.as_ref(), cli.out2.as_ref()) {
        (Some(in2), Some(out2)) => p.run_pe(&cli.in1, in2, &cli.out1, out2),
        (None, None) => {
            if cfg.loc == UmiLoc::Read2 {
                return Err(RsomicsError::ConfigError(
                    "--umi_loc read2 requires PE input (--in2/--out2)".into(),
                ));
            }
            p.run_se(&cli.in1, &cli.out1)
        }
        _ => Err(RsomicsError::ConfigError(
            "--in2 and --out2 must be supplied together for PE mode".into(),
        )),
    }
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
        "[OPTIONS] --umi_len <N> --in1 <PATH> --out1 <PATH>",
        "[OPTIONS] --umi_len <N> --in1 <R1> --in2 <R2> --out1 <O1> --out2 <O2>   (PE)",
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
                    description: "R1 input (gz autodetect by magic bytes)",
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
            title: "UMI",
            flags: &[
                FlagSpec {
                    short: None,
                    long: "umi_loc",
                    aliases: &["umi-loc"],
                    value: Some("<read1|read2>"),
                    type_hint: Some("String"),
                    required: false,
                    default: Some("read1"),
                    description: "Which read carries the inline 5' UMI",
                    why_default: Some("read1 — the dominant inline-UMI layout"),
                },
                FlagSpec {
                    short: None,
                    long: "umi_len",
                    aliases: &["umi-len"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: true,
                    default: None,
                    description: "UMI length in bases (clamped to read length)",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "umi_skip",
                    aliases: &["umi-skip"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("0"),
                    description: "Extra 5' bases removed after the UMI",
                    why_default: Some("fastp default"),
                },
                FlagSpec {
                    short: None,
                    long: "umi_prefix",
                    aliases: &["umi-prefix"],
                    value: Some("<s>"),
                    type_hint: Some("String"),
                    required: false,
                    default: Some("\"\""),
                    description: "Name prefix; joined to the UMI with `_`",
                    why_default: Some("fastp default — no prefix"),
                },
                FlagSpec {
                    short: None,
                    long: "umi_delim",
                    aliases: &["umi-delim"],
                    value: Some("<c>"),
                    type_hint: Some("String"),
                    required: false,
                    default: Some(":"),
                    description: "Read-name / UMI separator (single byte)",
                    why_default: Some("fastp default delimiter"),
                },
            ],
        },
        Section {
            title: "OUTPUT",
            flags: &[
                FlagSpec {
                    short: None,
                    long: "compression",
                    aliases: &["compression-level"],
                    value: Some("<lvl>"),
                    type_hint: Some("i32"),
                    required: false,
                    default: Some("4"),
                    description: "libdeflate gz compression level 1-12 for .gz output",
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
            description: "SE: extract an 8 bp inline UMI from R1 into the read name",
            command: "rsomics-fastq-umi -i in.fq.gz -o out.fq.gz --umi_len 8",
        },
        Example {
            description: "PE: UMI on R1, also skip 1 linker base, prefix the tag",
            command: "rsomics-fastq-umi -i r1.fq.gz -I r2.fq.gz -o o1.fq.gz -O o2.fq.gz --umi_len 8 --umi_skip 1 --umi_prefix UMI",
        },
        Example {
            description: "PE: UMI carried on R2, JSON report",
            command: "rsomics-fastq-umi -i r1.fq.gz -I r2.fq.gz -o o1.fq.gz -O o2.fq.gz --umi_loc read2 --umi_len 10 --json | jq .result",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-umi/0.1/#json-output-schema"),
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
