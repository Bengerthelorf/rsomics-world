use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_fastq_dedup::{DedupConfig, DedupMode, DedupReport, run_se};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-dedup", disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'i', long = "in1")]
    in1: PathBuf,
    #[arg(short = 'o', long = "out1")]
    out1: PathBuf,
    #[arg(long = "mode", default_value = "bin")]
    mode: String,
    #[arg(short = 'k', long = "kmer-size", default_value_t = 12)]
    kmer_size: usize,
    #[arg(long = "tail-offset", default_value_t = 0)]
    tail_offset: usize,
    #[command(flatten)]
    pub common: CommonFlags,
}

fn parse_mode(s: &str) -> Result<DedupMode> {
    match s {
        "bin" => Ok(DedupMode::KmerBin),
        "full" => Ok(DedupMode::FullSeq),
        other => Err(RsomicsError::ConfigError(format!(
            "unknown --mode {other:?}, expected `bin` or `full`"
        ))),
    }
}

impl Cli {
    pub fn execute(&self) -> Result<DedupReport> {
        let cfg = DedupConfig {
            mode: parse_mode(&self.mode)?,
            k: self.kmer_size,
            tail_offset: self.tail_offset,
        };
        run_se(&self.in1, &self.out1, cfg)
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Sequence-based FASTQ deduplication (fastp -D / seqkit rmdup compat).",
    origin: Some(Origin {
        upstream: "fastp + seqkit",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/bty560"),
    }),
    usage_lines: &["[OPTIONS] -i <FASTQ> -o <FASTQ>"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "in1",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Input FASTQ (gz/bz2/xz/zst autodetect)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "out1",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Output FASTQ (.gz writes gzipped)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "mode",
                aliases: &[],
                value: Some("<bin|full>"),
                type_hint: Some("enum"),
                required: false,
                default: Some("bin"),
                description: "bin = fastp -D hash-bin (fast, lossy); full = seqkit rmdup full-seq (exact)",
                why_default: None,
            },
            FlagSpec {
                short: Some('k'),
                long: "kmer-size",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("12"),
                description: "k-mer length for hash-bin mode (fastp default 12)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "tail-offset",
                aliases: &[],
                value: Some("<n>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("0"),
                description: "Offset from 3' end to place the fingerprint window (0 = last k bases)",
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
    }],
    examples: &[
        Example {
            description: "Fast hash-bin dedup (fastp -D)",
            command: "rsomics-fastq-dedup -i in.fq.gz -o out.fq.gz",
        },
        Example {
            description: "Exact full-sequence dedup (seqkit rmdup -s)",
            command: "rsomics-fastq-dedup -i in.fq -o out.fq --mode full",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-dedup/0.1/#json-output-schema"),
};
