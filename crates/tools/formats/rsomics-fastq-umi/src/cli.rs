use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_fastq_umi::{Pipeline, UmiConfig, UmiLoc, UmiReport};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str = "FASTQ inline-UMI extract + stamp: move the 5' UMI into the read name (per-function partition of fastp).";

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-umi", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
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

    /// UMI location (fastp `--umi_loc`): `read1` / `read2` (5' of that read's
    /// sequence), `index1` / `index2` (read-name index field), `per_index`,
    /// `per_read`. `read2` and `index2` require PE input.
    #[arg(long = "umi_loc", alias = "umi-loc", default_value = "read1")]
    umi_loc: String,

    /// UMI length in bases (fastp `--umi_len`); required for the sequence
    /// locations (`read1`/`read2`/`per_read`), unused by the index locations.
    #[arg(long = "umi_len", alias = "umi-len")]
    umi_len: Option<usize>,

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
    pub common: CommonFlags,
}

impl Cli {
    fn build_config(&self) -> Result<UmiConfig> {
        let loc = match self.umi_loc.as_str() {
            "read1" => UmiLoc::Read1,
            "read2" => UmiLoc::Read2,
            "index1" => UmiLoc::Index1,
            "index2" => UmiLoc::Index2,
            "per_index" | "per-index" => UmiLoc::PerIndex,
            "per_read" | "per-read" => UmiLoc::PerRead,
            other => {
                return Err(RsomicsError::ConfigError(format!(
                    "--umi_loc must be one of read1|read2|index1|index2|per_index|per_read, got {other:?}"
                )));
            }
        };
        let seq_loc = matches!(loc, UmiLoc::Read1 | UmiLoc::Read2 | UmiLoc::PerRead);
        let len = match self.umi_len {
            Some(n) if n > 0 => n,
            _ if seq_loc => {
                return Err(RsomicsError::ConfigError(
                    "--umi_len > 0 is required for --umi_loc read1/read2/per_read".into(),
                ));
            }
            _ => 0,
        };
        let delim = self.umi_delim.as_bytes();
        if delim.len() != 1 {
            return Err(RsomicsError::ConfigError(
                "--umi_delim must be a single byte".into(),
            ));
        }
        Ok(UmiConfig {
            loc,
            len,
            skip: self.umi_skip,
            prefix: self.umi_prefix.as_bytes().to_vec(),
            delimiter: delim[0],
        })
    }

    pub fn execute(&self) -> Result<UmiReport> {
        let cfg = self.build_config()?;
        let p = Pipeline::new(&cfg, self.compression);
        match (self.in2.as_ref(), self.out2.as_ref()) {
            (Some(in2), Some(out2)) => p.run_pe(&self.in1, in2, &self.out1, out2),
            (None, None) => {
                if matches!(
                    cfg.loc,
                    UmiLoc::Read2 | UmiLoc::Index2 | UmiLoc::PerIndex | UmiLoc::PerRead
                ) {
                    return Err(RsomicsError::ConfigError(
                        "--umi_loc read2/index2/per_index/per_read requires PE input (--in2/--out2)"
                            .into(),
                    ));
                }
                p.run_se(&self.in1, &self.out1)
            }
            _ => Err(RsomicsError::ConfigError(
                "--in2 and --out2 must be supplied together for PE mode".into(),
            )),
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
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
        "[OPTIONS] --umi_loc index1 --in1 <PATH> --out1 <PATH>   (read-name index, no --umi_len)",
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
                    value: Some("<loc>"),
                    type_hint: Some("String"),
                    required: false,
                    default: Some("read1"),
                    description: "read1|read2|index1|index2|per_index|per_read (read2/index2 need PE)",
                    why_default: Some("read1 — the dominant inline-UMI layout"),
                },
                FlagSpec {
                    short: None,
                    long: "umi_len",
                    aliases: &["umi-len"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: None,
                    description: "UMI length in bases (required for read1/read2/per_read; clamped to read length)",
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
#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    /// clap's `debug_assert` validates the whole arg graph (unique shorts
    /// incl. the flattened `CommonFlags`, no id clashes). It only fires
    /// when the binary parses, so without this test a CLI-definition error
    /// is invisible to `cargo test` and lib unit tests.
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
