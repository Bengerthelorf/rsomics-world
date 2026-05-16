use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_fastq_split::{Pipeline, SplitConfig, SplitMode, SplitReport};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str =
    "Split a FASTQ into N files or by line count (per-function partition of fastp).";

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-split", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// R1 input. `.fq` / `.fq.gz` autodetected by magic bytes.
    #[arg(short = 'i', long = "in1", alias = "in-1")]
    in1: PathBuf,

    /// R1 output base. Split files are `<digits>.<basename>` in its directory
    /// (e.g. `out.fq.gz` → `0001.out.fq.gz`). `.gz` keeps per-file compression.
    #[arg(short = 'o', long = "out1", alias = "out-1")]
    out1: PathBuf,

    /// R2 input (PE mode).
    #[arg(short = 'I', long = "in2", alias = "in-2")]
    in2: Option<PathBuf>,

    /// R2 output base (PE mode).
    #[arg(short = 'O', long = "out2", alias = "out-2")]
    out2: Option<PathBuf>,

    /// Split into exactly N files (exact read count: `ceil(total/N)` reads
    /// each). Mutually exclusive with `--split_by_lines`.
    #[arg(short = 's', long = "split", alias = "split-number")]
    split: Option<usize>,

    /// Split so each file holds L lines (a FASTQ record is 4 lines, so L must
    /// be a multiple of 4). Mutually exclusive with `--split`.
    #[arg(long = "split_by_lines", alias = "split-by-lines")]
    split_by_lines: Option<usize>,

    /// Zero-pad width of the numeric file prefix (fastp `--split_prefix_digits`).
    #[arg(
        long = "split_prefix_digits",
        alias = "split-prefix-digits",
        default_value_t = 4
    )]
    split_prefix_digits: usize,

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
    fn build_config(&self) -> Result<SplitConfig> {
        let mode = match (self.split, self.split_by_lines) {
            (Some(_), Some(_)) => {
                return Err(RsomicsError::ConfigError(
                    "--split and --split_by_lines are mutually exclusive".into(),
                ));
            }
            (Some(n), None) => SplitMode::ByNumber(n),
            (None, Some(l)) => SplitMode::ByLines(l),
            (None, None) => {
                return Err(RsomicsError::ConfigError(
                    "one of --split <N> or --split_by_lines <L> is required".into(),
                ));
            }
        };
        if self.split_prefix_digits == 0 {
            return Err(RsomicsError::ConfigError(
                "--split_prefix_digits must be > 0".into(),
            ));
        }
        Ok(SplitConfig {
            mode,
            digits: self.split_prefix_digits,
        })
    }

    pub fn execute(&self) -> Result<SplitReport> {
        let cfg = self.build_config()?;
        let p = Pipeline::new(&cfg, self.compression);
        match (self.in2.as_ref(), self.out2.as_ref()) {
            (Some(in2), Some(out2)) => p.run_pe(&self.in1, in2, &self.out1, out2),
            (None, None) => p.run_se(&self.in1, &self.out1),
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
        "--split <N> --in1 <PATH> --out1 <PATH>",
        "--split_by_lines <L> --in1 <R1> --in2 <R2> --out1 <O1> --out2 <O2>   (PE)",
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
                    description: "R1 output base → <digits>.<basename> per file",
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
                    description: "R2 output base (PE mode)",
                    why_default: None,
                },
            ],
        },
        Section {
            title: "SPLIT MODE (exactly one required)",
            flags: &[
                FlagSpec {
                    short: Some('s'),
                    long: "split",
                    aliases: &["split-number"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: None,
                    description: "Into exactly N files (exact ceil(total/N) per file)",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "split_by_lines",
                    aliases: &["split-by-lines"],
                    value: Some("<l>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: None,
                    description: "L lines per file (multiple of 4); byte-equal to fastp 0.20.1",
                    why_default: None,
                },
                FlagSpec {
                    short: None,
                    long: "split_prefix_digits",
                    aliases: &["split-prefix-digits"],
                    value: Some("<d>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("4"),
                    description: "Zero-pad width of the numeric file prefix",
                    why_default: Some("fastp default"),
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
                    why_default: Some("fastp default"),
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
            description: "Split into 8 files (exact count) for parallel downstream",
            command: "rsomics-fastq-split -i in.fq.gz -o out.fq.gz --split 8",
        },
        Example {
            description: "4000 lines (= 1000 reads) per file, fastp-0.20.1-equal",
            command: "rsomics-fastq-split -i in.fq.gz -o out.fq.gz --split_by_lines 4000",
        },
        Example {
            description: "PE split in lockstep, JSON report",
            command: "rsomics-fastq-split -i r1.fq.gz -I r2.fq.gz -o o1.fq.gz -O o2.fq.gz --split 16 --json | jq .result",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-split/0.1/#json-output-schema"),
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
