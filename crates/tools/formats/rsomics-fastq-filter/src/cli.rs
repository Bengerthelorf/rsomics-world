use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_fastq_filter::{FilterConfig, FilterReport, Pipeline};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str = "FASTQ per-read quality + length filter (pass/fail whole reads; per-function partition of fastp).";

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-filter", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
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

    /// Minimum Phred quality score for a base to be considered "qualified".
    /// A base with Phred score strictly less than this threshold is counted as
    /// low-quality. Default 15 (fastp default, `qualifiedQual = '0'` in ASCII = Q15).
    #[arg(
        long = "qualified_quality_phred",
        alias = "qualified-quality-phred",
        default_value_t = 15
    )]
    qualified_quality_phred: u16,

    /// Maximum percentage of low-quality bases permitted per read.
    /// A read fails iff `low_qual_count * 100 > unqualified_percent_limit * len`.
    /// Default 40 (fastp default).
    #[arg(
        long = "unqualified_percent_limit",
        alias = "unqualified-percent-limit",
        default_value_t = 40,
        value_parser = clap::value_parser!(u8).range(0..=100),
    )]
    unqualified_percent_limit: u8,

    /// Reads with more than this many N bases are discarded.
    /// Check: `n_count > n_base_limit` (strict greater-than). Default 5 (fastp default).
    #[arg(long = "n_base_limit", alias = "n-base-limit", default_value_t = 5)]
    n_base_limit: usize,

    /// Minimum read length. Reads shorter than this are discarded.
    /// Default 15 (fastp default).
    #[arg(
        short = 'l',
        long = "length_required",
        alias = "length-required",
        default_value_t = 15
    )]
    length_required: usize,

    /// Maximum read length. Reads longer than this are discarded. 0 = no upper bound.
    /// Default 0 (fastp default).
    #[arg(long = "length_limit", alias = "length-limit", default_value_t = 0)]
    length_limit: usize,

    /// Disable quality filtering entirely (overrides `--qualified_quality_phred` /
    /// `--unqualified_percent_limit` / `--n_base_limit`).
    #[arg(
        short = 'Q',
        long = "disable_quality_filtering",
        alias = "disable-quality-filtering"
    )]
    disable_quality_filtering: bool,

    /// Disable length filtering entirely (overrides `-l` / `--length_limit`).
    #[arg(
        short = 'L',
        long = "disable_length_filtering",
        alias = "disable-length-filtering"
    )]
    disable_length_filtering: bool,

    /// libdeflate gzip compression level for `.gz` output. Default 4 (fastp default).
    /// 1 = fastest/largest, 12 = slowest/smallest.
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
    fn build_config(&self) -> FilterConfig {
        FilterConfig {
            qual_threshold_ascii: self.qualified_quality_phred + 33,
            unqualified_percent_limit: self.unqualified_percent_limit,
            n_base_limit: self.n_base_limit,
            required_length: self.length_required,
            max_length: self.length_limit,
            quality_enabled: !self.disable_quality_filtering,
            length_enabled: !self.disable_length_filtering,
        }
    }

    pub fn execute(&self) -> Result<FilterReport> {
        let cfg = self.build_config();
        let p = Pipeline::new(&cfg, self.compression);

        let report = match (self.in2.as_ref(), self.out2.as_ref()) {
            (Some(in2), Some(out2)) => p.run_pe(&self.in1, in2, &self.out1, out2)?,
            (None, None) => p.run_se(&self.in1, &self.out1)?,
            _ => {
                return Err(RsomicsError::ConfigError(
                    "--in2 and --out2 must be supplied together for PE mode".into(),
                ));
            }
        };

        if !self.common.json {
            let dropped = report.reads_in - report.reads_out;
            if dropped > 0 {
                eprintln!("warning: {dropped} reads dropped (failed quality/length filter)");
            }
        }

        Ok(report)
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
            title: "QUALITY FILTER",
            flags: &[
                FlagSpec {
                    short: None,
                    long: "qualified_quality_phred",
                    aliases: &["qualified-quality-phred"],
                    value: Some("<n>"),
                    type_hint: Some("u16"),
                    required: false,
                    default: Some("15"),
                    description: "Phred threshold; bases with Phred < this are low-quality",
                    why_default: Some("fastp default (qualifiedQual='0' ASCII = Q15 Phred+33)"),
                },
                FlagSpec {
                    short: None,
                    long: "unqualified_percent_limit",
                    aliases: &["unqualified-percent-limit"],
                    value: Some("<pct>"),
                    type_hint: Some("u8"),
                    required: false,
                    default: Some("40"),
                    description: "Max % low-qual bases; read fails if low_qual_count > limit*len/100",
                    why_default: Some("fastp default"),
                },
                FlagSpec {
                    short: None,
                    long: "n_base_limit",
                    aliases: &["n-base-limit"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("5"),
                    description: "Max N bases per read (strict >); both N and n count",
                    why_default: Some("fastp default"),
                },
                FlagSpec {
                    short: Some('Q'),
                    long: "disable_quality_filtering",
                    aliases: &["disable-quality-filtering"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Disable quality filter entirely",
                    why_default: None,
                },
            ],
        },
        Section {
            title: "LENGTH FILTER",
            flags: &[
                FlagSpec {
                    short: Some('l'),
                    long: "length_required",
                    aliases: &["length-required"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("15"),
                    description: "Minimum read length; reads shorter are discarded",
                    why_default: Some("fastp default"),
                },
                FlagSpec {
                    short: None,
                    long: "length_limit",
                    aliases: &["length-limit"],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("0"),
                    description: "Maximum read length; 0 = no upper bound",
                    why_default: Some("fastp default — maxLength=0"),
                },
                FlagSpec {
                    short: Some('L'),
                    long: "disable_length_filtering",
                    aliases: &["disable-length-filtering"],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Disable length filter entirely",
                    why_default: None,
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
            description: "SE quality + length filter with fastp defaults",
            command: "rsomics-fastq-filter -i in.fq.gz -o out.fq.gz",
        },
        Example {
            description: "Strict filter: Q20 threshold, max 20% low-qual, min 50 bp",
            command: "rsomics-fastq-filter -i in.fq.gz -o out.fq.gz --qualified_quality_phred 20 --unqualified_percent_limit 20 -l 50",
        },
        Example {
            description: "PE mode: pair dropped if either mate fails, parallel libdeflate output",
            command: "rsomics-fastq-filter -i r1.fq.gz -I r2.fq.gz -o r1.filt.fq.gz -O r2.filt.fq.gz",
        },
        Example {
            description: "Length-only filter (disable quality), JSON report",
            command: "rsomics-fastq-filter -i in.fq.gz -o out.fq.gz -Q -l 30 --json | jq .result",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-filter/0.1/#json-output-schema"),
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
