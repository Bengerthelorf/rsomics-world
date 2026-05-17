use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_fastq_correct::{CorrectConfig, CorrectReport, Pipeline};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

const TAGLINE: &str =
    "FASTQ k-mer-spectrum substitution-error correction (independent Rust port of BFC).";

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-correct", version, about, long_about = None, disable_help_flag = true)]
#[allow(clippy::struct_excessive_bools)]
pub struct Cli {
    /// R1 input. `.fq` / `.fq.gz` autodetected by magic bytes.
    #[arg(short = 'i', long = "in1", alias = "in-1")]
    in1: PathBuf,

    /// R1 output. `.gz` suffix triggers parallel libdeflate compression.
    #[arg(short = 'o', long = "out1", alias = "out-1")]
    out1: PathBuf,

    /// R2 input (PE mode — BFC corrects each read independently; the pair
    /// is not jointly corrected, matching BFC's per-read model).
    #[arg(short = 'I', long = "in2", alias = "in-2")]
    in2: Option<PathBuf>,

    /// R2 output (PE mode).
    #[arg(short = 'O', long = "out2", alias = "out-2")]
    out2: Option<PathBuf>,

    /// K-mer length (BFC `-k`). Must be odd (the canonical-strand pick
    /// uses the middle base) and ≤ 63.
    #[arg(short = 'k', long = "kmer", default_value_t = 33)]
    kmer: usize,

    /// Minimum k-mer coverage to treat a k-mer as trusted/solid
    /// (BFC `-c`).
    #[arg(short = 'c', long = "min-cov", default_value_t = 3)]
    min_cov: i32,

    /// Suppress a second correction within this many bp of the previous
    /// one (BFC `-w`, `win_multi_ec`).
    #[arg(short = 'w', long = "max-ec-window", default_value_t = 10)]
    max_ec_window: i32,

    /// Phred quality at/above which a base is "high-quality" for the
    /// penalty model (BFC `-q`).
    // No short: `-q` is reserved by CommonFlags (`--quiet`) across the
    // rsomics-* family; BFC's `-q` maps to `--qual-threshold` long-only.
    #[arg(long = "qual-threshold", default_value_t = 20)]
    qual_threshold: u8,

    /// Skip correction — pass reads through unchanged (BFC `-E`). The
    /// k-mer count pass still runs (so `--drop-unique-kmer` still works).
    #[arg(short = 'E', long = "no-correct")]
    no_correct: bool,

    /// Discard reads BFC cannot correct (no solid k-mer, uncorrectable
    /// N, too many failures) instead of emitting them unchanged
    /// (BFC `-D`).
    #[arg(short = 'D', long = "discard-uncorrectable")]
    discard_uncorrectable: bool,

    /// Drop any read containing a k-mer seen only once (BFC `-1`).
    #[arg(short = '1', long = "drop-unique-kmer")]
    drop_unique_kmer: bool,

    /// Emit FASTA (drop qualities) instead of FASTQ (BFC `-Q`).
    #[arg(short = 'Q', long = "fasta-out")]
    fasta_out: bool,

    /// libdeflate gzip compression level for `.gz` output (1-12).
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
    fn build_config(&self) -> Result<CorrectConfig> {
        if self.kmer < 11 || self.kmer > 63 || self.kmer.is_multiple_of(2) {
            return Err(RsomicsError::ConfigError(format!(
                "--kmer must be odd and in 11..=63 (BFC contract), got {}",
                self.kmer
            )));
        }
        if self.min_cov < 1 {
            return Err(RsomicsError::ConfigError("--min-cov must be ≥ 1".into()));
        }
        Ok(CorrectConfig {
            k: self.kmer,
            min_cov: self.min_cov,
            win_multi_ec: self.max_ec_window,
            qual_threshold: self.qual_threshold,
            drop_unique_kmer: self.drop_unique_kmer,
            discard_uncorrectable: self.discard_uncorrectable,
            fasta_out: self.fasta_out,
            ..CorrectConfig::default()
        })
    }

    pub fn execute(&self) -> Result<CorrectReport> {
        if self.no_correct && !self.drop_unique_kmer {
            return Err(RsomicsError::ConfigError(
                "--no-correct without --drop-unique-kmer is a pass-through; nothing to do".into(),
            ));
        }
        let cfg = self.build_config()?;
        let p = Pipeline::new(&cfg, self.compression);
        match (self.in2.as_ref(), self.out2.as_ref()) {
            (Some(in2), Some(out2)) => {
                let mut a = p.run(&self.in1, &self.out1)?;
                let b = p.run(in2, out2)?;
                a.reads_in += b.reads_in;
                a.reads_out += b.reads_out;
                a.reads_dropped += b.reads_dropped;
                a.bases_in += b.bases_in;
                a.bases_corrected += b.bases_corrected;
                Ok(a)
            }
            (None, None) => p.run(&self.in1, &self.out1),
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
        upstream: "BFC",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btv592"),
    }),
    usage_lines: &[
        "[OPTIONS] --in1 <PATH> --out1 <PATH>",
        "[OPTIONS] --in1 <R1> --in2 <R2> --out1 <O1> --out2 <O2>   (PE, per-read)",
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
                    description: "R2 input (PE; corrected independently per BFC)",
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
            title: "CORRECTION",
            flags: &[
                FlagSpec {
                    short: Some('k'),
                    long: "kmer",
                    aliases: &[],
                    value: Some("<n>"),
                    type_hint: Some("usize"),
                    required: false,
                    default: Some("33"),
                    description: "K-mer length (odd, 11..=63)",
                    why_default: Some("BFC default -k 33"),
                },
                FlagSpec {
                    short: Some('c'),
                    long: "min-cov",
                    aliases: &[],
                    value: Some("<n>"),
                    type_hint: Some("i32"),
                    required: false,
                    default: Some("3"),
                    description: "Min k-mer coverage for a trusted/solid k-mer",
                    why_default: Some("BFC default -c 3"),
                },
                FlagSpec {
                    short: Some('w'),
                    long: "max-ec-window",
                    aliases: &[],
                    value: Some("<n>"),
                    type_hint: Some("i32"),
                    required: false,
                    default: Some("10"),
                    description: "Suppress a 2nd correction within N bp of the last",
                    why_default: Some("BFC default -w 10"),
                },
                FlagSpec {
                    short: None,
                    long: "qual-threshold",
                    aliases: &[],
                    value: Some("<n>"),
                    type_hint: Some("u8"),
                    required: false,
                    default: Some("20"),
                    description: "Phred ≥ this is high-quality in the penalty model",
                    why_default: Some("BFC default -q 20 (long-only; -q is --quiet)"),
                },
                FlagSpec {
                    short: Some('E'),
                    long: "no-correct",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Skip correction (count pass only; use with -1)",
                    why_default: None,
                },
            ],
        },
        Section {
            title: "OUTPUT POLICY",
            flags: &[
                FlagSpec {
                    short: Some('D'),
                    long: "discard-uncorrectable",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Drop uncorrectable reads instead of passing through",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('1'),
                    long: "drop-unique-kmer",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Drop reads containing a singleton k-mer",
                    why_default: None,
                },
                FlagSpec {
                    short: Some('Q'),
                    long: "fasta-out",
                    aliases: &[],
                    value: None,
                    type_hint: Some("bool"),
                    required: false,
                    default: Some("false"),
                    description: "Emit FASTA (drop qualities)",
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
                    description: "libdeflate gz level 1-12 for .gz output",
                    why_default: Some("ratio/speed trade-off"),
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
            description: "Correct an Illumina FASTQ with BFC defaults (k=33)",
            command: "rsomics-fastq-correct -i reads.fq.gz -o corrected.fq.gz",
        },
        Example {
            description: "Stricter trust threshold, discard uncorrectable, JSON",
            command: "rsomics-fastq-correct -i in.fq -o out.fq -c 5 -D --json | jq .result",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-correct/0.1/#json-output-schema"),
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    // clap debug_assert only fires when the binary parses; without this test a CLI error
    // is invisible to `cargo test` and only surfaces at runtime.
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
