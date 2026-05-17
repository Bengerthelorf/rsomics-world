use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Context, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fastq_stats::{Config, FastqStats, compute_stats, render_pretty, render_tabular};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-stats", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// FASTQ file(s). Gzip / bzip2 / xz / zstd inputs are auto-detected
    /// by extension or magic bytes. Stdin (`-`) is not yet supported.
    #[arg(required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Emit extended statistics: `Q1` / `Q2` / `Q3`, `sum_gap`, `N50`,
    /// `N50_num`, `Q20(%)`, `Q30(%)`, `AvgQual`, `GC(%)`, `sum_n`.
    #[arg(short = 'a', long = "all")]
    all: bool,

    /// Tab-separated machine-readable output. One header line, one row per
    /// input. Disables thousand-separator commas.
    #[arg(short = 'T', long = "tabular")]
    tabular: bool,

    /// Characters counted as gap when `--all` is set. Default matches
    /// seqkit: hyphen, space, period.
    #[arg(short = 'G', long = "gap-letters", default_value = "- .")]
    gap_letters: String,

    /// FASTQ quality encoding: `sanger` / `illumina-1.8+` (offset 33) or
    /// `solexa` / `illumina-1.3+` / `illumina-1.5+` (offset 64).
    #[arg(short = 'E', long = "fq-encoding", default_value = "sanger")]
    fq_encoding: String,

    /// Show only each file's basename in the `file` column (seqkit -b).
    #[arg(short = 'b', long = "basename")]
    basename: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

fn encoding_offset(name: &str) -> Result<u8> {
    match name {
        "sanger" | "illumina-1.8+" => Ok(33),
        "solexa" | "illumina-1.3+" | "illumina-1.5+" => Ok(64),
        other => Err(RsomicsError::InvalidInput(format!(
            "unknown --fq-encoding {other:?} (expected sanger / illumina-1.8+ / solexa / illumina-1.3+ / illumina-1.5+)"
        ))),
    }
}

impl Cli {
    pub fn execute(&self) -> Result<Vec<FastqStats>> {
        let cfg = Config {
            extended: self.all,
            gap_letters: self.gap_letters.as_bytes().to_vec(),
            qual_offset: encoding_offset(&self.fq_encoding)?,
            basename: self.basename,
        };
        let mut results = Vec::with_capacity(self.inputs.len());
        for input in &self.inputs {
            if input.as_os_str() == "-" {
                return Err(RsomicsError::InvalidInput(
                    "stdin (`-`) input not yet supported; pass a file path".into(),
                ));
            }
            let stats = compute_stats(input, &cfg)
                .rs_with_context(|| format!("computing stats for {}", input.display()))?;
            results.push(stats);
        }

        if !self.common.json {
            emit_stdout(&results, self.tabular);
        }
        Ok(results)
    }
}

fn emit_stdout(results: &[FastqStats], tabular: bool) {
    if tabular {
        for (i, s) in results.iter().enumerate() {
            let rendered = render_tabular(s);
            if i == 0 {
                print!("{rendered}");
            } else {
                for line in rendered.lines().skip(1) {
                    println!("{line}");
                }
            }
        }
    } else {
        for s in results {
            print!("{}", render_pretty(s));
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Per-file statistics for FASTQ inputs (Rust port of seqkit stats — FASTQ).",
    origin: Some(Origin {
        upstream: "seqkit",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1371/journal.pone.0163962"),
    }),
    usage_lines: &["[OPTIONS] <INPUTS>..."],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('a'),
                long: "all",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit extended stats (Q1/Q2/Q3, N50, Q20/Q30%, AvgQual, GC%, gaps)",
                why_default: Some(
                    "seqkit's `--all` opts in to the extended columns; off by default",
                ),
            },
            FlagSpec {
                short: Some('T'),
                long: "tabular",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Tab-separated machine-readable output",
                why_default: None,
            },
            FlagSpec {
                short: Some('G'),
                long: "gap-letters",
                aliases: &[],
                value: Some("<CHARS>"),
                type_hint: Some("String"),
                required: false,
                default: Some("- ."),
                description: "Characters counted as gap when --all is set",
                why_default: Some("matches seqkit's default gap letter set"),
            },
            FlagSpec {
                short: Some('E'),
                long: "fq-encoding",
                aliases: &[],
                value: Some("<NAME>"),
                type_hint: Some("String"),
                required: false,
                default: Some("sanger"),
                description: "Quality encoding (sanger/illumina-1.8+ = 33, solexa/1.3+/1.5+ = 64)",
                why_default: Some("sanger (offset 33) is the universal modern Illumina encoding"),
            },
            FlagSpec {
                short: Some('b'),
                long: "basename",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Show only the file's basename in the `file` column",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "json",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Emit AI-friendly JSON envelope on stdout",
                why_default: None,
            },
            FlagSpec {
                short: Some('t'),
                long: "threads",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: None,
                description: "Worker thread count (default: available cores)",
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
            description: "Default stats",
            command: "rsomics-fastq-stats reads.fq",
        },
        Example {
            description: "Extended, tabular, on gzip input",
            command: "rsomics-fastq-stats --tabular --all reads.fastq.gz",
        },
        Example {
            description: "JSON envelope through jq",
            command: "rsomics-fastq-stats --json reads.fq | jq .result",
        },
    ],
    json_result_schema_doc: Some("https://docs.rs/rsomics-fastq-stats/0.1/#json-output-schema"),
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    // debug_assert validates the full arg graph including flattened CommonFlags; only fires at binary parse time.
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
