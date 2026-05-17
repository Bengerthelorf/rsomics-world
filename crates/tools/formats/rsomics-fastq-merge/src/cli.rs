use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use needletail::parse_fastx_file;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_fastq_merge::{analyze, correct, merge};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-fastq-merge", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Read1 FASTQ (gz/bz2/xz/zst auto-detected).
    #[arg(long = "in1", short = 'i')]
    in1: PathBuf,
    /// Read2 FASTQ.
    #[arg(long = "in2", short = 'I')]
    in2: PathBuf,
    /// Merged-reads output (`-` = stdout).
    #[arg(long = "merged_out", short = 'm', default_value = "-")]
    merged_out: String,
    /// Minimum overlap length to call a pair overlapping (fastp default 30).
    #[arg(long = "overlap_len_require", default_value_t = 30)]
    overlap_len_require: usize,
    /// Max mismatches allowed in the overlap (fastp default 5).
    #[arg(long = "overlap_diff_limit", default_value_t = 5)]
    overlap_diff_limit: usize,
    /// Max mismatch percentage allowed in the overlap (fastp default 20).
    #[arg(long = "overlap_diff_percent_limit", default_value_t = 20)]
    overlap_diff_percent_limit: u32,
    /// Correct low-quality overlap mismatches from the high-quality mate
    /// before merging (fastp `--correction`).
    #[arg(long = "correction")]
    correction: bool,
    /// Also emit reads whose pair did not overlap (fastp
    /// `--include_unmerged`): read1 then read2, both as-is.
    #[arg(long = "include_unmerged")]
    include_unmerged: bool,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(&self) -> Result<()> {
        let mut r1 = parse_fastx_file(&self.in1).map_err(|e| {
            RsomicsError::InvalidInput(format!("opening {}: {e}", self.in1.display()))
        })?;
        let mut r2 = parse_fastx_file(&self.in2).map_err(|e| {
            RsomicsError::InvalidInput(format!("opening {}: {e}", self.in2.display()))
        })?;

        let mut out: Box<dyn Write> = if self.merged_out == "-" {
            Box::new(BufWriter::new(io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(&self.merged_out).map_err(RsomicsError::Io)?,
            ))
        };
        let diff_pct = f64::from(self.overlap_diff_percent_limit) / 100.0;
        let mut pair_no: u64 = 0;

        loop {
            let (next1, next2) = (r1.next(), r2.next());
            match (next1, next2) {
                (Some(rec1), Some(rec2)) => {
                    pair_no += 1;
                    let ra = rec1.map_err(|e| {
                        RsomicsError::InvalidInput(format!("{}: {e}", self.in1.display()))
                    })?;
                    let rb = rec2.map_err(|e| {
                        RsomicsError::InvalidInput(format!("{}: {e}", self.in2.display()))
                    })?;
                    let id1 = ra.id().to_vec();
                    let id2 = rb.id().to_vec();
                    let mut s1 = ra.seq().to_vec();
                    let mut q1 = ra
                        .qual()
                        .ok_or_else(|| {
                            RsomicsError::InvalidInput(format!(
                                "{}: record {pair_no} has no quality — not FASTQ",
                                self.in1.display()
                            ))
                        })?
                        .to_vec();
                    let mut s2 = rb.seq().to_vec();
                    let mut q2 = rb
                        .qual()
                        .ok_or_else(|| {
                            RsomicsError::InvalidInput(format!(
                                "{}: record {pair_no} has no quality — not FASTQ",
                                self.in2.display()
                            ))
                        })?
                        .to_vec();

                    let ov = analyze(
                        &s1,
                        &s2,
                        self.overlap_diff_limit,
                        self.overlap_len_require,
                        diff_pct,
                    );
                    if self.correction {
                        correct(&mut s1, &mut q1, &mut s2, &mut q2, &ov);
                    }
                    if ov.overlapped {
                        let merged = merge(&s1, &q1, &s2, &q2, &ov)
                            .expect("ov.overlapped ⇒ merge yields a read");
                        write_record(
                            &mut out,
                            &format!(
                                "{} merged_{}_{}",
                                String::from_utf8_lossy(&id1),
                                merged.len1,
                                merged.len2
                            ),
                            &merged.seq,
                            &merged.qual,
                        )?;
                    } else if self.include_unmerged {
                        // fastp writes both mates to the merged stream (original orientation) so the file stays a self-consistent read set.
                        write_record(&mut out, &String::from_utf8_lossy(&id1), &s1, &q1)?;
                        write_record(&mut out, &String::from_utf8_lossy(&id2), &s2, &q2)?;
                    }
                }
                (None, None) => break,
                _ => {
                    return Err(RsomicsError::InvalidInput(
                        "in1 and in2 have different read counts (not properly paired)".into(),
                    ));
                }
            }
        }
        out.flush().map_err(RsomicsError::Io)
    }
}

fn write_record(w: &mut dyn Write, name: &str, seq: &[u8], qual: &[u8]) -> Result<()> {
    w.write_all(b"@").map_err(RsomicsError::Io)?;
    w.write_all(name.as_bytes()).map_err(RsomicsError::Io)?;
    w.write_all(b"\n").map_err(RsomicsError::Io)?;
    w.write_all(seq).map_err(RsomicsError::Io)?;
    w.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
    w.write_all(qual).map_err(RsomicsError::Io)?;
    w.write_all(b"\n").map_err(RsomicsError::Io)
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Merge overlapping paired-end reads into consensus reads (Rust port of fastp --merge).",
    origin: Some(Origin {
        upstream: "fastp",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/bty560"),
    }),
    usage_lines: &["--in1 R1.fq --in2 R2.fq [--correction] [-m merged.fq]"],
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
                description: "Read1 FASTQ",
                why_default: None,
            },
            FlagSpec {
                short: Some('I'),
                long: "in2",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Read2 FASTQ",
                why_default: None,
            },
            FlagSpec {
                short: Some('m'),
                long: "merged_out",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("-"),
                description: "Merged output (default stdout)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "overlap_len_require",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("30"),
                description: "Min overlap length",
                why_default: Some("fastp default"),
            },
            FlagSpec {
                short: None,
                long: "overlap_diff_limit",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("5"),
                description: "Max mismatches in overlap",
                why_default: Some("fastp default"),
            },
            FlagSpec {
                short: None,
                long: "overlap_diff_percent_limit",
                aliases: &[],
                value: Some("<PCT>"),
                type_hint: Some("u32"),
                required: false,
                default: Some("20"),
                description: "Max mismatch % in overlap",
                why_default: Some("fastp default"),
            },
            FlagSpec {
                short: None,
                long: "correction",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Quality-correct overlap mismatches before merge",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "include_unmerged",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Also emit non-overlapping reads",
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
            description: "Merge with quality correction",
            command: "rsomics-fastq-merge --in1 R1.fq.gz --in2 R2.fq.gz --correction -m merged.fq",
        },
        Example {
            description: "Stream merged reads to a pipe",
            command: "rsomics-fastq-merge -i R1.fq -I R2.fq | gzip > merged.fq.gz",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    // clap debug_assert fires only when a binary parses; without this test, CLI-definition errors (duplicate shorts, id clashes) are invisible to `cargo test`.
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
