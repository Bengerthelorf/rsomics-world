use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;

use clap::Parser;
use needletail::parse_fastx_file;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bbduk::{Config, KTrim, MAX_HDIST, MAX_K, QTrim, RefKmers, process};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bbduk", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// Read1 FASTQ (gz/bz2/xz/zst auto-detected).
    #[arg(long = "in1", short = 'i')]
    in1: PathBuf,
    /// Read2 FASTQ (paired-end). Reads are processed as pairs.
    #[arg(long = "in2", short = 'I')]
    in2: Option<PathBuf>,
    /// Read1 output (`-` = stdout).
    #[arg(long = "out1", short = 'o', default_value = "-")]
    out1: String,
    /// Read2 output (required iff `--in2`).
    #[arg(long = "out2", short = 'O')]
    out2: Option<String>,
    /// Matched (removed-as-contaminant) reads output. kfilter mode only.
    #[arg(long = "outm")]
    outm: Option<String>,
    /// Reference FASTA(s) of contaminant/adapter sequences (BBDuk `ref=`).
    #[arg(long = "ref", short = 'r')]
    refs: Vec<PathBuf>,
    /// Comma-separated literal reference sequences (BBDuk `literal=`).
    #[arg(long = "literal")]
    literal: Option<String>,
    /// K-mer length (BBDuk default 27; max 31).
    #[arg(long = "k", default_value_t = 27)]
    k: usize,
    /// Use shorter ref k-mers down to this length at read tips for partial
    /// adapters (BBDuk `mink`; 0 = disabled).
    #[arg(long = "mink", default_value_t = 0)]
    mink: usize,
    /// Index k-mers within this Hamming distance of each ref k-mer
    /// (BBDuk `hdist`; 0 = exact; max 3).
    #[arg(long = "hdist", default_value_t = 0)]
    hdist: usize,
    /// Hamming distance for the short (`--mink`) tip k-mers (BBDuk
    /// `hdist2`; default 0 = exact tips, independent of `--hdist`; max 3).
    #[arg(long = "hdist2", default_value_t = 0)]
    hdist2: usize,
    /// Do NOT also match reverse complements (BBDuk `rcomp=f`; default on).
    #[arg(long = "no-rcomp")]
    no_rcomp: bool,
    /// Do NOT wildcard the middle k-mer base (BBDuk `maskmiddle=f`;
    /// default on, and BBDuk-style auto-disabled when `--mink>0`).
    #[arg(long = "no-maskmiddle")]
    no_maskmiddle: bool,
    /// A read is a contaminant at ≥ this many shared k-mers (BBDuk `mkh`).
    #[arg(long = "minkmerhits", visible_alias = "mkh", default_value_t = 1)]
    min_kmer_hits: usize,
    /// …or at ≥ this fraction of its k-mers (BBDuk `mkf`; 0 = unused).
    #[arg(long = "minkmerfraction", visible_alias = "mkf", default_value_t = 0.0)]
    min_kmer_fraction: f64,
    /// K-mer trim mode: `f` (filter, default) | `r` (3') | `l` (5').
    #[arg(long = "ktrim", default_value = "f")]
    ktrim: String,
    /// Quality trim mode: `f` (default) | `r` | `l` | `rl`.
    #[arg(long = "qtrim", default_value = "f")]
    qtrim: String,
    /// Phred quality-trim threshold (BBDuk `trimq`).
    #[arg(long = "trimq", default_value_t = 6)]
    trimq: u8,
    /// Discard reads shorter than this after trimming (BBDuk `minlength`).
    #[arg(long = "minlength", visible_alias = "minlen", default_value_t = 10)]
    min_length: usize,

    #[command(flatten)]
    pub common: CommonFlags,
}

fn parse_ktrim(s: &str) -> Result<KTrim> {
    match s {
        "f" | "false" => Ok(KTrim::None),
        "r" => Ok(KTrim::Right),
        "l" => Ok(KTrim::Left),
        _ => Err(RsomicsError::InvalidInput(format!(
            "--ktrim must be f|r|l (got {s:?})"
        ))),
    }
}

fn parse_qtrim(s: &str) -> Result<QTrim> {
    match s {
        "f" | "false" => Ok(QTrim::None),
        "r" => Ok(QTrim::Right),
        "l" => Ok(QTrim::Left),
        "rl" | "lr" => Ok(QTrim::Both),
        _ => Err(RsomicsError::InvalidInput(format!(
            "--qtrim must be f|r|l|rl (got {s:?})"
        ))),
    }
}

fn writer(path: &str) -> Result<Box<dyn Write>> {
    if path == "-" {
        Ok(Box::new(BufWriter::new(io::stdout().lock())))
    } else {
        Ok(Box::new(BufWriter::new(
            File::create(path).map_err(RsomicsError::Io)?,
        )))
    }
}

impl Cli {
    fn config(&self) -> Result<Config> {
        if !(1..=MAX_K).contains(&self.k) {
            return Err(RsomicsError::InvalidInput(format!(
                "--k must be in 1..={MAX_K} (2-bit codec / rcomp limit); got {}",
                self.k
            )));
        }
        if self.hdist > MAX_HDIST || self.hdist2 > MAX_HDIST {
            return Err(RsomicsError::InvalidInput(format!(
                "--hdist/--hdist2 > {MAX_HDIST} expands the k-mer index past memory; \
                 got hdist={} hdist2={}",
                self.hdist, self.hdist2
            )));
        }
        if self.mink > self.k {
            return Err(RsomicsError::InvalidInput(format!(
                "--mink ({}) must be ≤ --k ({})",
                self.mink, self.k
            )));
        }
        Ok(Config {
            k: self.k,
            mink: self.mink,
            hdist: self.hdist,
            hdist2: self.hdist2,
            rcomp: !self.no_rcomp,
            maskmiddle: !self.no_maskmiddle,
            min_kmer_hits: self.min_kmer_hits,
            min_kmer_fraction: self.min_kmer_fraction,
            ktrim: parse_ktrim(&self.ktrim)?,
            qtrim: parse_qtrim(&self.qtrim)?,
            trimq: self.trimq,
            min_length: self.min_length,
            qual_offset: 33,
        })
    }

    fn load_refs(&self) -> Result<Vec<Vec<u8>>> {
        let mut out = Vec::new();
        for p in &self.refs {
            let mut rdr = parse_fastx_file(p).map_err(|e| {
                RsomicsError::InvalidInput(format!("opening ref {}: {e}", p.display()))
            })?;
            while let Some(rec) = rdr.next() {
                let rec = rec
                    .map_err(|e| RsomicsError::InvalidInput(format!("ref {}: {e}", p.display())))?;
                out.push(rec.seq().into_owned());
            }
        }
        if let Some(lit) = &self.literal {
            for s in lit.split(',').filter(|s| !s.is_empty()) {
                out.push(s.as_bytes().to_vec());
            }
        }
        Ok(out)
    }

    pub fn execute(&self) -> Result<()> {
        let cfg = self.config()?;
        let ref_seqs = self.load_refs()?;
        // --ktrim r|l with no reference is a user mistake — fail loud; kfilter with no ref is a BBDuk-faithful no-op (quality-trim-only) and is allowed
        if cfg.ktrim != KTrim::None && ref_seqs.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "--ktrim r|l needs a reference: pass --ref FILE and/or --literal SEQ".into(),
            ));
        }
        let refs = RefKmers::build(ref_seqs.iter().map(Vec::as_slice), &cfg);

        if self.in2.is_some() != self.out2.is_some() {
            return Err(RsomicsError::InvalidInput(
                "--in2 and --out2 must be given together (paired-end)".into(),
            ));
        }

        let mut r1 = parse_fastx_file(&self.in1).map_err(|e| {
            RsomicsError::InvalidInput(format!("opening {}: {e}", self.in1.display()))
        })?;
        let mut w1 = writer(&self.out1)?;
        let mut wm = self.outm.as_deref().map(writer).transpose()?;

        if let Some(in2) = &self.in2 {
            let mut r2 = parse_fastx_file(in2).map_err(|e| {
                RsomicsError::InvalidInput(format!("opening {}: {e}", in2.display()))
            })?;
            let mut w2 = writer(self.out2.as_deref().expect("paired ⇒ out2 set"))?;
            loop {
                match (r1.next(), r2.next()) {
                    (Some(a), Some(b)) => {
                        let a = a.map_err(|e| {
                            RsomicsError::InvalidInput(format!("{}: {e}", self.in1.display()))
                        })?;
                        let b = b.map_err(|e| {
                            RsomicsError::InvalidInput(format!("{}: {e}", in2.display()))
                        })?;
                        let q1 = a.qual().ok_or_else(|| no_qual(&self.in1))?;
                        let q2 = b.qual().ok_or_else(|| no_qual(in2))?;
                        let (s1, s2) = (a.raw_seq(), b.raw_seq());
                        // BBDuk removeifeitherbad=t: a pair is kept/dropped as a unit.
                        match (process(s1, q1, &refs, &cfg), process(s2, q2, &refs, &cfg)) {
                            (Some(t1), Some(t2)) => {
                                write_fq(&mut w1, a.id(), &s1[t1.0..t1.1], &q1[t1.0..t1.1])?;
                                write_fq(&mut w2, b.id(), &s2[t2.0..t2.1], &q2[t2.0..t2.1])?;
                            }
                            _ => {
                                if let Some(wm) = wm.as_mut() {
                                    write_fq(wm, a.id(), s1, q1)?;
                                    write_fq(wm, b.id(), s2, q2)?;
                                }
                            }
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
            w2.flush().map_err(RsomicsError::Io)?;
        } else {
            while let Some(rec) = r1.next() {
                let rec = rec.map_err(|e| {
                    RsomicsError::InvalidInput(format!("{}: {e}", self.in1.display()))
                })?;
                let q = rec.qual().ok_or_else(|| no_qual(&self.in1))?;
                let s = rec.raw_seq();
                if let Some((a, e)) = process(s, q, &refs, &cfg) {
                    write_fq(&mut w1, rec.id(), &s[a..e], &q[a..e])?;
                } else if let Some(wm) = wm.as_mut() {
                    write_fq(wm, rec.id(), s, q)?;
                }
            }
        }

        if let Some(mut wm) = wm {
            wm.flush().map_err(RsomicsError::Io)?;
        }
        w1.flush().map_err(RsomicsError::Io)
    }
}

fn no_qual(path: &std::path::Path) -> RsomicsError {
    RsomicsError::InvalidInput(format!(
        "{}: a record has no quality — not FASTQ",
        path.display()
    ))
}

fn write_fq(w: &mut dyn Write, id: &[u8], seq: &[u8], qual: &[u8]) -> Result<()> {
    w.write_all(b"@").map_err(RsomicsError::Io)?;
    w.write_all(id).map_err(RsomicsError::Io)?;
    w.write_all(b"\n").map_err(RsomicsError::Io)?;
    w.write_all(seq).map_err(RsomicsError::Io)?;
    w.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
    w.write_all(qual).map_err(RsomicsError::Io)?;
    w.write_all(b"\n").map_err(RsomicsError::Io)
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "K-mer contaminant removal + adapter/quality trimming (clean-room Rust BBDuk).",
    origin: Some(Origin {
        upstream: "BBDuk (BBTools)",
        upstream_license: "BBTools license (free, redistribution-restricted)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: None,
    }),
    usage_lines: &[
        "--in1 R.fq --ref adapters.fa --ktrim r --out1 clean.fq",
        "-i R1.fq -I R2.fq -r contam.fa -o o1.fq -O o2.fq",
    ],
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
                required: false,
                default: None,
                description: "Read2 FASTQ (paired-end)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "out1",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("-"),
                description: "Read1 output (default stdout)",
                why_default: None,
            },
            FlagSpec {
                short: Some('O'),
                long: "out2",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Read2 output (required iff --in2)",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "outm",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Matched/removed reads output",
                why_default: None,
            },
            FlagSpec {
                short: Some('r'),
                long: "ref",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: None,
                description: "Reference FASTA of contaminants/adapters",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "literal",
                aliases: &[],
                value: Some("<seqs>"),
                type_hint: Some("String"),
                required: false,
                default: None,
                description: "Comma-separated literal reference seqs",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "k",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("27"),
                description: "K-mer length (max 31)",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "mink",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("0"),
                description: "Shorter k-mers at read tips (0=off)",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "hdist",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("0"),
                description: "Hamming distance for k-mer match (max 3)",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "hdist2",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("0"),
                description: "Hamming distance for short tip k-mers (max 3)",
                why_default: Some("BBDuk default (exact tips)"),
            },
            FlagSpec {
                short: None,
                long: "no-rcomp",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Disable reverse-complement matching",
                why_default: Some("BBDuk rcomp=t"),
            },
            FlagSpec {
                short: None,
                long: "no-maskmiddle",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: Some("false"),
                description: "Disable middle-base wildcard",
                why_default: Some("BBDuk maskmiddle=t"),
            },
            FlagSpec {
                short: None,
                long: "minkmerhits",
                aliases: &["mkh"],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("1"),
                description: "Min shared k-mers to call contaminant",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "minkmerfraction",
                aliases: &["mkf"],
                value: Some("<F>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.0"),
                description: "…or min fraction of k-mers (0=unused)",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "ktrim",
                aliases: &[],
                value: Some("<f|r|l>"),
                type_hint: Some("String"),
                required: false,
                default: Some("f"),
                description: "K-mer trim mode (f=filter)",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "qtrim",
                aliases: &[],
                value: Some("<f|r|l|rl>"),
                type_hint: Some("String"),
                required: false,
                default: Some("f"),
                description: "Quality trim mode",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "trimq",
                aliases: &[],
                value: Some("<Q>"),
                type_hint: Some("u8"),
                required: false,
                default: Some("6"),
                description: "Phred quality-trim threshold",
                why_default: Some("BBDuk default"),
            },
            FlagSpec {
                short: None,
                long: "minlength",
                aliases: &["minlen"],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("10"),
                description: "Discard reads shorter than this post-trim",
                why_default: Some("BBDuk default"),
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
            description: "Right-trim Illumina adapters (k-mer + partial-tip)",
            command: "rsomics-bbduk -i R.fq.gz -r adapters.fa --ktrim r --mink 11 --hdist 1 -o clean.fq",
        },
        Example {
            description: "Remove PhiX-contaminated read pairs",
            command: "rsomics-bbduk -i R1.fq -I R2.fq -r phix.fa -o o1.fq -O o2.fq --outm phix.fq",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    // clap debug_assert validates the whole arg graph (shorts, flattened CommonFlags, alias clashes) — fires only on binary parse, so a lib-only test suite misses a CLI-definition error
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
