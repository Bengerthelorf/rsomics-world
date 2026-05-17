use std::path::PathBuf;

use clap::Parser;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_fqgz::ChunkedWriter;
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_seqio::{OwnedRecord, open_fastq};

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

const CHUNK_RECORDS: usize = 8192;

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

        let threads = self.common.threads.unwrap_or(0);
        if threads > 0 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .ok();
        }

        let serial = threads == 1;
        if let Some(in2) = &self.in2 {
            self.run_pe(&refs, &cfg, in2)?;
        } else if serial {
            self.run_se_serial(&refs, &cfg)?;
        } else {
            self.run_se(&refs, &cfg)?;
        }
        Ok(())
    }

    fn run_se_serial(&self, refs: &RefKmers, cfg: &Config) -> Result<()> {
        use std::io::{BufWriter, Write};

        let mut reader = open_fastq(&self.in1)?;

        let stdout_mode = self.out1 == "-";
        let mut stdout_out =
            stdout_mode.then(|| BufWriter::with_capacity(256 * 1024, std::io::stdout().lock()));
        let mut file_out = if stdout_mode {
            None
        } else {
            Some(BufWriter::with_capacity(
                256 * 1024,
                std::fs::File::create(&self.out1).map_err(RsomicsError::Io)?,
            ))
        };

        for result in reader.by_ref() {
            let rec = result?;
            let (id, seq, qual) = (&rec.id, &rec.seq, &rec.qual);

            if let Some((s, e)) = process(&seq, &qual, refs, cfg) {
                if let Some(out) = stdout_out.as_mut() {
                    out.write_all(b"@").map_err(RsomicsError::Io)?;
                    out.write_all(id).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n").map_err(RsomicsError::Io)?;
                    out.write_all(&seq[s..e]).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
                    out.write_all(&qual[s..e]).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n").map_err(RsomicsError::Io)?;
                } else if let Some(out) = file_out.as_mut() {
                    out.write_all(b"@").map_err(RsomicsError::Io)?;
                    out.write_all(id).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n").map_err(RsomicsError::Io)?;
                    out.write_all(&seq[s..e]).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
                    out.write_all(&qual[s..e]).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n").map_err(RsomicsError::Io)?;
                }
            }
        }
        if let Some(out) = stdout_out.as_mut() {
            out.flush().map_err(RsomicsError::Io)?;
        }
        if let Some(out) = file_out.as_mut() {
            out.flush().map_err(RsomicsError::Io)?;
        }
        Ok(())
    }

    fn run_se(&self, refs: &RefKmers, cfg: &Config) -> Result<()> {
        use std::io::{BufWriter, Write};

        let mut reader = open_fastq(&self.in1)?;
        let has_outm = self.outm.is_some();
        let mut chunk: Vec<OwnedRecord> = Vec::with_capacity(CHUNK_RECORDS);
        let mut done = false;

        if self.out1 == "-" {
            let mut out = BufWriter::new(std::io::stdout().lock());
            while !done {
                chunk.clear();
                while chunk.len() < CHUNK_RECORDS {
                    let Some(rec) = reader.next() else {
                        done = true;
                        break;
                    };
                    chunk.push(rec?);
                }
                if chunk.is_empty() {
                    break;
                }
                let results: Vec<Option<OwnedRecord>> = chunk
                    .par_drain(..)
                    .map(|rec| match process(&rec.seq, &rec.qual, refs, cfg) {
                        Some((s, e)) => Some(OwnedRecord {
                            id: rec.id,
                            seq: rec.seq[s..e].to_vec(),
                            qual: rec.qual[s..e].to_vec(),
                        }),
                        None => None,
                    })
                    .collect();
                for rec in results.into_iter().flatten() {
                    out.write_all(b"@").map_err(RsomicsError::Io)?;
                    out.write_all(&rec.id).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n").map_err(RsomicsError::Io)?;
                    out.write_all(&rec.seq).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n+\n").map_err(RsomicsError::Io)?;
                    out.write_all(&rec.qual).map_err(RsomicsError::Io)?;
                    out.write_all(b"\n").map_err(RsomicsError::Io)?;
                }
            }
            out.flush().map_err(RsomicsError::Io)?;
        } else {
            let mut w1 = ChunkedWriter::create(std::path::Path::new(&self.out1), 4)?;
            let mut wm = self
                .outm
                .as_ref()
                .map(|p| ChunkedWriter::create(std::path::Path::new(p), 4))
                .transpose()?;
            while !done {
                chunk.clear();
                while chunk.len() < CHUNK_RECORDS {
                    let Some(rec) = reader.next() else {
                        done = true;
                        break;
                    };
                    chunk.push(rec?);
                }
                if chunk.is_empty() {
                    break;
                }
                let results: Vec<(Option<OwnedRecord>, Option<OwnedRecord>)> = chunk
                    .par_drain(..)
                    .map(|rec| match process(&rec.seq, &rec.qual, refs, cfg) {
                        Some((s, e)) => (
                            Some(OwnedRecord {
                                id: rec.id,
                                seq: rec.seq[s..e].to_vec(),
                                qual: rec.qual[s..e].to_vec(),
                            }),
                            None,
                        ),
                        None => (None, if has_outm { Some(rec) } else { None }),
                    })
                    .collect();
                for (pass, fail) in results {
                    if let Some(rec) = pass {
                        w1.write_record(&rec.id, &rec.seq, &rec.qual)?;
                    } else if let (Some(rec), Some(wm)) = (fail, wm.as_mut()) {
                        wm.write_record(&rec.id, &rec.seq, &rec.qual)?;
                    }
                }
            }
            w1.finalize()?;
            if let Some(wm) = wm {
                wm.finalize()?;
            }
        }
        Ok(())
    }

    fn run_pe(&self, refs: &RefKmers, cfg: &Config, in2: &std::path::Path) -> Result<()> {
        let mut r1 = open_fastq(&self.in1)?;
        let mut r2 = open_fastq(in2)?;
        let out2 = self.out2.as_ref().expect("paired ⇒ out2 set");
        let mut w1 = ChunkedWriter::create(std::path::Path::new(&self.out1), 4)?;
        let mut w2 = ChunkedWriter::create(std::path::Path::new(out2), 4)?;

        let mut chunk: Vec<(OwnedRecord, OwnedRecord)> = Vec::with_capacity(CHUNK_RECORDS);
        let mut done = false;

        while !done {
            chunk.clear();
            while chunk.len() < CHUNK_RECORDS {
                match (r1.next(), r2.next()) {
                    (Some(a), Some(b)) => chunk.push((a?, b?)),
                    (None, None) => {
                        done = true;
                        break;
                    }
                    _ => {
                        return Err(RsomicsError::InvalidInput(
                            "in1 and in2 have different read counts (not properly paired)".into(),
                        ));
                    }
                }
            }
            if chunk.is_empty() {
                break;
            }

            // BBDuk removeifeitherbad=t: pair kept/dropped as unit
            let results: Vec<Option<(OwnedRecord, OwnedRecord)>> = chunk
                .par_drain(..)
                .map(|(a, b)| {
                    match (
                        process(&a.seq, &a.qual, refs, cfg),
                        process(&b.seq, &b.qual, refs, cfg),
                    ) {
                        (Some(t1), Some(t2)) => Some((
                            OwnedRecord {
                                id: a.id,
                                seq: a.seq[t1.0..t1.1].to_vec(),
                                qual: a.qual[t1.0..t1.1].to_vec(),
                            },
                            OwnedRecord {
                                id: b.id,
                                seq: b.seq[t2.0..t2.1].to_vec(),
                                qual: b.qual[t2.0..t2.1].to_vec(),
                            },
                        )),
                        _ => None,
                    }
                })
                .collect();

            for (a, b) in results.into_iter().flatten() {
                w1.write_record(&a.id, &a.seq, &a.qual)?;
                w2.write_record(&b.id, &b.seq, &b.qual)?;
            }
        }
        w1.finalize()?;
        w2.finalize()?;
        Ok(())
    }
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
