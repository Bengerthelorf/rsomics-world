use std::fs::File;
use std::io::{self, BufWriter, Read};
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};
use rsomics_intervals::{IntervalSet, bed};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-sort", disable_help_flag = true)]
pub struct Cli {
    #[arg(short = 'i', long = "input", default_value = "-")]
    input: String,
    #[arg(short = 'o', long = "output", default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(&self) -> Result<()> {
        let intervals = if self.input == "-" {
            let mut buf = Vec::new();
            io::stdin()
                .read_to_end(&mut buf)
                .map_err(RsomicsError::Io)?;
            bed::read_bytes(&buf)?
        } else {
            bed::read(File::open(PathBuf::from(&self.input)).map_err(RsomicsError::Io)?)?
        };
        let mut set: IntervalSet = intervals.into_iter().collect();
        set.sort();
        let writer: Box<dyn io::Write> = if self.output == "-" {
            Box::new(BufWriter::new(io::stdout().lock()))
        } else {
            Box::new(BufWriter::new(
                File::create(PathBuf::from(&self.output)).map_err(RsomicsError::Io)?,
            ))
        };
        bed::write_bed3(writer, set.iter().cloned())
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Sort BED records by chromosome (lexicographic) then start position (bedtools sort equivalent).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -i <BED>", "[OPTIONS] < input.bed"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('i'),
                long: "input",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("-"),
                description: "Input BED (default stdin)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("-"),
                description: "Output BED (default stdout)",
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
    examples: &[Example {
        description: "Sort a BED file",
        command: "rsomics-bed-sort -i unsorted.bed -o sorted.bed",
    }],
    json_result_schema_doc: None,
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
