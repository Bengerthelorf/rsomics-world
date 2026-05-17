mod adapter_content;
mod basic_stats;
mod dup_levels;
mod kmer_content;
mod n_content;
mod overrepresented;
mod per_base_content;
mod per_base_quality;
mod per_seq_gc;
mod per_seq_quality;
mod per_tile_quality;
mod seq_length;

pub use adapter_content::AdapterContent;
pub use basic_stats::BasicStats;
pub use dup_levels::DuplicationLevels;
pub use kmer_content::KmerContent;
pub use n_content::PerBaseNContent;
pub use overrepresented::OverrepresentedSeqs;
pub use per_base_content::PerBaseContent;
pub use per_base_quality::PerBaseQuality;
pub use per_seq_gc::PerSeqGc;
pub use per_seq_quality::PerSeqQuality;
pub use per_tile_quality::PerTileQuality;
pub use seq_length::SeqLengthDistribution;

pub struct Record<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ModuleStatus {
    Pass,
    Warn,
    Fail,
}

impl ModuleStatus {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Pass => "PASS",
            Self::Warn => "WARN",
            Self::Fail => "FAIL",
        }
    }

    #[must_use]
    pub const fn as_data_token(self) -> &'static str {
        match self {
            Self::Pass => "pass",
            Self::Warn => "warn",
            Self::Fail => "fail",
        }
    }
}

pub trait QcModule {
    // verbatim FastQC name — used as the `>>` header in fastqc_data.txt and the column in summary.txt
    fn name(&self) -> &'static str;

    fn process(&mut self, rec: &Record);

    fn finalize(&mut self);

    fn status(&self) -> ModuleStatus;

    // appends the fastqc_data.txt body between `>>name` and `>>END_MODULE`
    fn write_data(&self, out: &mut String);
}

// modules in FastQC's report order; `filename` becomes the Basic Statistics `Filename` field
#[must_use]
pub fn default_modules(filename: &str) -> Vec<Box<dyn QcModule>> {
    vec![
        Box::new(BasicStats::new(filename)),
        Box::new(PerBaseQuality::new()),
        Box::new(PerTileQuality::new()),
        Box::new(PerSeqQuality::new()),
        Box::new(PerBaseContent::new()),
        Box::new(PerSeqGc::new()),
        Box::new(PerBaseNContent::new()),
        Box::new(SeqLengthDistribution::new()),
        Box::new(DuplicationLevels::new()),
        Box::new(OverrepresentedSeqs::new()),
        Box::new(AdapterContent::new()),
        Box::new(KmerContent::new()),
    ]
}
