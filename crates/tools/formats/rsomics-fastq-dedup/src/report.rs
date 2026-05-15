use serde::Serialize;

#[derive(Debug, Default, Clone, Serialize)]
pub struct DedupReport {
    pub mode: Option<&'static str>,
    pub input: Option<String>,
    pub output: Option<String>,
    pub reads_in: u64,
    pub reads_out: u64,
    pub duplicates_removed: u64,
}
