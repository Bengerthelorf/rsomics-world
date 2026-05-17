#![allow(clippy::cast_possible_truncation)]

use std::io;
use std::path::Path;

use noodles::vcf;
use rsomics_common::{Result, RsomicsError};

#[derive(Default)]
pub struct FilterConfig {
    pub min_qual: Option<f32>,
    pub pass_only: bool,
    pub regions: Vec<String>,
    pub include_expr: Option<String>,
    pub exclude_expr: Option<String>,
}

pub struct FilterStats {
    pub total: u64,
    pub passed: u64,
}

pub fn filter_vcf(
    input: &Path,
    output: &mut dyn io::Write,
    cfg: &FilterConfig,
) -> Result<FilterStats> {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;

    let header = reader
        .read_header()
        .map_err(|e| RsomicsError::InvalidInput(format!("reading VCF header: {e}")))?;

    let mut writer = vcf::io::Writer::new(output);
    writer
        .write_header(&header)
        .map_err(|e| RsomicsError::InvalidInput(format!("writing header: {e}")))?;

    let mut stats = FilterStats {
        total: 0,
        passed: 0,
    };

    for result in reader.records() {
        let record =
            result.map_err(|e| RsomicsError::InvalidInput(format!("reading VCF record: {e}")))?;
        stats.total += 1;

        if let Some(min_q) = cfg.min_qual {
            let dominated = record
                .quality_score()
                .map(|r| r.map_or(false, |q| q < min_q))
                .unwrap_or(false);
            if dominated {
                continue;
            }
        }

        writer
            .write_record(&header, &record)
            .map_err(|e| RsomicsError::InvalidInput(format!("writing record: {e}")))?;
        stats.passed += 1;
    }

    Ok(stats)
}
