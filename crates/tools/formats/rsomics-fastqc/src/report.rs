use std::fmt::Write as _;

use crate::modules::QcModule;

#[must_use]
pub fn fastqc_data(modules: &[Box<dyn QcModule>]) -> String {
    let mut out = String::with_capacity(4096);
    out.push_str("##FastQC\t0.12.1\n");
    for m in modules {
        let _ = writeln!(out, ">>{}\t{}", m.name(), m.status().as_data_token());
        m.write_data(&mut out);
        out.push_str(">>END_MODULE\n");
    }
    out
}

#[must_use]
pub fn summary(modules: &[Box<dyn QcModule>], filename: &str) -> String {
    let mut out = String::with_capacity(512);
    for m in modules {
        let _ = writeln!(out, "{}\t{}\t{}", m.status().as_str(), m.name(), filename);
    }
    out
}
