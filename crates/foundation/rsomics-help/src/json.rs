//! Structured JSON representation of a binary's CLI surface.
//!
//! Emitted by `--help --json`. The goal is an AI consumer reads this
//! one document and can synthesise a correct invocation or config
//! file without parsing `--help` ANSI.
//!
//! Schema version: `1.0`. Bump MAJOR when removing / renaming fields,
//! MINOR for additive optional fields. Consumers should accept any
//! MINOR within their pinned MAJOR.

use serde::Serialize;

pub const HELP_SCHEMA_VERSION: &str = "1.0";

#[derive(Debug, Serialize)]
pub struct HelpJson<'a> {
    pub schema_version: &'static str,
    pub tool: &'a str,
    pub tool_version: &'a str,
    pub tagline: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub origin: Option<Origin<'a>>,
    pub flag_groups: Vec<FlagGroup<'a>>,
    pub examples: Vec<Example<'a>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub json_result_schema_doc: Option<&'a str>,
}

impl<'a> HelpJson<'a> {
    #[must_use]
    pub fn new(tool: &'a str, version: &'a str, tagline: &'a str) -> Self {
        Self {
            schema_version: HELP_SCHEMA_VERSION,
            tool,
            tool_version: version,
            tagline,
            origin: None,
            flag_groups: Vec::new(),
            examples: Vec::new(),
            json_result_schema_doc: None,
        }
    }
}

#[derive(Debug, Serialize)]
pub struct Origin<'a> {
    pub upstream: &'a str,
    pub upstream_license: &'a str,
    pub our_license: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub paper_doi: Option<&'a str>,
}

#[derive(Debug, Serialize)]
pub struct FlagGroup<'a> {
    pub title: &'a str,
    pub flags: Vec<FlagSpec<'a>>,
}

#[derive(Debug, Serialize)]
pub struct FlagSpec<'a> {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub short: Option<char>,
    pub long: &'a str,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub aliases: Vec<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<&'a str>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub type_hint: Option<&'a str>,
    pub required: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub default: Option<&'a str>,
    pub description: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub why_default: Option<&'a str>,
}

#[derive(Debug, Serialize)]
pub struct Example<'a> {
    pub description: &'a str,
    pub command: &'a str,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimal_help_serialises() {
        let h = HelpJson::new("rsomics-x", "0.1.0", "test");
        let s = serde_json::to_string(&h).unwrap();
        assert!(s.contains("\"schema_version\":\"1.0\""));
        assert!(s.contains("\"tool\":\"rsomics-x\""));
    }

    #[test]
    fn flag_spec_serialises_only_present_fields() {
        let f = FlagSpec {
            short: Some('i'),
            long: "in1",
            aliases: vec!["in-1"],
            value: Some("<path>"),
            type_hint: Some("PathBuf"),
            required: true,
            default: None,
            description: "R1 input",
            why_default: None,
        };
        let s = serde_json::to_string(&f).unwrap();
        assert!(s.contains("\"short\":\"i\""));
        assert!(s.contains("\"aliases\":[\"in-1\"]"));
        assert!(!s.contains("\"default\""));
        assert!(!s.contains("\"why_default\""));
    }
}
