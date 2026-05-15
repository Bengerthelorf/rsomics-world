
use serde::Serialize;

pub const HELP_SCHEMA_VERSION: &str = "1.0";

#[derive(Debug, Serialize)]
pub struct HelpSpec<'a> {
    #[serde(rename = "tool")]
    pub name: &'a str,
    #[serde(rename = "tool_version")]
    pub version: &'a str,
    pub tagline: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub origin: Option<Origin<'a>>,
    pub usage_lines: &'a [&'a str],
    pub sections: &'a [Section<'a>],
    pub examples: &'a [Example<'a>],
    #[serde(skip_serializing_if = "Option::is_none")]
    pub json_result_schema_doc: Option<&'a str>,
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
pub struct Section<'a> {
    pub title: &'a str,
    pub flags: &'a [FlagSpec<'a>],
}

#[derive(Debug, Serialize)]
pub struct FlagSpec<'a> {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub short: Option<char>,
    pub long: &'a str,
    #[serde(skip_serializing_if = "slice_is_empty")]
    pub aliases: &'a [&'a str],
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

#[allow(clippy::ref_option)]
fn slice_is_empty<T>(s: &&[T]) -> bool {
    s.is_empty()
}
