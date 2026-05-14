//! Machine-readable JSON envelope for AI / scripting consumers.
//!
//! Every `rsomics-*` binary supports `--json` (workspace-wide bool from
//! [`crate::CommonFlags`]). When set, the binary emits a stable, versioned
//! JSON object on success (to **stdout**) or on error (to **stderr**),
//! mutually exclusive — stdout stays parseable for the success case.
//!
//! Envelope shape (success):
//!
//! ```json
//! {
//!   "schema_version": "1.0",
//!   "tool": "rsomics-fastp",
//!   "tool_version": "0.2.0",
//!   "status": "ok",
//!   "result": { ...tool-specific summary... }
//! }
//! ```
//!
//! Envelope shape (error):
//!
//! ```json
//! {
//!   "schema_version": "1.0",
//!   "tool": "rsomics-fastp",
//!   "tool_version": "0.2.0",
//!   "status": "error",
//!   "error": { "kind": "InvalidInput", "message": "..." },
//!   "exit_code": 1
//! }
//! ```
//!
//! ## Schema versioning
//!
//! `schema_version` is `MAJOR.MINOR`. MINOR bumps add optional fields; MAJOR
//! bumps remove or rename fields. Consumers should accept any MINOR within
//! their pinned MAJOR. The current value is [`SCHEMA_VERSION`].
//!
//! ## Why stdout/stderr split
//!
//! Success JSON to stdout means `rsomics-fastp ... --json | jq '.result'`
//! works directly. Error JSON to stderr means a streaming consumer's stdout
//! parser doesn't have to decide between "is this a partial result or an
//! error?" — errors never appear in the stdout stream.

use std::io::Write;

use serde::Serialize;

use crate::error::RsomicsError;
use crate::exit::ExitCode;

/// Current envelope schema version. Bump MINOR for additive fields,
/// MAJOR for breaking changes (field removal / rename / type change).
pub const SCHEMA_VERSION: &str = "1.0";

/// Tool identity, baked into every envelope. Tools construct one of these
/// from `env!("CARGO_PKG_NAME")` + `env!("CARGO_PKG_VERSION")` and pass it
/// to [`crate::run`].
#[derive(Debug, Clone, Copy)]
pub struct ToolMeta {
    pub name: &'static str,
    pub version: &'static str,
}

#[derive(Serialize)]
struct OkEnvelope<'a, T: Serialize> {
    schema_version: &'static str,
    tool: &'a str,
    tool_version: &'a str,
    status: &'static str,
    result: &'a T,
}

#[derive(Serialize)]
struct ErrorEnvelope<'a> {
    schema_version: &'static str,
    tool: &'a str,
    tool_version: &'a str,
    status: &'static str,
    error: ErrorBody<'a>,
    exit_code: u8,
}

#[derive(Serialize)]
struct ErrorBody<'a> {
    kind: &'static str,
    message: &'a str,
}

/// Emit a success envelope wrapping `result` to stdout. Adds a trailing
/// newline so consumers using `read_until('\n')` see a complete line.
pub fn emit_ok<T: Serialize>(meta: &ToolMeta, result: &T) {
    let env = OkEnvelope {
        schema_version: SCHEMA_VERSION,
        tool: meta.name,
        tool_version: meta.version,
        status: "ok",
        result,
    };
    let _ = serde_json::to_writer(std::io::stdout().lock(), &env);
    let _ = writeln!(std::io::stdout().lock());
}

/// Emit an error envelope to stderr. The `kind` field uses one of the
/// stable variant names (`Io`, `InvalidInput`, `ConfigError`,
/// `UpstreamError`) so scripts can dispatch without parsing the message.
pub fn emit_error(meta: &ToolMeta, err: &RsomicsError) {
    let kind = match err {
        RsomicsError::Io(_) => "Io",
        RsomicsError::InvalidInput(_) => "InvalidInput",
        RsomicsError::ConfigError(_) => "ConfigError",
        RsomicsError::UpstreamError(_) => "UpstreamError",
    };
    let message = err.to_string();
    let exit_code = ExitCode::from(err) as u8;
    let env = ErrorEnvelope {
        schema_version: SCHEMA_VERSION,
        tool: meta.name,
        tool_version: meta.version,
        status: "error",
        error: ErrorBody {
            kind,
            message: &message,
        },
        exit_code,
    };
    let _ = serde_json::to_writer(std::io::stderr().lock(), &env);
    let _ = writeln!(std::io::stderr().lock());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Serialize)]
    struct Summary {
        total: u64,
        passed: u64,
    }

    #[test]
    fn ok_envelope_serializes_with_expected_keys() {
        let meta = ToolMeta {
            name: "rsomics-test",
            version: "0.0.0",
        };
        let s = Summary {
            total: 100,
            passed: 90,
        };
        let env = OkEnvelope {
            schema_version: SCHEMA_VERSION,
            tool: meta.name,
            tool_version: meta.version,
            status: "ok",
            result: &s,
        };
        let v: serde_json::Value =
            serde_json::from_str(&serde_json::to_string(&env).expect("ser")).expect("de");
        assert_eq!(v["schema_version"], "1.0");
        assert_eq!(v["tool"], "rsomics-test");
        assert_eq!(v["status"], "ok");
        assert_eq!(v["result"]["total"], 100);
        assert_eq!(v["result"]["passed"], 90);
    }

    #[test]
    fn error_envelope_has_kind_message_and_exit_code() {
        let meta = ToolMeta {
            name: "rsomics-test",
            version: "0.0.0",
        };
        let err = RsomicsError::InvalidInput("bad header".into());
        let exit = ExitCode::from(&err) as u8;
        let body = ErrorEnvelope {
            schema_version: SCHEMA_VERSION,
            tool: meta.name,
            tool_version: meta.version,
            status: "error",
            error: ErrorBody {
                kind: "InvalidInput",
                message: "invalid input: bad header",
            },
            exit_code: exit,
        };
        let v: serde_json::Value =
            serde_json::from_str(&serde_json::to_string(&body).expect("ser")).expect("de");
        assert_eq!(v["status"], "error");
        assert_eq!(v["error"]["kind"], "InvalidInput");
        assert!(
            v["error"]["message"]
                .as_str()
                .unwrap()
                .contains("bad header")
        );
        assert_eq!(v["exit_code"], exit);
    }
}
