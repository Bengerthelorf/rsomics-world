pub mod error;
pub mod exit;
pub mod flags;
pub mod json;
pub mod log;
pub mod runner;
#[cfg(feature = "test-support")]
pub mod test_support;

pub use error::{Context, Result, RsomicsError};
pub use exit::ExitCode;
pub use flags::CommonFlags;
pub use json::{SCHEMA_VERSION, ToolMeta};
pub use log::StderrLog;
pub use runner::run;
