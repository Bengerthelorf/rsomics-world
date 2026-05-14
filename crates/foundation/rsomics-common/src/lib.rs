pub mod error;
pub mod exit;
pub mod flags;
pub mod log;
pub mod runner;

pub use error::{Result, RsomicsError};
pub use exit::ExitCode;
pub use flags::CommonFlags;
pub use log::StderrLog;
pub use runner::run;
