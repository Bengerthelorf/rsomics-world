pub mod compute;
pub mod output;

pub use compute::{Config, FastqStats, SeqType, compute_stats};
pub use output::{render_pretty, render_tabular};
