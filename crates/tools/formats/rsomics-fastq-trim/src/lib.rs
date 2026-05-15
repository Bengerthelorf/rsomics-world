pub mod adapter;
pub mod fixed;
pub mod overlap;
pub(crate) mod parallel_gz;
pub mod pipeline;
pub mod polyx;

pub use adapter::{AdapterConfig, find_adapter_3p};
pub use fixed::{FixedTrimConfig, apply_fixed};
pub use overlap::{OverlapConfig, OverlapResult, analyze as analyze_overlap, reverse_complement};
pub use pipeline::{Pipeline, PipelineConfig, TrimReport};
pub use polyx::{PolyXConfig, find_dominant_polyx_3p, find_polyx_3p};
