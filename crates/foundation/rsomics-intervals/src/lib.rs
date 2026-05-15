
pub mod algebra;
pub mod bed;
pub mod index;
pub mod interval;
pub mod set;

pub use algebra::{complement, coverage_bases, intersect, merge, subtract};
pub use index::IntervalIndex;
pub use interval::{Interval, IntervalError, Strand};
pub use set::IntervalSet;
