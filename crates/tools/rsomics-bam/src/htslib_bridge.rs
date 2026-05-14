//! Local bridge between `rust-htslib`'s error type and our own
//! [`RsomicsError`]. Lives in rsomics-bam (not rsomics-common) per the
//! "promote to Layer A only when 2+ B crates need it" rule — until a
//! second htslib-using crate exists (rsomics-bcftools, an htslib-based
//! variant tool, etc.), the conversion is local.

use rsomics_common::RsomicsError;

/// Map any `rust-htslib` error into [`RsomicsError::UpstreamError`]. The
/// resulting message preserves htslib's full Display chain so the
/// underlying C-side reason survives one layer of context.
///
/// Takes the error by value so `.map_err(from_htslib)` works directly —
/// [`std::result::Result::map_err`] passes the error owned, not by
/// reference. Don't refactor this signature without also touching
/// [`HtsResultExt::rs`].
#[must_use]
#[allow(clippy::needless_pass_by_value)]
pub fn from_htslib(e: rust_htslib::errors::Error) -> RsomicsError {
    RsomicsError::UpstreamError(format!("htslib: {e}"))
}

/// Extension trait that gives any `Result<T, rust_htslib::errors::Error>`
/// a `.rs()?` that converts to [`rsomics_common::Result`]. Avoids
/// peppering call sites with `.map_err(htslib_bridge::from_htslib)?`.
pub trait HtsResultExt<T> {
    /// Convert to our `Result<T>`.
    ///
    /// # Errors
    ///
    /// Returns whatever `rust-htslib` error was inside, wrapped as
    /// [`RsomicsError::UpstreamError`].
    fn rs(self) -> rsomics_common::Result<T>;
}

impl<T> HtsResultExt<T> for std::result::Result<T, rust_htslib::errors::Error> {
    fn rs(self) -> rsomics_common::Result<T> {
        self.map_err(from_htslib)
    }
}
