use std::sync::atomic::{AtomicU64, Ordering};
use std::thread;

use clap::Args;
use rayon::ThreadPoolBuilder;

use crate::error::{Result, RsomicsError};

/// Flag block every `rsomics-*` tool flattens into its own `clap::Parser`
/// struct via `#[command(flatten)]`. Holding these in one place keeps short
/// names, help text, and semantics consistent across the family — a user
/// running `rsomics-bam`'s `--threads` learns nothing new vs. `rsomics-fastp`.
#[derive(Debug, Clone, Args)]
pub struct CommonFlags {
    /// Number of worker threads to use. Defaults to available parallelism.
    #[arg(short = 't', long = "threads", global = true)]
    pub threads: Option<usize>,

    /// Emit machine-readable JSON to stdout where applicable.
    #[arg(long = "json", global = true, default_value_t = false)]
    pub json: bool,

    /// Suppress all non-error stderr output.
    #[arg(short = 'q', long = "quiet", global = true, default_value_t = false)]
    pub quiet: bool,

    /// Increase stderr verbosity (debug-level diagnostics).
    #[arg(short = 'v', long = "verbose", global = true, default_value_t = false)]
    pub verbose: bool,

    /// Deterministic seed for any RNG-driven step (sampling, shuffles).
    #[arg(long = "seed", global = true)]
    pub seed: Option<u64>,
}

impl CommonFlags {
    /// Effective thread count. Falls back to `std::thread::available_parallelism`
    /// then to 1 if even that fails.
    #[must_use]
    pub fn thread_count(&self) -> usize {
        self.threads
            .or_else(|| thread::available_parallelism().ok().map(Into::into))
            .unwrap_or(1)
    }

    /// Configure the global rayon pool to match `thread_count()`. Idempotent
    /// inside a single process: calling twice returns `Ok(())` the second
    /// time because rayon won't let us replace an already-installed pool.
    ///
    /// # Errors
    ///
    /// Returns `ConfigError` if rayon refused the build for a reason other
    /// than "pool already initialised".
    pub fn install_rayon_pool(&self) -> Result<()> {
        match ThreadPoolBuilder::new()
            .num_threads(self.thread_count())
            .build_global()
        {
            Ok(()) => Ok(()),
            Err(e) => {
                let msg = e.to_string();
                if msg.contains("already") {
                    Ok(())
                } else {
                    Err(RsomicsError::ConfigError(format!(
                        "rayon pool init failed: {msg}"
                    )))
                }
            }
        }
    }

    /// Seed value for downstream RNGs. If `--seed` was supplied, that
    /// value; otherwise a high-entropy value drawn from the OS once per
    /// process (memoised so repeated calls within one run agree).
    #[must_use]
    pub fn seed_rng(&self) -> u64 {
        static FRESH_SEED: AtomicU64 = AtomicU64::new(0);
        if let Some(s) = self.seed {
            return s;
        }
        let cached = FRESH_SEED.load(Ordering::Relaxed);
        if cached != 0 {
            return cached;
        }
        let s = fresh_os_seed();
        FRESH_SEED
            .compare_exchange(0, s, Ordering::Relaxed, Ordering::Relaxed)
            .map_or_else(|existing| existing, |_| s)
    }

    #[must_use]
    pub fn is_quiet(&self) -> bool {
        self.quiet
    }

    #[must_use]
    pub fn is_json(&self) -> bool {
        self.json
    }

    #[must_use]
    pub fn is_verbose(&self) -> bool {
        self.verbose
    }
}

/// Fold the process id, monotonic time, and a small bit of address-space
/// entropy into a non-zero `u64`. The output isn't cryptographic but is good
/// enough to seed downstream RNGs without pulling in a fresh `getrandom`
/// dependency at this layer.
#[allow(clippy::cast_possible_truncation)]
fn fresh_os_seed() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos() as u64);
    let pid = u64::from(std::process::id());
    let mix = nanos.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ pid.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    if mix == 0 { 1 } else { mix }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser, Debug)]
    struct Cli {
        #[command(flatten)]
        common: CommonFlags,
    }

    #[test]
    fn defaults_when_no_flags() {
        let cli = Cli::parse_from(["test"]);
        assert!(cli.common.threads.is_none());
        assert!(!cli.common.json);
        assert!(!cli.common.quiet);
        assert!(!cli.common.verbose);
        assert!(cli.common.seed.is_none());
        assert!(cli.common.thread_count() >= 1);
    }

    #[test]
    fn short_forms_parse() {
        let cli = Cli::parse_from(["test", "-t", "4", "-q", "-v"]);
        assert_eq!(cli.common.threads, Some(4));
        assert!(cli.common.quiet);
        assert!(cli.common.verbose);
        assert_eq!(cli.common.thread_count(), 4);
    }

    #[test]
    fn explicit_seed_is_used_verbatim() {
        let cli = Cli::parse_from(["test", "--seed", "42"]);
        assert_eq!(cli.common.seed_rng(), 42);
    }

    #[test]
    fn unseeded_runs_produce_stable_seed_within_process() {
        let cli = Cli::parse_from(["test"]);
        let a = cli.common.seed_rng();
        let b = cli.common.seed_rng();
        assert_eq!(a, b, "seed should be memoised inside a single process");
        assert_ne!(a, 0, "fresh seed must be non-zero");
    }

    #[test]
    fn long_forms_parse_json_flag() {
        let cli = Cli::parse_from(["test", "--json"]);
        assert!(cli.common.is_json());
    }
}
