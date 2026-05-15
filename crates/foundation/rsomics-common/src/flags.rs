use std::sync::OnceLock;
use std::thread;

use clap::Args;
#[cfg(feature = "rayon")]
use rayon::ThreadPoolBuilder;

use crate::error::Result;
#[cfg(feature = "rayon")]
use crate::error::RsomicsError;

/// Common flags flattened into every tool's top-level `clap::Parser` via
/// `#[command(flatten)]`. Each flag carries `global = true` so clap
/// propagates it through subcommands — flatten at the top level only.
#[derive(Debug, Clone, Args)]
pub struct CommonFlags {
    /// Number of worker threads to use. Defaults to available parallelism.
    /// Sets the global rayon pool size; whether a tool actually parallelises
    /// over it depends on the tool — check its `--help` for which paths
    /// consume the pool.
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
    #[must_use]
    pub fn thread_count(&self) -> usize {
        self.threads
            .or_else(|| thread::available_parallelism().ok().map(Into::into))
            .unwrap_or(1)
    }

    #[cfg(feature = "rayon")]
    #[allow(clippy::missing_errors_doc)]
    pub fn install_rayon_pool(&self) -> Result<()> {
        let want = self.thread_count();
        if ThreadPoolBuilder::new()
            .num_threads(want)
            .build_global()
            .is_err()
        {
            let active = rayon::current_num_threads();
            if active != want {
                return Err(RsomicsError::ConfigError(format!(
                    "rayon pool already initialised with {active} threads; \
                     cannot reconfigure to {want}"
                )));
            }
        }
        Ok(())
    }

    #[cfg(not(feature = "rayon"))]
    #[allow(clippy::unused_self)]
    pub fn install_rayon_pool(&self) -> Result<()> {
        Ok(())
    }

    /// `--seed 0` is preserved verbatim — explicit zero is not treated as "no seed".
    #[must_use]
    pub fn seed_rng(&self) -> u64 {
        static FRESH_SEED: OnceLock<u64> = OnceLock::new();
        if let Some(s) = self.seed {
            return s;
        }
        *FRESH_SEED.get_or_init(fresh_os_seed)
    }
}

#[allow(clippy::cast_possible_truncation)]
fn fresh_os_seed() -> u64 {
    use std::time::{SystemTime, UNIX_EPOCH};
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos() as u64);
    let pid = u64::from(std::process::id());
    nanos.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ pid.wrapping_mul(0xBF58_476D_1CE4_E5B9)
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
    }

    #[test]
    fn explicit_seed_zero_is_preserved() {
        let cli = Cli::parse_from(["test", "--seed", "0"]);
        assert_eq!(cli.common.seed_rng(), 0, "--seed 0 must round-trip");
    }

    #[test]
    fn long_forms_parse_json_flag() {
        let cli = Cli::parse_from(["test", "--json"]);
        assert!(cli.common.json);
    }
}
