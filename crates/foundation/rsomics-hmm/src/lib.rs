#![allow(
    clippy::cast_precision_loss,
    clippy::missing_errors_doc,
    clippy::missing_panics_doc,
    clippy::similar_names,
    clippy::needless_range_loop,
    clippy::explicit_iter_loop
)]

#[derive(Debug, Clone)]
pub struct Hmm {
    pub n_states: usize,
    pub n_symbols: usize,
    pub pi: Vec<f64>,
    pub trans: Vec<f64>, // row-major: trans[s_from * n_states + s_to]
    pub emit: Vec<f64>,  // row-major: emit[s * n_symbols + sym]
}

#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum HmmError {
    #[error("empty observation sequence")]
    Empty,
    #[error("symbol {sym} out of range (n_symbols={n_symbols})")]
    BadSymbol { sym: usize, n_symbols: usize },
    #[error("model shape mismatch: {detail}")]
    Shape { detail: String },
}

pub type Result<T> = std::result::Result<T, HmmError>;

impl Hmm {
    pub fn new(pi: Vec<f64>, trans: Vec<f64>, emit: Vec<f64>, n_symbols: usize) -> Result<Self> {
        let n_states = pi.len();
        if n_states == 0 {
            return Err(HmmError::Shape {
                detail: "n_states must be > 0".into(),
            });
        }
        if trans.len() != n_states * n_states {
            return Err(HmmError::Shape {
                detail: format!(
                    "trans len {} ≠ n_states * n_states ({})",
                    trans.len(),
                    n_states * n_states
                ),
            });
        }
        if emit.len() != n_states * n_symbols {
            return Err(HmmError::Shape {
                detail: format!(
                    "emit len {} ≠ n_states * n_symbols ({})",
                    emit.len(),
                    n_states * n_symbols
                ),
            });
        }
        Ok(Self {
            n_states,
            n_symbols,
            pi,
            trans,
            emit,
        })
    }

    pub fn viterbi(&self, obs: &[usize]) -> Result<Vec<usize>> {
        if obs.is_empty() {
            return Err(HmmError::Empty);
        }
        for (t, &o) in obs.iter().enumerate() {
            if o >= self.n_symbols {
                return Err(HmmError::BadSymbol {
                    sym: o,
                    n_symbols: self.n_symbols,
                });
            }
            let _ = t;
        }
        let n = self.n_states;
        let t_len = obs.len();
        let mut delta = vec![f64::NEG_INFINITY; n * t_len];
        let mut psi = vec![0_usize; n * t_len];
        let ln = |x: f64| if x > 0.0 { x.ln() } else { f64::NEG_INFINITY };

        for s in 0..n {
            delta[s] = ln(self.pi[s]) + ln(self.emit[s * self.n_symbols + obs[0]]);
        }

        for t in 1..t_len {
            for s_to in 0..n {
                let mut best = f64::NEG_INFINITY;
                let mut argmax = 0;
                for s_from in 0..n {
                    let v = delta[(t - 1) * n + s_from] + ln(self.trans[s_from * n + s_to]);
                    if v > best {
                        best = v;
                        argmax = s_from;
                    }
                }
                delta[t * n + s_to] = best + ln(self.emit[s_to * self.n_symbols + obs[t]]);
                psi[t * n + s_to] = argmax;
            }
        }

        let mut last_best = 0;
        let mut last_val = f64::NEG_INFINITY;
        for s in 0..n {
            let v = delta[(t_len - 1) * n + s];
            if v > last_val {
                last_val = v;
                last_best = s;
            }
        }

        let mut path = vec![0_usize; t_len];
        path[t_len - 1] = last_best;
        for t in (0..t_len - 1).rev() {
            path[t] = psi[(t + 1) * n + path[t + 1]];
        }
        Ok(path)
    }

    pub fn log_likelihood(&self, obs: &[usize]) -> Result<f64> {
        if obs.is_empty() {
            return Err(HmmError::Empty);
        }
        for &o in obs {
            if o >= self.n_symbols {
                return Err(HmmError::BadSymbol {
                    sym: o,
                    n_symbols: self.n_symbols,
                });
            }
        }
        let n = self.n_states;
        let mut alpha = vec![0.0_f64; n];
        for s in 0..n {
            alpha[s] = self.pi[s] * self.emit[s * self.n_symbols + obs[0]];
        }
        let mut sum: f64 = alpha.iter().sum();
        if sum == 0.0 {
            return Ok(f64::NEG_INFINITY);
        }
        let mut log_scale = sum.ln();
        for a in alpha.iter_mut() {
            *a /= sum;
        }

        for &o in &obs[1..] {
            let mut next = vec![0.0_f64; n];
            for s_to in 0..n {
                let mut acc = 0.0;
                for s_from in 0..n {
                    acc += alpha[s_from] * self.trans[s_from * n + s_to];
                }
                next[s_to] = acc * self.emit[s_to * self.n_symbols + o];
            }
            sum = next.iter().sum();
            if sum == 0.0 {
                return Ok(f64::NEG_INFINITY);
            }
            log_scale += sum.ln();
            for a in next.iter_mut() {
                *a /= sum;
            }
            alpha = next;
        }
        Ok(log_scale)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn weather_hmm() -> Hmm {
        // Russell-Norvig umbrella HMM; 0=Rain 1=Sun, 0=Umbrella 1=NoUmbrella.
        let pi = vec![0.5, 0.5];
        let trans = vec![0.7, 0.3, 0.3, 0.7];
        let emit = vec![0.9, 0.1, 0.2, 0.8];
        Hmm::new(pi, trans, emit, 2).unwrap()
    }

    #[test]
    fn viterbi_umbrella_sequence_picks_rain_first_then_sun() {
        let hmm = weather_hmm();
        let path = hmm.viterbi(&[0, 1, 1]).unwrap();
        assert_eq!(path, vec![0, 1, 1]);
    }

    #[test]
    fn viterbi_all_umbrellas_stays_in_rain() {
        let hmm = weather_hmm();
        let path = hmm.viterbi(&[0, 0, 0]).unwrap();
        assert_eq!(path, vec![0, 0, 0]);
    }

    #[test]
    fn log_likelihood_is_finite_and_negative() {
        let hmm = weather_hmm();
        let ll = hmm.log_likelihood(&[0, 1, 1]).unwrap();
        assert!(ll.is_finite());
        assert!(ll < 0.0);
    }

    #[test]
    fn log_likelihood_increases_with_consistent_obs() {
        let hmm = weather_hmm();
        let ll_consistent = hmm.log_likelihood(&[0, 0, 0, 0]).unwrap();
        let ll_random = hmm.log_likelihood(&[0, 1, 0, 1]).unwrap();
        assert!(
            ll_consistent > ll_random,
            "{ll_consistent} should be > {ll_random}"
        );
    }

    #[test]
    fn shape_mismatch_rejected() {
        let r = Hmm::new(
            vec![0.5, 0.5],
            vec![0.7, 0.3, 0.3],
            vec![0.9, 0.1, 0.2, 0.8],
            2,
        );
        assert!(matches!(r, Err(HmmError::Shape { .. })));
    }

    #[test]
    fn bad_symbol_rejected() {
        let hmm = weather_hmm();
        assert!(matches!(
            hmm.viterbi(&[0, 5, 1]),
            Err(HmmError::BadSymbol { .. })
        ));
    }

    #[test]
    fn empty_obs_rejected() {
        let hmm = weather_hmm();
        assert!(matches!(hmm.viterbi(&[]), Err(HmmError::Empty)));
        assert!(matches!(hmm.log_likelihood(&[]), Err(HmmError::Empty)));
    }

    #[test]
    fn viterbi_path_length_matches_obs() {
        let hmm = weather_hmm();
        let path = hmm.viterbi(&[0, 1, 0, 1, 0, 1, 0]).unwrap();
        assert_eq!(path.len(), 7);
    }
}
