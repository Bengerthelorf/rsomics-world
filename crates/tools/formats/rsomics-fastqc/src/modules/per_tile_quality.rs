// u64 sums/counts → f64 only at the final mean/deviation stage.
#![allow(clippy::cast_precision_loss)]

use std::collections::HashMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

// Phred offset cancels in the mean difference; raw bytes used directly. Non-Illumina IDs ⇒ no tile data ⇒ Pass.
pub struct PerTileQuality {
    tile: HashMap<(u32, usize), (u64, u64)>,
    overall: Vec<(u64, u64)>,
}

impl PerTileQuality {
    #[must_use]
    pub fn new() -> Self {
        Self {
            tile: HashMap::new(),
            overall: Vec::new(),
        }
    }

    // Illumina 1.8+ (≥7 colon fields): tile = field[4]; classic pre-1.8 (5 fields): tile = field[2]
    fn parse_tile(id: &[u8]) -> Option<u32> {
        let head = id.split(|&b| b == b' ').next()?;
        let fields: Vec<&[u8]> = head.split(|&b| b == b':').collect();
        let tile = match fields.len() {
            n if n >= 7 => fields[4],
            5 => fields[2],
            _ => return None,
        };
        std::str::from_utf8(tile).ok()?.parse::<u32>().ok()
    }

    fn worst_deviation(&self) -> f64 {
        let mut worst = 0.0_f64;
        for (&(_, pos), &(s, c)) in &self.tile {
            if c == 0 {
                continue;
            }
            let (os, oc) = self.overall.get(pos).copied().unwrap_or((0, 0));
            if oc == 0 {
                continue;
            }
            let dev = s as f64 / c as f64 - os as f64 / oc as f64;
            if dev < worst {
                worst = dev;
            }
        }
        worst
    }
}

impl Default for PerTileQuality {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerTileQuality {
    fn name(&self) -> &'static str {
        "Per tile sequence quality"
    }

    fn process(&mut self, rec: &Record) {
        let Some(tile) = Self::parse_tile(rec.id) else {
            return;
        };
        if rec.qual.len() > self.overall.len() {
            self.overall.resize(rec.qual.len(), (0, 0));
        }
        for (pos, &q) in rec.qual.iter().enumerate() {
            let e = self.tile.entry((tile, pos)).or_insert((0, 0));
            e.0 += u64::from(q);
            e.1 += 1;
            self.overall[pos].0 += u64::from(q);
            self.overall[pos].1 += 1;
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let w = self.worst_deviation();
        if w < -5.0 {
            ModuleStatus::Fail
        } else if w < -2.0 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Tile\tBase\tMean\n");
        let mut keys: Vec<(u32, usize)> = self.tile.keys().copied().collect();
        keys.sort_unstable();
        for (tile, pos) in keys {
            let (s, c) = self.tile[&(tile, pos)];
            let (os, oc) = self.overall.get(pos).copied().unwrap_or((0, 0));
            if c == 0 || oc == 0 {
                continue;
            }
            let dev = s as f64 / c as f64 - os as f64 / oc as f64;
            let _ = writeln!(out, "{tile}\t{}\t{dev:.6}", pos + 1);
        }
    }
}
