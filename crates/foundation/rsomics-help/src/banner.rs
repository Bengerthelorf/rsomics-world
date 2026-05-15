use figlet_rs::FIGlet;

use crate::ansi::{Palette, rgb};

const FALLBACK_TERM_WIDTH: usize = 80;

#[derive(Debug, Clone, Copy)]
pub struct Gradient<'a> {
    pub stops: &'a [(u8, u8, u8)],
}

impl<'a> Gradient<'a> {
    #[must_use]
    pub const fn new(stops: &'a [(u8, u8, u8)]) -> Self {
        Self { stops }
    }

    #[must_use]
    pub const fn family_default() -> Gradient<'static> {
        Gradient::new(&Palette::FAMILY_GRADIENT)
    }

    #[must_use]
    pub fn at(&self, t: f32) -> (u8, u8, u8) {
        if self.stops.is_empty() {
            return (255, 255, 255);
        }
        if self.stops.len() == 1 {
            return self.stops[0];
        }
        let t = t.clamp(0.0, 1.0);
        let n_segs = self.stops.len() - 1;
        #[allow(clippy::cast_precision_loss)]
        let scaled = t * n_segs as f32;
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let idx = (scaled.floor() as usize).min(n_segs - 1);
        #[allow(clippy::cast_precision_loss)]
        let local_t = scaled - idx as f32;
        let lo = self.stops[idx];
        let hi = self.stops[idx + 1];
        let lerp = |from: u8, to: u8| -> u8 {
            let v = f32::from(from) + (f32::from(to) - f32::from(from)) * local_t;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let rounded = v.round().clamp(0.0, 255.0) as u8;
            rounded
        };
        (lerp(lo.0, hi.0), lerp(lo.1, hi.1), lerp(lo.2, hi.2))
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum Font {
    Slant,
    Small,
}

impl Font {
    fn load(self) -> FIGlet {
        match self {
            Self::Slant => FIGlet::slant(),
            Self::Small => FIGlet::small(),
        }
        .expect("built-in figlet font load is infallible")
    }
}

#[derive(Debug, Clone)]
pub struct Banner<'a> {
    pub text: String,
    pub gradient: Gradient<'a>,
}

impl Banner<'_> {
    #[must_use]
    pub fn family(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            gradient: Gradient::family_default(),
        }
    }

    #[must_use]
    pub fn render(&self, color: bool) -> String {
        let term_width =
            terminal_size::terminal_size().map_or(FALLBACK_TERM_WIDTH, |(w, _)| usize::from(w.0));
        for font in [Font::Slant, Font::Small] {
            if let Some(out) = self.try_render(font, color, term_width) {
                return out;
            }
        }
        String::new()
    }

    fn try_render(&self, font: Font, color: bool, term_width: usize) -> Option<String> {
        let figure = font.load().convert(&self.text)?.to_string();
        let lines: Vec<&str> = figure.lines().collect();
        if lines.is_empty() {
            return None;
        }
        let last_non_empty = lines.iter().rposition(|l| !l.trim().is_empty())?;
        let max_w = lines[..=last_non_empty]
            .iter()
            .map(|l| l.chars().count())
            .max()
            .unwrap_or(0);
        if max_w > term_width {
            return None;
        }
        let painted: Vec<String> = lines[..=last_non_empty]
            .iter()
            .enumerate()
            .map(|(i, line)| {
                let t = if last_non_empty == 0 {
                    0.0
                } else {
                    #[allow(clippy::cast_precision_loss)]
                    let t = i as f32 / last_non_empty as f32;
                    t
                };
                let c = self.gradient.at(t);
                rgb(color, c, line)
            })
            .collect();
        Some(painted.join("\n"))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_endpoints_and_midpoint() {
        let g = Gradient::new(&[(0, 0, 0), (255, 255, 255)]);
        assert_eq!(g.at(0.0), (0, 0, 0));
        assert_eq!(g.at(1.0), (255, 255, 255));
        let mid = g.at(0.5);
        assert!(mid.0 >= 127 && mid.0 <= 128);
    }

    #[test]
    fn gradient_multi_stop_hits_each_stop() {
        let g = Gradient::family_default();
        assert_eq!(g.at(0.0), Palette::FAMILY_GRADIENT[0]);
        assert_eq!(g.at(1.0 / 3.0), Palette::FAMILY_GRADIENT[1]);
        assert_eq!(g.at(2.0 / 3.0), Palette::FAMILY_GRADIENT[2]);
        assert_eq!(g.at(1.0), Palette::FAMILY_GRADIENT[3]);
    }

    #[test]
    fn short_banner_renders_with_color_when_requested() {
        let b = Banner::family("rs");
        let rich = b.render(true);
        let plain = b.render(false);
        assert!(!rich.is_empty());
        assert!(!plain.is_empty());
        assert!(rich.contains("\x1b["));
        assert!(!plain.contains("\x1b["));
    }
}
