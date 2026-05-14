//! figlet-rs banner with per-line linear gradient.
//!
//! The whole-family banner is auto-generated from the binary name —
//! no per-crate ASCII art files. Customisation lives in two axes:
//! the figlet font (defaulted to "slant" for the family) and the
//! gradient (defaulted to teal → violet, the family palette).

use figlet_rs::FIGlet;

use crate::ansi::{Palette, rgb};

/// Linear two-stop gradient in 24-bit RGB.
#[derive(Debug, Clone, Copy)]
pub struct Gradient {
    pub start: (u8, u8, u8),
    pub end: (u8, u8, u8),
}

impl Gradient {
    #[must_use]
    pub const fn new(start: (u8, u8, u8), end: (u8, u8, u8)) -> Self {
        Self { start, end }
    }

    /// Family-default: teal at the top → violet at the bottom.
    #[must_use]
    pub const fn family_default() -> Self {
        Self::new(Palette::TEAL, Palette::VIOLET)
    }

    /// Interpolate at `t` ∈ `[0.0, 1.0]`.
    #[must_use]
    pub fn at(&self, t: f32) -> (u8, u8, u8) {
        let t = t.clamp(0.0, 1.0);
        let lerp = |a: u8, b: u8| -> u8 {
            let v = f32::from(a) + (f32::from(b) - f32::from(a)) * t;
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            let rounded = v.round().clamp(0.0, 255.0) as u8;
            rounded
        };
        (
            lerp(self.start.0, self.end.0),
            lerp(self.start.1, self.end.1),
            lerp(self.start.2, self.end.2),
        )
    }
}

/// One of figlet-rs's built-in fonts.
#[derive(Debug, Clone, Copy)]
pub enum Font {
    Standard,
    Slant,
    Small,
    Big,
}

impl Font {
    fn load(self) -> FIGlet {
        match self {
            Self::Standard => FIGlet::standard(),
            Self::Slant => FIGlet::slant(),
            Self::Small => FIGlet::small(),
            Self::Big => FIGlet::big(),
        }
        .expect("built-in figlet font load is infallible")
    }
}

/// Banner spec: what to render, in what font, with what gradient.
#[derive(Debug, Clone)]
pub struct Banner {
    pub text: String,
    pub font: Font,
    pub gradient: Gradient,
}

impl Banner {
    /// Family-default constructor: slant font, family gradient.
    #[must_use]
    pub fn family(text: impl Into<String>) -> Self {
        Self {
            text: text.into(),
            font: Font::Slant,
            gradient: Gradient::family_default(),
        }
    }

    /// Render the banner. When `color` is `false`, ANSI escapes are
    /// omitted (plain ASCII art only).
    #[must_use]
    pub fn render(&self, color: bool) -> String {
        let figure = self
            .font
            .load()
            .convert(&self.text)
            .map_or_else(|| self.text.clone(), |f| f.to_string());
        let lines: Vec<&str> = figure.lines().collect();
        if lines.is_empty() {
            return String::new();
        }
        // Drop trailing empty lines that figlet pads with.
        let last_non_empty = lines
            .iter()
            .rposition(|l| !l.trim().is_empty())
            .unwrap_or(lines.len() - 1);
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
        painted.join("\n")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gradient_endpoints() {
        let g = Gradient::new((0, 0, 0), (255, 255, 255));
        assert_eq!(g.at(0.0), (0, 0, 0));
        assert_eq!(g.at(1.0), (255, 255, 255));
        let mid = g.at(0.5);
        assert!(mid.0 >= 127 && mid.0 <= 128);
    }

    #[test]
    fn banner_renders_some_lines() {
        let b = Banner::family("rs");
        let rich = b.render(true);
        let plain = b.render(false);
        assert!(!rich.is_empty());
        assert!(!plain.is_empty());
        // Color version must contain at least one escape sequence.
        assert!(rich.contains("\x1b["));
        // Plain version must contain no escape sequences.
        assert!(!plain.contains("\x1b["));
    }
}
