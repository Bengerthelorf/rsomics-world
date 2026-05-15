//! ANSI escape primitives + `NO_COLOR` detection.
//!
//! 24-bit foreground colors via `\x1b[38;2;R;G;Bm` (universally supported
//! in modern terminals; degrades visibly but harmlessly elsewhere). The
//! family palette is bound here so every binary's `--help` lands on the
//! same hues.

pub const RESET: &str = "\x1b[0m";
pub const BOLD: &str = "\x1b[1m";
pub const DIM: &str = "\x1b[2m";
pub const ITALIC: &str = "\x1b[3m";
pub const UNDERLINE: &str = "\x1b[4m";

/// Returns `true` if the user has set `NO_COLOR` to any value. Per
/// <https://no-color.org/>, that's an environment-wide opt-out.
#[must_use]
pub fn no_color_env() -> bool {
    std::env::var_os("NO_COLOR").is_some()
}

/// Family color palette. TEAL anchors taglines / section colors; the
/// banner uses a separate 4-stop progression that mirrors pikpaktui's
/// bright-cyan → cyan → bright-blue → bright-magenta hue walk.
pub struct Palette;

impl Palette {
    pub const TEAL: (u8, u8, u8) = (96, 218, 220);
    pub const SKY: (u8, u8, u8) = (96, 165, 250);
    pub const VIOLET: (u8, u8, u8) = (167, 139, 250);
    pub const MAGENTA: (u8, u8, u8) = (217, 119, 233);
    pub const GREEN: (u8, u8, u8) = (74, 222, 128);
    pub const YELLOW: (u8, u8, u8) = (250, 204, 21);
    pub const SLATE: (u8, u8, u8) = (148, 163, 184);

    /// Banner gradient stops — RGB approximations of the ANSI-16 sequence
    /// `\x1b[96m → [36m → [94m → [95m` used by pikpaktui's `dMP` art.
    pub const FAMILY_GRADIENT: [(u8, u8, u8); 4] = [
        (85, 255, 255),
        (0, 205, 205),
        (95, 135, 255),
        (255, 95, 255),
    ];
}

/// Format `s` with a 24-bit RGB foreground. When `color` is `false`,
/// returns `s` unchanged (no escapes).
#[must_use]
pub fn rgb(color: bool, rgb: (u8, u8, u8), s: &str) -> String {
    if !color {
        return s.to_string();
    }
    let (r, g, b) = rgb;
    format!("\x1b[38;2;{r};{g};{b}m{s}{RESET}")
}

/// Bold-style wrapper. No-op when `color` is `false`.
#[must_use]
pub fn bold(color: bool, s: &str) -> String {
    if color {
        format!("{BOLD}{s}{RESET}")
    } else {
        s.to_string()
    }
}

/// Dim-style wrapper. No-op when `color` is `false`.
#[must_use]
pub fn dim(color: bool, s: &str) -> String {
    if color {
        format!("{DIM}{s}{RESET}")
    } else {
        s.to_string()
    }
}
