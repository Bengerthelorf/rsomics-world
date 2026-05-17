pub(crate) const RESET: &str = "\x1b[0m";
pub(crate) const BOLD: &str = "\x1b[1m";
pub(crate) const DIM: &str = "\x1b[2m";

#[must_use]
pub fn no_color_env() -> bool {
    std::env::var_os("NO_COLOR").is_some()
}

pub struct Palette;

impl Palette {
    pub const TEAL: (u8, u8, u8) = (96, 218, 220);
    pub const SKY: (u8, u8, u8) = (96, 165, 250);
    pub const VIOLET: (u8, u8, u8) = (167, 139, 250);
    pub const MAGENTA: (u8, u8, u8) = (217, 119, 233);
    pub const GREEN: (u8, u8, u8) = (74, 222, 128);
    pub const YELLOW: (u8, u8, u8) = (250, 204, 21);
    pub const SLATE: (u8, u8, u8) = (148, 163, 184);

    // RGB approximations of ANSI-16 \x1b[96m → [36m → [94m → [95m
    pub const FAMILY_GRADIENT: [(u8, u8, u8); 4] = [
        (85, 255, 255),
        (0, 205, 205),
        (95, 135, 255),
        (255, 95, 255),
    ];
}

#[must_use]
pub fn rgb(color: bool, rgb: (u8, u8, u8), s: &str) -> String {
    if !color {
        return s.to_string();
    }
    let (r, g, b) = rgb;
    format!("\x1b[38;2;{r};{g};{b}m{s}{RESET}")
}

#[must_use]
pub fn bold(color: bool, s: &str) -> String {
    if color {
        format!("{BOLD}{s}{RESET}")
    } else {
        s.to_string()
    }
}

#[must_use]
pub fn dim(color: bool, s: &str) -> String {
    if color {
        format!("{DIM}{s}{RESET}")
    } else {
        s.to_string()
    }
}
