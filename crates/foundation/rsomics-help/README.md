# rsomics-help

Layer-A foundation crate: family-wide `--help` renderer for every
`rsomics-*` CLI.

Binaries declare a `HelpSpec` data literal — name, version, tagline,
usage lines, flag tables, examples — and hand it to `render(&spec, mode)`.
This crate owns all formatting decisions: figlet banner with
width-adaptive font + pikpaktui-style cyan→magenta gradient, ANSI section
headers, three output modes (Rich / Plain / JSON), and the
`detect_mode` priority chain (NO_COLOR → !isatty → --plain → --json).

```rust
use rsomics_help::{HelpSpec, Section, FlagSpec, Example, HelpMode, intercept_help, render};

const HELP: HelpSpec = HelpSpec {
    name: "rsomics-example",
    version: "0.1.0",
    tagline: "Example tool.",
    origin: None,
    usage_lines: &["[OPTIONS] <FILE>"],
    sections: &[Section { title: "OPTIONS", flags: &[/* … */] }],
    examples: &[Example { description: "basic", command: "rsomics-example a.fa" }],
    json_result_schema_doc: None,
};

fn main() -> std::process::ExitCode {
    let argv: Vec<String> = std::env::args().collect();
    if let Some(mode) = intercept_help(&argv) {
        render(&HELP, mode);
        return std::process::ExitCode::SUCCESS;
    }
    // … real CLI parsing and pipeline …
    std::process::ExitCode::SUCCESS
}
```

## Modes

| Mode | When | Output |
|---|---|---|
| `Rich` | TTY + colour OK | figlet banner + 24-bit RGB gradient + bold section headers |
| `Plain` | `--plain`, `NO_COLOR`, or stdout-isn't-a-tty | compact ASCII, em-dash tagline, ALL-CAPS sections, no colour |
| `Json` | `--json` | structured envelope with `schema_version` for AI consumers |

## External deps (4-quadrant classification)

- `figlet-rs` — Quadrant ① (pure Rust, no FFI).
- `terminal_size` — Quadrant ①.
- `serde`, `serde_json` — Quadrant ④.

License: MIT OR Apache-2.0.
