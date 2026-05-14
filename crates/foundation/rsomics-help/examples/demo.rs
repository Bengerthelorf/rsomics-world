use rsomics_help::{Banner, example_line, flag_row, section_header, tagline};

fn main() {
    let color = true;
    println!();
    println!("{}", Banner::family("rsomics-fastq-trim").render(color));
    println!();
    println!(
        "  {}",
        tagline(
            "rsomics-fastq-trim",
            "0.2.0",
            "FASTQ adapter / polyG / polyX / fixed-length trimming",
            color,
        )
    );
    println!();
    println!("{}", section_header("USAGE", color));
    println!("  rsomics-fastq-trim [OPTIONS] --in1 <PATH> --out1 <PATH>");
    println!();
    println!("{}", section_header("OPTIONS", color));
    println!(
        "{}",
        flag_row(Some('i'), "in1", Some("<path>"), "R1 input", color, 26)
    );
    println!(
        "{}",
        flag_row(Some('o'), "out1", Some("<path>"), "R1 output", color, 26)
    );
    println!(
        "{}",
        flag_row(
            None,
            "adapter_min_len",
            Some("<n>"),
            "min match (default 5, matches fastp)",
            color,
            26,
        )
    );
    println!(
        "{}",
        flag_row(
            None,
            "compression",
            Some("<lvl>"),
            "libdeflate level 1-12 (default 4)",
            color,
            26,
        )
    );
    println!();
    println!("{}", section_header("EXAMPLES", color));
    println!(
        "{}",
        example_line(
            "Adapter trim a single-end gz",
            "rsomics-fastq-trim -i in.fq.gz -o out.fq.gz",
            color,
        )
    );
    println!();
}
