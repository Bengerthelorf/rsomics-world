
use std::fmt::Write as _;

use crate::compute::{ExtendedStats, FastaStats};

const BASE_HEADERS: &[&str] = &[
    "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len",
];

const EXTENDED_HEADERS: &[&str] = &[
    "Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num", "GC(%)", "sum_n",
];

#[must_use]
pub fn render_tabular(s: &FastaStats) -> String {
    let mut out = String::with_capacity(256);
    write_tab_header(&mut out, s.extended.is_some());
    write_tab_row(&mut out, s);
    out
}

fn write_tab_header(out: &mut String, extended: bool) {
    let mut first = true;
    for h in BASE_HEADERS {
        if !first {
            out.push('\t');
        }
        first = false;
        out.push_str(h);
    }
    if extended {
        for h in EXTENDED_HEADERS {
            out.push('\t');
            out.push_str(h);
        }
    }
    out.push('\n');
}

fn write_tab_row(out: &mut String, s: &FastaStats) {
    let _ = write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}",
        s.file,
        s.format,
        s.seq_type.as_str(),
        s.num_seqs,
        s.sum_len,
        s.min_len,
        s.avg_len,
        s.max_len,
    );
    if let Some(e) = &s.extended {
        write_extended_row(out, e);
    }
    out.push('\n');
}

fn write_extended_row(out: &mut String, e: &ExtendedStats) {
    let _ = write!(
        out,
        "\t{:.0}\t{:.0}\t{:.0}\t{}\t{}\t{}\t{:.2}\t{}",
        e.q1, e.q2, e.q3, e.sum_gap, e.n50, e.l50, e.gc_percent, e.sum_n,
    );
}

#[must_use]
pub fn render_pretty(s: &FastaStats) -> String {
    let mut rows: Vec<Vec<String>> = Vec::with_capacity(2);

    let mut header: Vec<String> = BASE_HEADERS.iter().map(|h| (*h).to_string()).collect();
    let mut data: Vec<String> = vec![
        s.file.clone(),
        s.format.to_string(),
        s.seq_type.as_str().to_string(),
        humanize_u64(s.num_seqs),
        humanize_u64(s.sum_len),
        humanize_u64(s.min_len),
        format!("{:.1}", s.avg_len),
        humanize_u64(s.max_len),
    ];

    if let Some(e) = &s.extended {
        for h in EXTENDED_HEADERS {
            header.push((*h).to_string());
        }
        data.push(format!("{:.0}", e.q1));
        data.push(format!("{:.0}", e.q2));
        data.push(format!("{:.0}", e.q3));
        data.push(humanize_u64(e.sum_gap));
        data.push(humanize_u64(e.n50));
        data.push(humanize_u64(e.l50));
        data.push(format!("{:.2}", e.gc_percent));
        data.push(humanize_u64(e.sum_n));
    }

    rows.push(header);
    rows.push(data);
    render_columns(&rows)
}

fn render_columns(rows: &[Vec<String>]) -> String {
    let ncols = rows[0].len();
    let mut widths = vec![0usize; ncols];
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            widths[i] = widths[i].max(cell.len());
        }
    }
    let mut out = String::with_capacity(ncols * widths.iter().sum::<usize>() + 32);
    for row in rows {
        for (i, cell) in row.iter().enumerate() {
            if i > 0 {
                out.push_str("  ");
            }
            let _ = write!(out, "{:<width$}", cell, width = widths[i]);
        }
        out.push('\n');
    }
    out
}

fn humanize_u64(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    let first_chunk = s.len() % 3;
    let (head, tail) = s.split_at(first_chunk);
    out.push_str(head);
    for triplet in tail.as_bytes().chunks(3) {
        if !out.is_empty() {
            out.push(',');
        }
        out.push_str(std::str::from_utf8(triplet).expect("decimal digits are ASCII"));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::SeqType;

    fn sample() -> FastaStats {
        FastaStats {
            file: "tiny.fa".into(),
            format: "FASTA",
            seq_type: SeqType::Dna,
            num_seqs: 3,
            sum_len: 14,
            min_len: 4,
            max_len: 6,
            avg_len: 14.0 / 3.0,
            extended: None,
        }
    }

    #[test]
    fn humanize_basic() {
        assert_eq!(humanize_u64(0), "0");
        assert_eq!(humanize_u64(999), "999");
        assert_eq!(humanize_u64(1_000), "1,000");
        assert_eq!(humanize_u64(12_345), "12,345");
        assert_eq!(humanize_u64(1_234_567), "1,234,567");
    }

    #[test]
    fn tabular_basic_shape() {
        let out = render_tabular(&sample());
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(
            lines[0],
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len"
        );
        let cells: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(cells[0], "tiny.fa");
        assert_eq!(cells[1], "FASTA");
        assert_eq!(cells[2], "DNA");
        assert_eq!(cells[3], "3");
        assert_eq!(cells[4], "14");
        assert_eq!(cells[5], "4");
        assert_eq!(cells[6], "4.7");
        assert_eq!(cells[7], "6");
    }

    #[test]
    fn tabular_extended_appends_extra_columns() {
        let mut s = sample();
        s.extended = Some(ExtendedStats {
            q1: 4.0,
            q2: 4.0,
            q3: 6.0,
            sum_gap: 0,
            n50: 6,
            l50: 1,
            gc_percent: 50.0,
            sum_n: 2,
        });
        let out = render_tabular(&s);
        let lines: Vec<&str> = out.lines().collect();
        assert!(lines[0].ends_with("\tsum_n"));
        let cells: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(cells.len(), 16);
        assert_eq!(cells[14], "50.00");
        assert_eq!(cells[15], "2");
    }

    #[test]
    fn pretty_keeps_columns_aligned() {
        let out = render_pretty(&sample());
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].contains("file"));
        assert!(lines[1].contains('3'));
    }
}
