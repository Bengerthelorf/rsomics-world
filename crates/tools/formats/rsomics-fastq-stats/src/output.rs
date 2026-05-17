use std::fmt::Write as _;

use crate::compute::{ExtendedStats, FastqStats};

const BASE_HEADERS: &[&str] = &[
    "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len",
];

const EXTENDED_HEADERS: &[&str] = &[
    "Q1", "Q2", "Q3", "sum_gap", "N50", "N50_num", "Q20(%)", "Q30(%)", "AvgQual", "GC(%)", "sum_n",
];

#[must_use]
pub fn render_tabular(s: &FastqStats) -> String {
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

fn write_tab_row(out: &mut String, s: &FastqStats) {
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
        "\t{:.0}\t{:.0}\t{:.0}\t{}\t{}\t{}\t{:.0}\t{:.0}\t{:.2}\t{:.2}\t{}",
        e.q1,
        e.q2,
        e.q3,
        e.sum_gap,
        e.n50,
        e.l50,
        e.q20_percent,
        e.q30_percent,
        e.avg_qual,
        e.gc_percent,
        e.sum_n,
    );
}

#[must_use]
pub fn render_pretty(s: &FastqStats) -> String {
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
        data.push(format!("{:.0}", e.q20_percent));
        data.push(format!("{:.0}", e.q30_percent));
        data.push(format!("{:.2}", e.avg_qual));
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

    fn sample() -> FastqStats {
        FastqStats {
            file: "tiny.fq".into(),
            format: "FASTQ",
            seq_type: SeqType::Dna,
            num_seqs: 3,
            sum_len: 28,
            min_len: 6,
            max_len: 12,
            avg_len: 28.0 / 3.0,
            extended: None,
        }
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
        assert_eq!(cells[1], "FASTQ");
        assert_eq!(cells[6], "9.3");
    }

    #[test]
    fn tabular_extended_column_order_matches_seqkit() {
        let mut s = sample();
        s.extended = Some(ExtendedStats {
            q1: 8.0,
            q2: 10.0,
            q3: 11.0,
            sum_gap: 0,
            n50: 10,
            l50: 2,
            q20_percent: 92.86,
            q30_percent: 92.86,
            avg_qual: 13.45,
            gc_percent: 67.857,
            sum_n: 2,
        });
        let out = render_tabular(&s);
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(
            lines[0],
            "file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\tsum_n"
        );
        let cells: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(cells.len(), 19);
        assert_eq!(cells[14], "93"); // Q20(%) — %.0f of 92.86
        assert_eq!(cells[15], "93"); // Q30(%)
        assert_eq!(cells[16], "13.45"); // AvgQual — %.2f
        assert_eq!(cells[17], "67.86"); // GC(%) — %.2f
        assert_eq!(cells[18], "2"); // sum_n
    }

    #[test]
    fn pretty_keeps_columns_aligned() {
        let out = render_pretty(&sample());
        let lines: Vec<&str> = out.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(lines[0].contains("file"));
        assert!(lines[1].contains("FASTQ"));
    }
}
