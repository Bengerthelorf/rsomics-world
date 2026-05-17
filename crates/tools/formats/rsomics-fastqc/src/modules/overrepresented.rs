// u64 counts → f64 only at the final percentage stage.
#![allow(clippy::cast_precision_loss)]

use std::cmp::Reverse;
use std::collections::HashMap;
use std::fmt::Write as _;

use super::{ModuleStatus, QcModule, Record};

/// `FastQC` "Overrepresented sequences" — sequences that are >0.1% of the
/// total are listed; WARN if any sequence exceeds 0.1%, FAIL if any
/// exceeds 1% (clean-room `FastQC` contract). Tracking mirrors the
/// duplication module: only sequences first seen within the first 100 000
/// reads are kept, reads >75 bp keyed by their first 50 bp.
///
/// Each hit is labelled with the best match in `FastQC`'s public
/// `contaminant_list` (a hit must be ≥20 bp with ≤1 mismatch). The
/// pass/warn/fail decision is purely count-based and is exact regardless
/// of the label table.
pub struct OverrepresentedSeqs {
    counts: HashMap<Vec<u8>, u64>,
    seen_reads: u64,
    total: u64,
}

const TRACK_LIMIT: u64 = 100_000;
const KEY_TRUNC_OVER: usize = 75;
const KEY_LEN: usize = 50;

/// `(name, sequence)` — `FastQC`'s public `contaminant_list.txt` verbatim
/// (public config data, not GPL source), used to label the `Possible
/// Source` of an overrepresented sequence. A hit must be ≥20 bp with ≤1
/// mismatch. Names are kept exactly as `FastQC`'s file (including its
/// spelling) so the label column byte-matches. The pass/warn/fail decision
/// is count-based and independent of this table.
const CONTAMINANTS: &[(&str, &[u8])] = &[
    (
        "Illumina Single End Adapter 1",
        b"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "Illumina Single End Adapter 2",
        b"CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    ),
    (
        "Illumina Single End PCR Primer 1",
        b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Single End PCR Primer 2",
        b"CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
    ),
    (
        "Illumina Single End Sequencing Primer",
        b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End Adapter 1",
        b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End Adapter 2",
        b"GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
    ),
    (
        "Illumina Paried End PCR Primer 1",
        b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End PCR Primer 2",
        b"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    ),
    (
        "Illumina Paried End Sequencing Primer 1",
        b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Paired End Sequencing Primer 2",
        b"CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
    ),
    (
        "Illumina DpnII expression Adapter 1",
        b"ACAGGTTCAGAGTTCTACAGTCCGAC",
    ),
    (
        "Illumina DpnII expression Adapter 2",
        b"CAAGCAGAAGACGGCATACGA",
    ),
    (
        "Illumina DpnII expression PCR Primer 1",
        b"CAAGCAGAAGACGGCATACGA",
    ),
    (
        "Illumina DpnII expression PCR Primer 2",
        b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    ),
    (
        "Illumina DpnII expression Sequencing Primer",
        b"CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
    ),
    (
        "Illumina NlaIII expression Adapter 1",
        b"ACAGGTTCAGAGTTCTACAGTCCGACATG",
    ),
    (
        "Illumina NlaIII expression Adapter 2",
        b"CAAGCAGAAGACGGCATACGA",
    ),
    (
        "Illumina NlaIII expression PCR Primer 1",
        b"CAAGCAGAAGACGGCATACGA",
    ),
    (
        "Illumina NlaIII expression PCR Primer 2",
        b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    ),
    (
        "Illumina NlaIII expression Sequencing Primer",
        b"CCGACAGGTTCAGAGTTCTACAGTCCGACATG",
    ),
    (
        "Illumina Small RNA Adapter 1",
        b"GTTCAGAGTTCTACAGTCCGACGATC",
    ),
    ("Illumina Small RNA Adapter 2", b"TGGAATTCTCGGGTGCCAAGG"),
    ("Illumina Small RNA RT Primer", b"CAAGCAGAAGACGGCATACGA"),
    (
        "Illumina Small RNA PCR Primer 2",
        b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    ),
    (
        "Illumina Small RNA Sequencing Primer",
        b"CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
    ),
    ("Illumina Multiplexing Adapter 1", b"GATCGGAAGAGCACACGTCT"),
    (
        "Illumina Multiplexing Adapter 2",
        b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Multiplexing PCR Primer 1.01",
        b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Multiplexing PCR Primer 2.01",
        b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    ),
    (
        "Illumina Multiplexing Read1 Sequencing Primer",
        b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "Illumina Multiplexing Index Sequencing Primer",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
    ),
    (
        "Illumina Multiplexing Read2 Sequencing Primer",
        b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    ),
    (
        "Illumina PCR Primer Index 1",
        b"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 2",
        b"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 3",
        b"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 4",
        b"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 5",
        b"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 6",
        b"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 7",
        b"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 8",
        b"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 9",
        b"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 10",
        b"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 11",
        b"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC",
    ),
    (
        "Illumina PCR Primer Index 12",
        b"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC",
    ),
    (
        "Illumina DpnII Gex Adapter 1",
        b"GATCGTCGGACTGTAGAACTCTGAAC",
    ),
    (
        "Illumina DpnII Gex Adapter 1.01",
        b"ACAGGTTCAGAGTTCTACAGTCCGAC",
    ),
    ("Illumina DpnII Gex Adapter 2", b"CAAGCAGAAGACGGCATACGA"),
    ("Illumina DpnII Gex Adapter 2.01", b"TCGTATGCCGTCTTCTGCTTG"),
    ("Illumina DpnII Gex PCR Primer 1", b"CAAGCAGAAGACGGCATACGA"),
    (
        "Illumina DpnII Gex PCR Primer 2",
        b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    ),
    (
        "Illumina DpnII Gex Sequencing Primer",
        b"CGACAGGTTCAGAGTTCTACAGTCCGACGATC",
    ),
    ("Illumina NlaIII Gex Adapter 1.01", b"TCGGACTGTAGAACTCTGAAC"),
    (
        "Illumina NlaIII Gex Adapter 1.02",
        b"ACAGGTTCAGAGTTCTACAGTCCGACATG",
    ),
    ("Illumina NlaIII Gex Adapter 2.01", b"CAAGCAGAAGACGGCATACGA"),
    ("Illumina NlaIII Gex Adapter 2.02", b"TCGTATGCCGTCTTCTGCTTG"),
    ("Illumina NlaIII Gex PCR Primer 1", b"CAAGCAGAAGACGGCATACGA"),
    (
        "Illumina NlaIII Gex PCR Primer 2",
        b"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
    ),
    (
        "Illumina NlaIII Gex Sequencing Primer",
        b"CCGACAGGTTCAGAGTTCTACAGTCCGACATG",
    ),
    ("Illumina 5p RNA Adapter", b"GTTCAGAGTTCTACAGTCCGACGATC"),
    ("Illumina RNA Adapter1", b"TGGAATTCTCGGGTGCCAAGG"),
    (
        "Illumina Small RNA 3p Adapter 1",
        b"ATCTCGTATGCCGTCTTCTGCTTG",
    ),
    ("Illumina Small RNA PCR Primer 1", b"CAAGCAGAAGACGGCATACGA"),
    (
        "TruSeq Universal Adapter",
        b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    ),
    (
        "TruSeq Adapter, Index 1",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 2",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 3",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 4",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 5",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 6",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 7",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 8",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 9",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 10",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 11",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 12",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 13",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 14",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 15",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 16",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 18",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 19",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 20",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 21",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 22",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 23",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 25",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG",
    ),
    (
        "TruSeq Adapter, Index 27",
        b"GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG",
    ),
    ("Illumina RNA RT Primer", b"GCCTTGGCACCCGAGAATTCCA"),
    (
        "Illumina RNA PCR Primer",
        b"AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA",
    ),
    (
        "RNA PCR Primer, Index 1",
        b"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 2",
        b"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 3",
        b"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 4",
        b"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 5",
        b"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 6",
        b"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 7",
        b"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 8",
        b"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 9",
        b"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 10",
        b"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 11",
        b"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 12",
        b"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 13",
        b"CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 14",
        b"CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 15",
        b"CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 16",
        b"CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 17",
        b"CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 18",
        b"CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 19",
        b"CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 20",
        b"CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 21",
        b"CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 22",
        b"CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 23",
        b"CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 24",
        b"CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 25",
        b"CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 26",
        b"CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 27",
        b"CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 28",
        b"CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 29",
        b"CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 30",
        b"CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 31",
        b"CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 32",
        b"CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 33",
        b"CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 34",
        b"CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 35",
        b"CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 36",
        b"CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 37",
        b"CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 38",
        b"CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 39",
        b"CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 40",
        b"CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 41",
        b"CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 42",
        b"CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 43",
        b"CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 44",
        b"CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 45",
        b"CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 46",
        b"CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 47",
        b"CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    (
        "RNA PCR Primer, Index 48",
        b"CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
    ),
    ("ABI Dynabead EcoP Oligo", b"CTGATCTAGAGGTACCGGATCCCAGCAGT"),
    ("ABI Solid3 Adapter A", b"CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"),
    (
        "ABI Solid3 Adapter B",
        b"CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT",
    ),
    ("ABI Solid3 5' AMP Primer", b"CCACTACGCCTCCGCTTTCCTCTCTATG"),
    ("ABI Solid3 3' AMP Primer", b"CTGCCCCGGGTTCCTCATTCT"),
    ("ABI Solid3 EF1 alpha Sense Primer", b"CATGTGTGTTGAGAGCTTC"),
    (
        "ABI Solid3 EF1 alpha Antisense Primer",
        b"GAAAACCAAAGTGGTCCAC",
    ),
    ("ABI Solid3 GAPDH Forward Primer", b"TTAGCACCCCTGGCCAAGG"),
    ("ABI Solid3 GAPDH Reverse Primer", b"CTTACTCCTTGGAGGCCATG"),
    (
        "Clontech Universal Primer Mix Short",
        b"CTAATACGACTCACTATAGGGC",
    ),
    (
        "Clontech Universal Primer Mix Long",
        b"CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT",
    ),
    (
        "Clontech SMARTer II A Oligonucleotide",
        b"AAGCAGTGGTATCAACGCAGAGTAC",
    ),
    (
        "Clontech SMART CDS Primer II A",
        b"AAGCAGTGGTATCAACGCAGAGTACT",
    ),
];

fn matches_contaminant(seq: &[u8]) -> Option<&'static str> {
    // FastQC rule: slide the shorter sequence along the longer; a hit is a
    // window with ≤1 mismatch. The overlap must be ≥20 bp, except a
    // contaminant shorter than 20 bp (the 12 bp adapters) must match in
    // full. `short` is whichever of (read, contaminant) is shorter, so its
    // whole length is the overlap.
    for &(name, cont) in CONTAMINANTS {
        let (short, long) = if seq.len() <= cont.len() {
            (seq, cont)
        } else {
            (cont, seq)
        };
        let required = cont.len().min(20);
        if short.is_empty() || short.len() < required {
            continue;
        }
        for start in 0..=(long.len() - short.len()) {
            let window = &long[start..start + short.len()];
            let mismatches = window.iter().zip(short).filter(|(a, b)| a != b).count();
            if mismatches <= 1 {
                return Some(name);
            }
        }
    }
    None
}

impl OverrepresentedSeqs {
    #[must_use]
    pub fn new() -> Self {
        Self {
            counts: HashMap::new(),
            seen_reads: 0,
            total: 0,
        }
    }

    fn key(seq: &[u8]) -> &[u8] {
        if seq.len() > KEY_TRUNC_OVER {
            &seq[..KEY_LEN]
        } else {
            seq
        }
    }

    fn worst_fraction(&self) -> f64 {
        if self.total == 0 {
            return 0.0;
        }
        let max = self.counts.values().copied().max().unwrap_or(0);
        max as f64 / self.total as f64
    }
}

impl Default for OverrepresentedSeqs {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for OverrepresentedSeqs {
    fn name(&self) -> &'static str {
        "Overrepresented sequences"
    }

    fn process(&mut self, rec: &Record) {
        self.seen_reads += 1;
        self.total += 1;
        let key = Self::key(rec.seq);
        if let Some(c) = self.counts.get_mut(key) {
            *c += 1;
        } else if self.seen_reads <= TRACK_LIMIT {
            self.counts.insert(key.to_vec(), 1);
        }
    }

    fn finalize(&mut self) {}

    fn status(&self) -> ModuleStatus {
        let f = self.worst_fraction();
        if f > 0.01 {
            ModuleStatus::Fail
        } else if f > 0.001 {
            ModuleStatus::Warn
        } else {
            ModuleStatus::Pass
        }
    }

    fn write_data(&self, out: &mut String) {
        out.push_str("#Sequence\tCount\tPercentage\tPossible Source\n");
        if self.total == 0 {
            return;
        }
        let mut rows: Vec<(&Vec<u8>, u64)> = self
            .counts
            .iter()
            .filter(|(_, c)| **c as f64 / self.total as f64 > 0.001)
            .map(|(s, c)| (s, *c))
            .collect();
        rows.sort_by_key(|r| Reverse(r.1));
        for (seq, count) in rows {
            let pct = count as f64 / self.total as f64 * 100.0;
            let source = matches_contaminant(seq).unwrap_or("No Hit");
            let _ = writeln!(
                out,
                "{}\t{count}\t{pct:.6}\t{source}",
                String::from_utf8_lossy(seq)
            );
        }
    }
}
