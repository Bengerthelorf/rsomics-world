#!/usr/bin/env python3
"""Deterministic synthetic fixtures for perfgate. Tier-3 (HDD, not git).

    mkfixture.py fasta   OUT n_records seq_len
    mkfixture.py fastq   OUT n_reads read_len     # 3' TruSeq adapter on ~60%
    mkfixture.py fastqgz OUT n_reads read_len     # same bytes, gzip (mtime=0)
    mkfixture.py bed     OUT n_intervals n_chroms  # sorted, overlapping

Fixed seed → byte-identical across runs, so a fixture's sha256 is a
stable identity recorded by perfgate.
"""
import gzip
import random
import sys

KIND, OUT = sys.argv[1], sys.argv[2]
A, B = int(sys.argv[3]), int(sys.argv[4])
# Optional 5th arg = seed, so a two-input tool gets distinct a/b fixtures.
SEED = int(sys.argv[5]) if len(sys.argv) > 5 else 0x00C0FFEE
random.seed(SEED)
ACGT = b"ACGT"
ADAPTER = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

if KIND == "fasta":
    with open(OUT, "wb") as f:
        for i in range(A):
            f.write(f">contig_{i}\n".encode())
            f.write(bytes(ACGT[random.getrandbits(2)] for _ in range(B)) + b"\n")

elif KIND == "fastq":
    with open(OUT, "wb") as f:
        for i in range(A):
            insert = random.randint(B // 3, B)
            seq = bytes(ACGT[random.getrandbits(2)] for _ in range(insert))
            if random.random() < 0.6:
                seq = (seq + ADAPTER)[:B]
            seq = (seq + bytes(ACGT[random.getrandbits(2)]
                               for _ in range(B)))[:B]
            q = bytes(33 + min(40, 20 + random.randint(-8, 15))
                      for _ in range(B))
            f.write(b"@r%d\n%s\n+\n%s\n" % (i, seq, q))

elif KIND == "fastqgz":
    # Byte-identical uncompressed content to `fastq` for the same seed/params,
    # then single-member gzip. mtime=0 keeps the header fixed so the fixture
    # sha256 is a stable identity. level 6 = typical real .fastq.gz.
    buf = bytearray()
    for i in range(A):
        insert = random.randint(B // 3, B)
        seq = bytes(ACGT[random.getrandbits(2)] for _ in range(insert))
        if random.random() < 0.6:
            seq = (seq + ADAPTER)[:B]
        seq = (seq + bytes(ACGT[random.getrandbits(2)]
                           for _ in range(B)))[:B]
        q = bytes(33 + min(40, 20 + random.randint(-8, 15))
                  for _ in range(B))
        buf += b"@r%d\n%s\n+\n%s\n" % (i, seq, q)
    with open(OUT, "wb") as f:
        f.write(gzip.compress(bytes(buf), compresslevel=6, mtime=0))

elif KIND == "bed":
    rows = []
    for _ in range(A):
        c = random.randint(1, B)
        s = random.randint(0, 250_000_000)
        rows.append((f"chr{c}", s, s + random.randint(50, 5000)))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(OUT, "w") as f:
        f.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in rows)

elif KIND == "genome":
    # A = n_chroms. Lexicographic chrom order matches the `bed` kind's
    # sort so `bedtools complement -g` accepts the pairing. Size exceeds
    # the bed kind's max coordinate (250M + 5k).
    chroms = sorted(f"chr{i}" for i in range(1, A + 1))
    with open(OUT, "w") as f:
        f.writelines(f"{c}\t300000000\n" for c in chroms)

elif KIND == "bam":
    # A = n_records, B = read_len. Generates a coordinate-sorted BAM via
    # samtools (must be on PATH). Paired-end, 2 chroms, ~50% properly paired.
    import subprocess, tempfile, os
    n_chroms = 2
    sam_lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        *[f"@SQ\tSN:chr{c}\tLN:300000000" for c in range(1, n_chroms + 1)],
    ]
    for i in range(A):
        chrom = f"chr{random.randint(1, n_chroms)}"
        pos = random.randint(1, 299_000_000)
        seq = "".join(chr(ACGT[random.getrandbits(2)]) for _ in range(B))
        qual = "".join(chr(33 + min(40, 20 + random.randint(-8, 15))) for _ in range(B))
        flag = 99 if random.random() < 0.5 else 0
        sam_lines.append(f"r{i}\t{flag}\t{chrom}\t{pos}\t60\t{B}M\t*\t0\t0\t{seq}\t{qual}")
    sam_content = "\n".join(sam_lines) + "\n"
    with tempfile.NamedTemporaryFile(suffix=".sam", mode="w", delete=False) as tmp:
        tmp.write(sam_content)
        tmp_name = tmp.name
    subprocess.run(["samtools", "view", "-bS", tmp_name, "-o", OUT], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", OUT, "-o", OUT], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "index", OUT], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.unlink(tmp_name)

elif KIND == "vcf":
    # A = n_variants, B = n_samples. Generates a minimal VCF with random SNPs.
    chroms = [f"chr{c}" for c in range(1, 3)]
    with open(OUT, "w") as f:
        f.write("##fileformat=VCFv4.3\n")
        for c in chroms:
            f.write(f"##contig=<ID={c},length=300000000>\n")
        samples = [f"SAMPLE{s}" for s in range(1, B + 1)]
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in samples:
            f.write(f"\t{s}")
        f.write("\n")
        for i in range(A):
            chrom = chroms[random.randint(0, len(chroms) - 1)]
            pos = random.randint(1, 299_000_000)
            ref = chr(ACGT[random.getrandbits(2)])
            alt = chr(ACGT[random.getrandbits(2)])
            while alt == ref:
                alt = chr(ACGT[random.getrandbits(2)])
            qual = random.randint(1, 99)
            filt = "PASS" if qual > 30 else "LowQual"
            gts = "\t".join(
                f"{random.randint(0,1)}/{random.randint(0,1)}:{random.randint(5,50)}"
                for _ in samples
            )
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filt}\t.\tGT:DP\t{gts}\n")

else:
    sys.exit(f"unknown kind {KIND}")
print(f"wrote {OUT}")
