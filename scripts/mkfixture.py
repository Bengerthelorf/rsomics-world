#!/usr/bin/env python3
"""Deterministic synthetic fixtures for perfgate. Tier-3 (HDD, not git).

    mkfixture.py fasta OUT n_records seq_len
    mkfixture.py fastq OUT n_reads read_len      # 3' TruSeq adapter on ~60%
    mkfixture.py bed   OUT n_intervals n_chroms   # sorted, overlapping

Fixed seed → byte-identical across runs, so a fixture's sha256 is a
stable identity recorded by perfgate.
"""
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

elif KIND == "bed":
    rows = []
    for _ in range(A):
        c = random.randint(1, B)
        s = random.randint(0, 250_000_000)
        rows.append((f"chr{c}", s, s + random.randint(50, 5000)))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(OUT, "w") as f:
        f.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in rows)

else:
    sys.exit(f"unknown kind {KIND}")
print(f"wrote {OUT}")
