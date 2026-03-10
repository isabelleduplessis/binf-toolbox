#!/usr/bin/env python3
import sys
import pysam
import statistics

def main():
    if len(sys.argv) != 2:
        print("Usage: python bam_readlen_stats.py <input.bam>", file=sys.stderr)
        sys.exit(1)

    bam_path = sys.argv[1]

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        lengths = sorted(
            read.query_length for read in bam.fetch()
            if not read.is_unmapped
            and not read.is_secondary
            and not read.is_supplementary
        )

    if not lengths:
        print("No mapped reads found.", file=sys.stderr)
        sys.exit(1)

    n = len(lengths)
    q1     = lengths[int(n * 0.25)]
    median = statistics.median(lengths)
    q3     = lengths[int(n * 0.75)]

    print("\t".join(["min", "q1", "median", "q3", "max", "mean"]))
    print("\t".join([
        str(lengths[0]),
        str(q1),
        str(median),
        str(q3),
        str(lengths[-1]),
        f"{sum(lengths)/n:.2f}"
    ]))

if __name__ == "__main__":
    main()
