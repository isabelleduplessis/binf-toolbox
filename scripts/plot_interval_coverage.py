#!/usr/bin/env python3
"""
plot_interval_coverage.py

Plot average depth of coverage across intervals in a BED file,
normalizing all intervals to 0-100% of their length (metagene-style).
Intervals are always plotted 5'->3': plus-strand intervals go start->end,
minus-strand intervals go end->start. BED file must have strand in column 6.

Requirements:
    pip install pysam numpy matplotlib pandas

Usage:
    python plot_interval_coverage.py -b reads.bam -i intervals.bed -o outdir/

Outputs (all written to outdir, named after the BAM file stem):
    <stem>_coverage.png
    <stem>_coverage.tsv
    <stem>_coverage_per_interval.tsv

Optional:
    -n / --bins         Number of bins across the normalized interval (default: 100)
    -m / --min-length   Minimum interval length to include (default: 1)
    -q / --mapq         Minimum mapping quality (default: 0)
    --title             Custom plot title
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pysam


def parse_bed(bed_file, min_length=1):
    """Parse BED file, returning list of (chrom, start, end, name, strand) tuples.
    Strand is required (column 6)."""
    intervals = []
    with open(bed_file) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                sys.exit(f"Error: BED file must have at least 6 columns (chrom, start, end, name, score, strand). "
                         f"Line {line_num} has only {len(fields)}. "
                         f"Use: awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,$3,\".\",\".\", $4}}' input.bed > fixed.bed")
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                print(f"Warning: skipping line {line_num} with non-integer coords", file=sys.stderr)
                continue
            if end <= start:
                print(f"Warning: skipping zero/negative-length interval at line {line_num}", file=sys.stderr)
                continue
            if (end - start) < min_length:
                continue
            name   = fields[3]
            strand = fields[5]
            if strand not in ("+", "-"):
                print(f"Warning: unrecognised strand '{strand}' at line {line_num}, skipping", file=sys.stderr)
                continue
            intervals.append((chrom, start, end, name, strand))
    return intervals


def get_coverage_array(bam, chrom, start, end, min_mapq=0):
    """Return a numpy array of per-base depth over [start, end)."""
    length = end - start
    cov = np.zeros(length, dtype=np.float32)
    try:
        for col in bam.pileup(chrom, start, end,
                               truncate=True,
                               min_mapping_quality=min_mapq,
                               ignore_overlaps=False,
                               flag_filter=0x904):  # exclude unmapped, secondary, supplementary
            pos = col.reference_pos - start
            if 0 <= pos < length:
                cov[pos] = col.nsegments
    except (ValueError, KeyError):
        pass
    return cov


def normalize_to_bins(cov_array, n_bins):
    """Resample a coverage array of arbitrary length to exactly n_bins values."""
    length = len(cov_array)
    if length == 0:
        return np.zeros(n_bins)
    x_in  = np.linspace(0, length - 1, length)
    x_out = np.linspace(0, length - 1, n_bins)
    return np.interp(x_out, x_in, cov_array)


def compute_matrix(bam_path, intervals, n_bins=100, min_mapq=0):
    """
    Build a (n_intervals x n_bins) coverage matrix.
    Each interval is oriented 5'->3' before binning:
      plus strand : start -> end  (genomic order)
      minus strand: end   -> start (reversed)
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    rows = []

    for i, (chrom, start, end, name, strand) in enumerate(intervals):
        if (i + 1) % 500 == 0 or i == 0:
            print(f"  Processing interval {i+1}/{len(intervals)}...", file=sys.stderr)

        cov = get_coverage_array(bam, chrom, start, end, min_mapq)

        if strand == "-":
            cov = cov[::-1]

        rows.append(normalize_to_bins(cov, n_bins))

    bam.close()

    if not rows:
        sys.exit("Error: no valid coverage data extracted. Check chromosome names match between BAM and BED.")

    print(f"  Done. Processed {len(rows)} intervals.", file=sys.stderr)
    return np.vstack(rows)


def export_tsv(matrix, n_bins, summary_path, per_interval_path):
    """
    Write two TSVs:
      1. Summary: one row per bin with mean, CI, std across all intervals.
      2. Per-interval: full matrix with one row per interval.
    """
    x    = np.linspace(0, 100, n_bins)
    n    = matrix.shape[0]
    mean = matrix.mean(axis=0)
    std  = matrix.std(axis=0)
    sem  = std / np.sqrt(n)

    summary_df = pd.DataFrame({
        "position_pct" : np.round(x, 4),
        "n_intervals"  : n,
        "mean_depth"   : np.round(mean,            6),
        "ci_lo_95"     : np.round(mean - 1.96*sem, 6),
        "ci_hi_95"     : np.round(mean + 1.96*sem, 6),
        "std_depth"    : np.round(std,             6),
    })

    bin_labels = [f"bin_{i+1}_pct{round(p,1)}" for i, p in enumerate(x)]
    matrix_df  = pd.DataFrame(matrix, columns=bin_labels)

    summary_df.to_csv(summary_path,      sep="\t", index=False)
    matrix_df.to_csv(per_interval_path,  sep="\t", index=False)

    print(f"Summary TSV saved  {summary_path}")
    print(f"Per-interval TSV saved {per_interval_path}")


def plot_coverage(matrix, n_bins, plot_path, title):
    """Single mean line, white background, clean axes."""
    x    = np.linspace(0, 100, n_bins)
    mean = matrix.mean(axis=0)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x, mean, color="steelblue", linewidth=1.5)

    ax.set_xlabel("Position within interval (%)")
    ax.set_ylabel("Mean depth of coverage")
    ax.set_title(title if title else f"Average coverage across {matrix.shape[0]} intervals")

    ax.set_xlim(0, 100)
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter("%g%%"))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(plot_path, dpi=180, bbox_inches="tight")
    print(f"Plot saved             → {plot_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Metagene-style average coverage plot from BAM + BED")
    ap.add_argument("-b", "--bam",        required=True, help="Indexed BAM file")
    ap.add_argument("-i", "--bed",        required=True, help="BED file with strand in column 6")
    ap.add_argument("-o", "--outdir",     required=True,
                    help="Output directory for all results files")
    ap.add_argument("-n", "--bins",       type=int, default=100,
                    help="Number of bins across normalized interval (default: 100)")
    ap.add_argument("-m", "--min-length", type=int, default=1,
                    help="Minimum interval length to include (default: 1)")
    ap.add_argument("-q", "--mapq",       type=int, default=0,
                    help="Minimum mapping quality (default: 0)")
    ap.add_argument("--title",           default=None,
                    help="Custom plot title")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    

    plot_path    = os.path.join(args.outdir, "coverage.png")
    tsv_path     = os.path.join(args.outdir, "coverage.tsv")
    tsv_per_path = os.path.join(args.outdir, "coverage_per_interval.tsv")

    print(f"Parsing BED file: {args.bed}", file=sys.stderr)
    intervals = parse_bed(args.bed, min_length=args.min_length)
    if not intervals:
        sys.exit("Error: no valid intervals found in BED file.")
    print(f"  Found {len(intervals)} intervals.", file=sys.stderr)

    print(f"Extracting coverage from BAM: {args.bam}", file=sys.stderr)
    matrix = compute_matrix(args.bam, intervals, n_bins=args.bins, min_mapq=args.mapq)

    export_tsv(matrix, args.bins, tsv_path, tsv_per_path)

    print("Rendering plot...", file=sys.stderr)
    plot_coverage(matrix, args.bins, plot_path, args.title)


if __name__ == "__main__":
    main()