"""
coveragemaster-compute-stats – Step 2 of reference construction.

Reads the concatenated total_ref file (produced by coveragemaster-build-ref
and chromosome concatenation) and emits one line per position:

    chr  pos  mean  std

Redirect to total_ref_m_std, which is the file passed as -r to coveragemaster.

Usage:
    coveragemaster-compute-stats total_ref > total_ref_m_std
    # or
    coveragemaster-compute-stats total_ref -o total_ref_m_std
"""
import argparse
import sys

import numpy as np


def compute_stats(in_path: str, out_path: str = None):
    """Collapse per-sample coverage columns to mean ± std per position.

    Fixes a bug in the original create_total_ref.py that skipped the first
    sample column (used ll[3:] instead of the correct ll[2:]).
    """
    out = open(out_path, "w") if out_path else sys.stdout
    try:
        with open(in_path) as fh:
            for line in fh:
                cols = line.strip().split()
                if len(cols) < 3:
                    continue
                chr_, pos = cols[0], cols[1]
                # Coverage values begin at index 2 (chr=0, pos=1, cov_s1=2, ...)
                values = list(map(float, cols[2:]))
                mn = float(np.mean(values))
                st = float(np.std(values))
                out.write(f"{chr_}\t{pos}\t{mn:.6f}\t{st:.6f}\n")
    finally:
        if out_path:
            out.close()


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="coveragemaster-compute-stats",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "total_ref",
        help="Concatenated per-chromosome file (output of coveragemaster-build-ref)",
    )
    parser.add_argument(
        "-o", "--out", default=None, metavar="FILE",
        help="Write output to FILE instead of stdout",
    )
    args = parser.parse_args(argv)
    compute_stats(args.total_ref, args.out)


if __name__ == "__main__":
    main()
