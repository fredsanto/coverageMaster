"""
coveragemaster-build-ref – Step 1 of reference construction.

Reads a folder of .cov + .report.txt files, normalises each sample's
coverage by (mapped_reads / 100 M), and writes one file per chromosome:

    <out_dir>/total_<CHR>.res   columns: chr  pos  cov_s1  cov_s2 ...

After this script finishes, concatenate the chromosomes and compute
mean/std with coveragemaster-compute-stats:

    for i in $(seq 1 22) X Y; do
        cat <out_dir>/total_chr$i.res >> total_ref
    done
    coveragemaster-compute-stats total_ref > total_ref_m_std
"""
import argparse
import os
import re
import sys
from collections import defaultdict
from glob import glob


_DEFAULT_CHRS = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]


def _extract_tot_reads(stats_text: str) -> float:
    try:
        tot = int(re.findall(r"([0-9]*) \+ 0 *mapped", stats_text)[0])
        dup = int(re.findall(r"([0-9]*) \+ 0 *duplicates", stats_text)[0])
        return (tot - dup) / 100e6
    except (IndexError, ValueError):
        return 0.0


def _dump(d, fout):
    keys = sorted(d.keys(), key=lambda k: int(k.split()[-1]))
    for k in keys:
        fout.write(f"{k}\t{chr(9).join(d[k])}\n")


def build_ref(cov_dir: str, out_dir: str, chromosomes: list, max_samples: int = 50):
    """Build per-chromosome reference files from a folder of .cov samples.

    Parameters
    ----------
    cov_dir:     Folder with .cov and matching .report.txt files.
    out_dir:     Destination for total_<CHR>.res output files.
    chromosomes: List of chromosome names to process.
    max_samples: Cap on the number of samples (default: 50).
    """
    os.makedirs(out_dir, exist_ok=True)

    cov_files  = sorted(glob(os.path.join(cov_dir, "*.cov")))[:max_samples]
    stat_files = sorted(glob(os.path.join(cov_dir, "*.report.txt")))[:max_samples]

    if not cov_files:
        sys.exit(f"No .cov files found in '{cov_dir}'")

    if len(cov_files) != len(stat_files):
        print(
            f"Warning: {len(cov_files)} .cov files but {len(stat_files)} .report.txt – "
            "pairing by sort order.",
            file=sys.stderr,
        )

    open_f = [open(f, "rb") for f in cov_files]
    stats  = [_extract_tot_reads(open(f).read()) for f in stat_files]

    for chrom in chromosomes:
        out_path = os.path.join(out_dir, f"total_{chrom}.res")
        print(f"  {chrom} → {out_path}", file=sys.stderr)
        d = defaultdict(list)

        for n, fh in enumerate(open_f):
            if stats[n] <= 0:
                continue
            fh.seek(0)
            inside = False
            while True:
                line = fh.readline().decode(errors="replace")
                if not line:
                    break
                parts = line.split()
                if len(parts) < 3:
                    continue
                chr_, pos, cov = parts[0], parts[1], parts[2]
                if chr_ == chrom:
                    cov_n = float(cov) / stats[n]
                    if cov_n <= 10000:
                        d[f"{chr_}\t{pos}"].append(f"{cov_n:.6f}")
                    inside = True
                elif inside:
                    fh.seek(-len(line.encode()), 1)
                    break

        with open(out_path, "w") as fout:
            _dump(d, fout)

    for fh in open_f:
        fh.close()

    print("Done.", file=sys.stderr)


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="coveragemaster-build-ref",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("cov_dir", help="Folder containing .cov and .report.txt files")
    parser.add_argument("out_dir", help="Output folder for total_<CHR>.res files")
    parser.add_argument(
        "--chr",
        dest="chromosomes",
        default=",".join(_DEFAULT_CHRS),
        metavar="LIST",
        help=(
            "Comma-separated chromosome names to process "
            f"(default: {','.join(_DEFAULT_CHRS)})"
        ),
    )
    parser.add_argument(
        "--max-samples", type=int, default=50, metavar="N",
        help="Maximum number of samples to use (default: 50)",
    )
    parser.add_argument(
        "--include-mt", action="store_true", default=False,
        help="Also process mitochondrial chromosome (MT/chrM)",
    )
    args = parser.parse_args(argv)

    chromosomes = [c.strip() for c in args.chromosomes.split(",") if c.strip()]
    if args.include_mt:
        chromosomes += ["MT", "chrM"]

    print(f"Processing {len(chromosomes)} chromosome(s) from '{args.cov_dir}'", file=sys.stderr)
    build_ref(args.cov_dir, args.out_dir, chromosomes, args.max_samples)


if __name__ == "__main__":
    main()
