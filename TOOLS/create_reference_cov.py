import sys
import re
import os
from collections import defaultdict
from glob import glob

MAX_SAMPLES = 50  # cap on how many samples go into the reference


def extract_tot_reads(stats):
    try:
        tot_targ_reads = int(re.findall("([0-9]*).+[0/9].mapped", stats)[0])
        tot_dup_reads = int(re.findall("([0-9]*) \+ 0 *duplicates", stats)[0])  # correct for dups
        return (tot_targ_reads - tot_dup_reads) / 100e6
    except Exception:
        return 0


def dump(d, fout):
    keys = sorted(list(d.keys()), key=lambda k: int(k.split()[-1]))
    for k in keys:
        fout.write("%s\t%s\n" % (k, '\t'.join(d[k])))


def sample_name(path, suffix):
    return os.path.basename(path)[:-len(suffix)]


def index_chromosomes(f):
    """One linear pass recording the byte offset where each chromosome's block
    begins, so later per-chromosome reads can seek directly to it instead of
    rescanning from byte 0 every time (this was the O(n_chromosomes * filesize)
    bottleneck - 25 chromosomes meant every sample file got read ~25x over)."""
    f.seek(0)
    offsets = {}
    prev_chr = None
    pos = f.tell()
    while True:
        line = f.readline()
        if not line:
            break
        try:
            tok = line.split(None, 1)[0].decode()
            chrom = "chr" + tok.replace("chr", "")
        except Exception:
            pos = f.tell()
            continue
        if chrom != prev_chr:
            offsets[chrom] = pos
            prev_chr = chrom
        pos = f.tell()
    return offsets


def pair_samples(cov_dir):
    """Match .cov files to their .report.txt by sample name (not glob order),
    so a mismatched/missing file can't silently pair a sample with the wrong stats."""
    cov_files = {}
    for f in glob(os.path.join(cov_dir, '*.cov')):
        suffix = '.bam.cov' if f.endswith('.bam.cov') else '.cov'
        cov_files[sample_name(f, suffix)] = f

    report_files = {}
    for f in glob(os.path.join(cov_dir, '*.report.txt')):
        suffix = '.bam.report.txt' if f.endswith('.bam.report.txt') else '.report.txt'
        report_files[sample_name(f, suffix)] = f

    paired = []
    for name in sorted(cov_files):
        if name in report_files:
            paired.append((name, cov_files[name], report_files[name]))
        else:
            print("SKIPPING sample '%s': .cov found but no matching .report.txt (%s)" % (name, cov_files[name]))
    for name in sorted(set(report_files) - set(cov_files)):
        print("SKIPPING sample '%s': .report.txt found but no matching .cov (%s)" % (name, report_files[name]))

    if len(paired) > MAX_SAMPLES:
        print("WARNING: %d paired samples found, capping to the first %d (MAX_SAMPLES). "
              "Increase MAX_SAMPLES in this script if you need more." % (len(paired), MAX_SAMPLES))
        paired = paired[:MAX_SAMPLES]

    return paired


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: python create_reference_cov.py <COV_folder> [out_dir]\n")
        sys.exit(1)

    cov_dir = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else cov_dir

    if not os.path.isdir(cov_dir):
        sys.stderr.write("FATAL: input COV folder not found or not a directory: %s\n" % cov_dir)
        sys.exit(1)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        print("Created output directory: %s" % out_dir)
    elif not os.path.isdir(out_dir):
        sys.stderr.write("FATAL: output path exists and is not a directory: %s\n" % out_dir)
        sys.exit(1)

    paired = pair_samples(cov_dir)
    if not paired:
        sys.stderr.write("FATAL: no sample had both a .cov and a matching .report.txt in %s\n" % cov_dir)
        sys.exit(1)

    names = [p[0] for p in paired]
    print("Paired %d samples: %s" % (len(paired), ", ".join(names)))

    open_f = [open(p[1], "rb") for p in paired]
    stats = [extract_tot_reads(open(p[2]).read()) for p in paired]

    zero_stat_samples = [names[n] for n, s in enumerate(stats) if s <= 0]
    if zero_stat_samples:
        print("WARNING: %d sample(s) had unreadable/zero read stats and will be skipped entirely: %s"
              % (len(zero_stat_samples), ", ".join(zero_stat_samples)))

    chr_list = ['chr%s' % str(i) for i in list(range(1, 23)) + ['X', 'Y']] + ["MT"]

    print("Indexing %d sample file(s) (one pass each)..." % len(open_f))
    indexes = []
    for name, f in zip(names, open_f):
        idx = index_chromosomes(f)
        indexes.append(idx)
        missing = [c for c in chr_list if c not in idx]
        if missing:
            print("  %s: %d chromosome block(s) found, missing: %s" % (name, len(idx), ", ".join(missing)))

    try:
        for CHR in chr_list:
            fout = open(os.path.join(out_dir, 'total_%s.res' % CHR), 'w')
            d = defaultdict(list)
            n_used = 0
            for n, file in enumerate(open_f):
                if stats[n] > 0 and CHR in indexes[n]:
                    n_used += 1
                    file.seek(indexes[n][CHR])
                    bad_lines = 0
                    while True:
                        l = file.readline().decode()
                        if not l:
                            break
                        try:
                            chr, pos, cov = l.split()
                            chr = "chr" + chr.replace("chr", "")
                        except Exception:
                            bad_lines += 1
                            continue
                        if CHR != chr:
                            break  # reached the next chromosome's block
                        cov_n = float(cov) / stats[n]  # normalized
                        d[chr + '\t' + pos].append(str(cov_n))
                    if bad_lines:
                        print("  %s: skipped %d malformed line(s) in %s" % (CHR, bad_lines, names[n]))
            print("%s: aggregated %d position(s) from %d sample(s)" % (CHR, len(d), n_used))
            dump(d, fout)
            fout.close()
    finally:
        for f in open_f:
            f.close()


if __name__ == "__main__":
    main()
