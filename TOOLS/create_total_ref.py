from numpy import mean, std
import sys
import os

if len(sys.argv) != 2:
    sys.stderr.write("Usage: python create_total_ref.py <total_ref_file> > total_ref_m_std\n")
    sys.exit(1)

in_path = sys.argv[1]
if not os.path.exists(in_path):
    sys.stderr.write("FATAL: input file not found: %s\n" % in_path)
    sys.exit(1)

n_lines = 0
n_skipped = 0
with open(in_path) as tot_n:
    for lineno, l in enumerate(tot_n, 1):
        ll = l.strip().split()
        if len(ll) < 3:
            sys.stderr.write("WARNING: line %d has no sample values, skipping: %r\n" % (lineno, l.rstrip()))
            n_skipped += 1
            continue
        chr, pos = ll[:2]
        values = ll[2:]  # all per-sample values (was ll[3:], which silently dropped the first sample)
        try:
            fvalues = list(map(float, values))
        except ValueError:
            sys.stderr.write("WARNING: line %d has non-numeric value(s), skipping: %r\n" % (lineno, l.rstrip()))
            n_skipped += 1
            continue
        mn = mean(fvalues)
        st = std(fvalues)
        print("%s\t%s\t%f\t%f" % (chr, pos, mn, st))
        n_lines += 1

sys.stderr.write("Done: %d position(s) written, %d line(s) skipped\n" % (n_lines, n_skipped))
