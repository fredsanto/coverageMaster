"""
coverageMaster library
Provides reference loading, gene coordinate resolution, signal helpers, and plotting.
"""
import gzip
import os
import pickle
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from numpy import array, isnan, isinf, median, diff, where
from sympy import Interval, Union

from coveragemaster.log import logreport

# ---------------------------------------------------------------------------
# Reference genome configuration
# ---------------------------------------------------------------------------

_GENOME_FILES = {
    "hg19": {
        "refseq":   "REFSEQ_hg19.chr.complete.txt.gz",
        "exon_bed": "hg19.exons.merged.bed",
    },
    "hg38": {
        "refseq":   "REFSEQ_hg38_HGMD3",
        "exon_bed": "hg38.exons.merged.bed",
    },
}


def _find_ref_file(filename: str) -> str:
    """Locate a bundled reference file.

    Search order:
    1. ``COVERAGEMASTER_REF_DIR`` environment variable.
    2. ``<package_dir>/REF/<filename>``    – pip-installed layout.
    3. ``<package_dir>/../REF/<filename>`` – git-checkout / editable install.
    """
    env_dir = os.environ.get("COVERAGEMASTER_REF_DIR")
    if env_dir:
        candidate = os.path.join(env_dir, filename)
        if os.path.exists(candidate):
            return candidate
        raise FileNotFoundError(
            f"COVERAGEMASTER_REF_DIR='{env_dir}' but '{filename}' not found there."
        )

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    for rel in ("REF", os.path.join("..", "REF")):
        candidate = os.path.normpath(os.path.join(pkg_dir, rel, filename))
        if os.path.exists(candidate):
            return candidate

    raise FileNotFoundError(
        f"Reference file '{filename}' not found.\n"
        f"Set COVERAGEMASTER_REF_DIR, or pass --refseq / --exon-bed explicitly."
    )


def load_references(genome: str = "hg38", refseq: str = None, exon_bed: str = None):
    """Load the transcript reference and exon BED for *genome*.

    Parameters
    ----------
    genome:
        ``'hg19'`` or ``'hg38'`` (default: ``'hg38'``).
    refseq:
        Path to a custom RefSeq transcript table (overrides the *genome* default).
    exon_bed:
        Path to a custom merged exon BED file (overrides the *genome* default).

    Returns
    -------
    reference : list[str]
        Lines of the RefSeq transcript table.
    exon_reference : dict[str, list[tuple]]
        Exon intervals keyed by chromosome string.
    """
    if genome not in _GENOME_FILES:
        raise ValueError(
            f"Unknown genome '{genome}'. Available: {list(_GENOME_FILES)}"
        )

    files = _GENOME_FILES[genome]
    refseq_path   = refseq   if refseq   else _find_ref_file(files["refseq"])
    exon_bed_path = exon_bed if exon_bed else _find_ref_file(files["exon_bed"])

    if refseq_path.endswith(".gz"):
        with gzip.open(refseq_path) as fh:
            reference = fh.read().decode().strip().split("\n")
    else:
        with open(refseq_path) as fh:
            reference = fh.read().strip().split("\n")

    exon_reference = defaultdict(list)
    with open(exon_bed_path) as fh:
        for line in fh.read().strip().split("\n"):
            parts = line.split()
            if len(parts) >= 3:
                exon_reference[parts[0]].append((parts[1], parts[2]))

    return reference, exon_reference


# ---------------------------------------------------------------------------
# Stats helper
# ---------------------------------------------------------------------------

def extract_tot_reads(stats: str) -> float:
    """Parse ``samtools flagstat`` text and return mapped_reads / 100 M."""
    tot = int(re.findall(r"([0-9]*) \+ 0 *mapped", stats)[0])
    return tot / 100e6


# ---------------------------------------------------------------------------
# WIG writer
# ---------------------------------------------------------------------------

def wig_writer(chr_, data):
    """Convert *data* (list of (pos, raw, norm) tuples) to UCSC variableStep WIG."""
    span = 20.0
    wig = ["track type=wiggle_0 name=p_full"]
    for snp in data:
        pos, _, value = snp
        c = chr_ if "MT" in chr_ else "chr" + chr_.replace("chr", "")
        wig.append(f"variableStep chrom={c} span={int(2 * span)}")
        for p in range(int(pos) - int(span), int(pos) + int(span)):
            slope = 1 if p == float(pos) else 0
            wig.append(f"{int(p)}\t{slope * value:.6f}")
    return wig


# ---------------------------------------------------------------------------
# Gene coordinate resolution
# ---------------------------------------------------------------------------

def _union(data):
    """Sympy Union of (start, end) interval pairs."""
    intervals = [Interval(b, e) for b, e in data]
    u = Union(*intervals)
    return [u] if isinstance(u, Interval) else list(u.args)


def exon_boundaries(chr_, start, end, exon_reference, nexons=10):
    """Extend a gene's boundaries by *nexons* flanking exons on each side."""
    offset = nexons
    START = False
    END = False
    for n, (estart, eend) in enumerate(exon_reference[chr_]):
        if int(start) <= int(estart) and not START:
            START = exon_reference[chr_][max(n - offset, 0)][0]
        if int(eend) >= int(end) and not END:
            idx = min(n + offset, len(exon_reference[chr_]) - 1)
            END = exon_reference[chr_][idx][1]
    return {"chr": chr_, "start": START, "end": END}


def gene_reference(genes, reference, exon_reference, LOGFILE, nexons=10, cache_path=None):
    """Resolve gene names to genomic coordinate dicts using the RefSeq table.

    Parameters
    ----------
    genes : list[str]
    reference : list[str]
        Lines from :func:`load_references`.
    exon_reference : dict
        Exon intervals from :func:`load_references`.
    LOGFILE : file-like
    nexons : int
        Extra flanking exons for the query window (default: 10).
    cache_path : str, optional
        Pickle cache path.  Loaded on hit; written on miss (only when nexons ≤ 10).

    Returns
    -------
    list[dict]
    """
    if cache_path and nexons <= 10 and os.path.exists(cache_path):
        cached = pickle.load(open(cache_path, "rb"))
        if cached:
            return cached

    gf = []
    for g in sorted(set(genes)):
        try:
            transcripts = [r for r in reference if g == r.split()[-1]]
            tx_X    = [t for t in transcripts if t.split()[1].replace("chr", "") == "X"]
            tx_Y    = [t for t in transcripts if t.split()[1].replace("chr", "") == "Y"]
            tx_auto = [t for t in transcripts if t not in tx_X + tx_Y]

            for txs in [tx_auto, tx_X, tx_Y]:
                if not txs:
                    continue
                maxend, minstart = 0, 9_999_999_999_999_999
                tmax, exonst, exonen = [], [], []
                for t in txs:
                    tt = t.split()
                    if len(tt[1]) < 6:          # skip alt chromosomes
                        exonst += tt[7].strip(",").split(",")
                        exonen += tt[8].strip(",").split(",")
                        tmax = tt if len(tt) > len(tmax) else tmax
                        maxend   = max(maxend,   int(tt[4]))
                        minstart = min(minstart, int(tt[3]))
                tmax[3] = str(minstart)
                tmax[4] = str(maxend)
                tmax[7] = ",".join(exonst)
                tmax[8] = ",".join(exonen)
                tfinal = "\t".join(tmax)
                if maxend - minstart < 3e6:
                    gf.append(tfinal)
                else:
                    logreport(f"Transcript too long, skipping: {g}", logfile=LOGFILE)
        except Exception:
            logreport(f"{g} not found in reference. Check gene name.", logfile=LOGFILE)

    gs = []
    for g in gf:
        cols = g.split("\t")
        exfirst, exlast = cols[3], cols[4]
        exstart = list(map(int, cols[7].split(","))) if len(cols[7].split(",")) > 1 else [int(cols[7])]
        exend   = list(map(int, cols[8].split(","))) if len(cols[8].split(",")) > 1 else [int(cols[8])]
        newex   = _union(list(zip(exstart, exend)))
        exstart = [str(x.start) for x in newex]
        exend   = [str(x.end)   for x in newex]

        try:
            replace_s = exstart.index([x for x in exstart if float(exfirst) <= float(x)][0]) - 1
        except (IndexError, ValueError):
            replace_s = len(exstart) - 1
        try:
            replace_e = exend.index([x for x in exend if float(exlast) >= float(x)][-1]) + 1
        except (IndexError, ValueError):
            replace_e = 0

        replace_s = max(replace_s, 0)
        replace_e = min(replace_e, len(exend) - 1)
        exstart[replace_s] = exfirst
        exend[replace_e]   = exlast
        expairs = list(zip(exstart[replace_s: replace_e + 1], exend[replace_s: replace_e + 1]))

        gdata = {
            "NM": cols[0], "chr": cols[1], "expairs": expairs,
            "gene": cols[-1], "start": exfirst, "end": exlast, "LONG": "",
        }
        gdata["query"] = exon_boundaries(
            gdata["chr"], gdata["start"], gdata["end"], exon_reference, nexons
        )
        gs.append(gdata)

    if cache_path and nexons <= 10 and gs and not os.path.exists(cache_path):
        try:
            pickle.dump(gs, open(cache_path, "wb"))
        except Exception:
            pass

    return gs


# ---------------------------------------------------------------------------
# Signal helpers
# ---------------------------------------------------------------------------

def calculate_exon(gene, pos):
    """Return ±1 for odd/even exon at *pos*, or 0 if intergenic."""
    try:
        exn = gene["expairs"].index(
            [x for x in gene["expairs"] if int(x[0]) < int(pos) < int(x[1])][0]
        )
        return 1 if exn % 2 else -1
    except (IndexError, ValueError):
        return 0


def convertHMM(approx, data_n):
    """Convert HMM state array to a list of (type, start, end, size) call tuples."""
    if not len(data_n):
        return "NODATAN"
    call = []
    dd = diff([1] + approx.tolist() + [1])
    transitions = where(dd != 0)[0]
    txp = [transitions[i: i + 2] for i in range(0, len(transitions), 2)]
    for n, el in enumerate(txp):
        if len(el) == 1:
            el = (txp[n - 1][1], el[0])
        if dd[el[0]] < dd[el[1]]:
            mut = "DEL"
        elif dd[el[0]] > dd[el[1]]:
            mut = "DUP"
        else:
            mut = "ERROR"
        try:
            size = int(data_n[el[1] - 1][0]) - int(data_n[el[0]][0])
            call.append((mut, data_n[el[0]][0], data_n[el[1] - 1][0], str(size)))
        except Exception:
            print(f"Problem converting HMM segment {el}", file=sys.stderr)
    return call


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _coverage_panel(ax, signal, csignal, approx, ref_exon_avg, stdM, stdm, xlim=None):
    """Render a coverage overview panel onto *ax*.

    Y-axis: 0 (bottom) → 1 (diploid centre) → 2 (top).
    A dashed grey reference line marks y = 1.
    """
    ax.plot(stdM, "r", linewidth=0.8)
    ax.plot(stdm, "r", linewidth=0.8)
    ax.axhline(y=1, color="gray", linestyle="--", linewidth=0.5)
    ax.plot(csignal,     color="#111111", linewidth=0.8, linestyle="--")
    ax.plot(signal,      color="b",       linewidth=1.0)
    ax.plot(approx,      "k",             linewidth=1.5)
    ax.plot(ref_exon_avg, "k",            linewidth=1.5)
    ax.set_ylim(0, 2)
    ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
    ax.set_yticklabels(["0", "0.5", "1", "1.5", "2"], fontsize=7)
    ax.set_xticks([])
    if xlim is not None:
        ax.set_xlim(*xlim)


def plotter(repository, output_px, qregions_failed, cgd=None, dgv_xplr=None):
    """Write the PDF report and CMcalls / CMreport text files."""
    if cgd is None:
        cgd = {}

    report    = open(output_px + ".CMreport", "w")
    calls_out = open(output_px + ".CMcalls",  "w")

    try:
        repository = sorted(repository, key=lambda e: max(e[-1]), reverse=True)
    except Exception:
        pass

    with PdfPages(output_px + ".CMpositives.pdf") as pp:
        for box in repository:
            if not box or box == "NO COVERAGE":
                continue

            gene, signal, csignal, enlight, stdm, stdM, approx, ref_exon_avg, ratio, call, Q = box
            cgd_inh = cgd.get(gene["gene"], "")

            # DGV overlaps
            overlaps = None
            if dgv_xplr is not None:
                overlaps = dgv_xplr.get_overlap(
                    chrom=gene["chr"],
                    alpha_0=int(gene["start"]),
                    beta_0=int(gene["end"]),
                )
                if overlaps is not None:
                    freq = dgv_xplr.get_freq(overlaps)
                    overlaps["Gain freq."] = freq["freq_g"].tolist()
                    overlaps["Loss freq."] = freq["freq_l"].tolist()
                    overlaps["Name."]      = freq["name"].tolist()

            # Text labels
            gene_txt = "{}\t{}\t{}".format(
                gene["gene"], cgd_inh, " ".join(f"{q:.2f}" for q in Q)
            )
            call_txt = ""
            try:
                for c in call:
                    size = int(c[2]) - int(c[1])
                    call_txt += f"{gene['chr']}\t{c[1]}\t{c[2]}\t{c[0]}\t{size}\t"
                calls_out.write(gene_txt + " " + call_txt)
                if overlaps is not None:
                    calls_out.write("Gains")
                    for gain, name in zip(overlaps["Gain freq."].values, overlaps["Name."].values):
                        calls_out.write(f"\t{name}\t{gain:.3f}\t")
                    calls_out.write("Losses")
                    for loss, name in zip(overlaps["Loss freq."].values, overlaps["Name."].values):
                        calls_out.write(f"\t{name}\t{loss:.3f}\t")
                calls_out.write("\n")
            except Exception:
                pass

            fig, axes = plt.subplots(4, 1, figsize=(7, 10))

            # Panel 1 – gene structure
            axes[0].plot(enlight, color="g", linewidth=1.0)
            axes[0].set_ylim(-3, 1)
            axes[0].set_xticks([])
            axes[0].set_yticks([])
            axes[0].set_title(gene_txt.replace("\t", " ") + "\n" + call_txt, fontsize=8)

            # Panel 2 – full coverage overview (fixed scale: 0–2, centre=1)
            _coverage_panel(axes[1], signal, csignal, approx, ref_exon_avg, stdM, stdm)

            # Panel 3 – ratio track
            axes[2].axhline(y=0, color="k", linewidth=0.8)
            axes[2].plot(ratio - 1, "k", linewidth=1.5)
            axes[2].set_ylim(-1, 1)
            axes[2].set_xticks([])
            axes[2].set_yticks([-1, -0.5, 0, 0.5, 1])
            axes[2].set_yticklabels(["-1", "-0.5", "0", "0.5", "1"], fontsize=7)

            # Panel 4 – zoomed call region
            enlight_arr = array(enlight)
            roi = where(((approx - 1) * (enlight_arr != 0)) != 0)[0]
            if len(roi):
                xlim = (int(min(roi)) - 100, int(max(roi)) + 100)
                _coverage_panel(axes[3], signal, csignal, approx, ref_exon_avg, stdM, stdm, xlim=xlim)
            else:
                axes[3].set_visible(False)

            report.write(f"{gene['gene']}\n")
            plt.tight_layout()
            plt.savefig(pp, format="pdf")
            plt.close()

    for q in qregions_failed:
        calls_out.write(f"{q['chr']}\t{q['start']}\t{q['end']}\tNO COVERAGE\n")

    report.close()
    calls_out.close()
