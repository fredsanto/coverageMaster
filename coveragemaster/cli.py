"""
coverageMaster – main entry point.

Standard usage:
    coveragemaster sample.cov sample.report.txt GENE \\
        -c control.cov -r total_ref_m_std -o sample.out
"""
import argparse
import os
import sys
import time
from functools import partial
from multiprocessing import Manager, Pool

import matplotlib
matplotlib.use("PDF")
from more_itertools import locate
from numpy import (
    array, concatenate, correlate, diff, isnan, isinf,
    median, ones, where, zeros,
)

from coveragemaster import __version__
from coveragemaster.dgv import DGVExplorer
from coveragemaster.hmm import HMM_long
from coveragemaster.lib import (
    calculate_exon,
    convertHMM,
    extract_tot_reads,
    gene_reference,
    load_references,
    plotter,
    wig_writer,
)
from coveragemaster.log import logreport
from coveragemaster.readcount import Regions

# ---------------------------------------------------------------------------
# Module-level worker globals (populated in main() before Pool creation)
# ---------------------------------------------------------------------------
_cr      = None   # Regions: test sample .cov
_cref    = None   # Regions: reference .cov (mean + std)
_ccont   = None   # Regions: current control .cov
_stats   = None   # float: test sample normalisation factor
_cstats  = None   # float: control normalisation factor
_normX   = 1.0    # chrX normalisation factor (test)
_normCX  = 1.0    # chrX normalisation factor (control)
_LOGFILE = None

_lev            = 5
_minlev         = 0
_wid            = 1.0
_offset         = 0.0
_FORCE          = False
_output_px      = ""
_callThreshold  = 0.0
_write_wig      = False


# ---------------------------------------------------------------------------
# Signal processor
# ---------------------------------------------------------------------------

def _signal_processor(gene, cr, cref, stats, signal_buffer=None, store=False):
    """Extract normalised coverage for *gene* from *cr* against *cref*."""
    logfile = _LOGFILE

    start_i = int(gene["end"]) - int(gene["start"])
    if not (start_i > 0) or gene["chr"] == "chrM":
        if signal_buffer and gene.get("gene") in signal_buffer:
            sb = signal_buffer[gene["gene"]]
            return sb["plot_region_n"], sb["ref_exon_avg"], sb["ref_exon_min"], sb["enlight"], sb["data_n"]
        return None, None, None, None, None

    if signal_buffer and gene.get("gene") in signal_buffer:
        sb = signal_buffer[gene["gene"]]
        return sb["plot_region_n"], sb["ref_exon_avg"], sb["ref_exon_min"], sb["enlight"], sb["data_n"]

    qregion = gene.get("query", gene)
    if "gene" not in qregion:
        qregion = dict(qregion)
        qregion["gene"] = f"{gene['chr']}:{gene['start']}-{gene['end']}"

    logreport(
        f"Processing {'Gene: ' + gene.get('gene', '') + ' - ' if 'gene' in gene else ''}"
        f"Region {qregion['chr']}:{qregion['start']}-{qregion['end']}",
        logfile=logfile,
    )

    size = int(qregion["end"]) - int(qregion["start"])
    step = 1 if size < 1e6 else 5

    pos_region    = []
    ref_exon_min  = []
    ref_exon_avg  = []

    gen_cr = cr.focus_single_prl(
        (qregion["chr"], int(qregion["start"]), int(qregion["end"])), step=1
    )

    for r in cref.focus_single_prl(
        (qregion["chr"], int(qregion["start"]), int(qregion["end"])), step=step
    ):
        rr = r.strip().split()
        chr_, pos = rr[0], rr[1]
        cov_cols = rr[2:]

        if len(cov_cols) > 3:
            vals = list(map(float, cov_cols))
            mcov = float(array(vals).mean())
            stdv = float(array(vals).std())
        elif len(cov_cols) > 1:
            mcov = float(cov_cols[0])
            stdv = float(cov_cols[1])
        else:
            mcov = float(cov_cols[0]) / stats
            stdv = 1.0

        LOOP = FREEZE = FREEZE_B = False
        count = 0
        gpos = pos
        cov = 0
        while not LOOP and not FREEZE:
            count += 1
            if count > 10000:
                logreport(f"Error: Loophole at {chr_} {pos}", logfile=logfile)
                break
            try:
                if not FREEZE_B:
                    _, gpos, cov = next(gen_cr).strip().split()
                if int(gpos) == int(pos) + 1:
                    FREEZE = True
                if int(gpos) > int(pos) + 1:
                    FREEZE_B = True
                    break
            except StopIteration:
                gpos, cov = pos, 0
                break
            LOOP = (pos == gpos)

        if int(gpos) < int(pos):
            FREEZE = FREEZE_B = False
        if pos == gpos:
            FREEZE = FREEZE_B = False
            if stdv == 0:
                stdv = 1.0
            ref_exon_min.append(stdv)
            ref_exon_avg.append(mcov)
            enlight_val = (
                -1.5 + calculate_exon(gene, pos)
                if int(gene["start"]) <= int(pos) <= int(gene["end"])
                else 0
            )
            pos_region.append((pos, float(cov) / stats, enlight_val))

    if not pos_region:
        return None, None, None, None, None

    plot_region  = array([p[1] for p in pos_region])
    enlight      = [p[2] for p in pos_region]
    ref_arr      = array(ref_exon_avg)

    plot_region_n = plot_region / ref_arr
    plot_region_n = array([
        1.0 if (isnan(v) or isinf(v)) else v for v in plot_region_n
    ])
    data_n = [(p[0], p[1], plot_region_n[i]) for i, p in enumerate(pos_region)]

    if store and signal_buffer is not None:
        signal_buffer[gene["gene"]] = {
            "plot_region_n": plot_region_n,
            "ref_exon_avg":  ref_arr,
            "ref_exon_min":  ref_exon_min,
            "enlight":       enlight,
            "data_n":        data_n,
            "infidx":        [],
            "noninfidx":     [],
        }

    return plot_region_n, ref_arr, ref_exon_min, enlight, data_n


# ---------------------------------------------------------------------------
# Per-gene CNV worker
# ---------------------------------------------------------------------------

def _process_coverage(terminal, gene, signal_buffer):
    time.sleep(0.1)

    try:
        _cr.f   = open(_cr.f.name)
        _cref.f = open(_cref.f.name)
        _ccont.f = open(_ccont.f.name)
    except Exception:
        raise RuntimeError("Invalid control/reference files")

    signal, ref_exon_avg, ref_exon_min, enlight, data_n = _signal_processor(
        gene, _cr, _cref, _stats, signal_buffer, store=True
    )
    if signal is None:
        return "NO COVERAGE"

    csignal, _, _, _, _ = _signal_processor(gene, _ccont, _cref, _cstats)

    if gene["chr"] == "chrX":
        signal  = signal / _normX
        csignal = csignal / _normCX + 1e-4
    else:
        signal = signal + _offset

    if sum(signal) == 0 and sum(csignal) == 0:
        return "NO COVERAGE"
    if sum(csignal) == 0:
        logreport(f"Warning: no coverage in control for {gene}", logfile=_LOGFILE)
        return "NO COVERAGE IN CONTROL"

    if len(signal) < len(csignal):
        logreport(f"Warning: signal < control for {gene} – trimming", logfile=_LOGFILE)
        csignal = csignal[: len(signal)]
    if len(signal) > len(csignal):
        logreport(f"Warning: signal > control for {gene} – padding", logfile=_LOGFILE)
        csignal = array(csignal.tolist() + [1.0] * (len(signal) - len(csignal)))

    if len(signal) != len(csignal):
        return None

    sd       = array(ref_exon_min) / array(ref_exon_avg)
    genethere = array(enlight) != 0

    stdM = 1 + _wid * sd
    stdm = 1 - _wid * sd
    unormsignal  = signal  - 1
    unormcsignal = csignal - 1

    informative = (abs(unormsignal) > 0.4) * (
        (unormcsignal < _wid * sd) * (unormsignal > _wid * sd) * (csignal < signal)
        + (unormsignal < -_wid * sd) * (csignal > signal)
    )
    infidx_org = where(informative)[0]

    ratiofull = signal / (0.001 + csignal)

    _t = signal_buffer[gene["gene"]]
    if len(_t["infidx"]) == 0:
        _t["infidx"] = infidx_org
    elif len(set(infidx_org).intersection(_t["infidx"])):
        _t["infidx"] = infidx_org

    win = 50
    _tmp = [
        range(max(0, i - win), min(len(signal), i + win))
        for i in _t["infidx"]
        if win < i < (len(signal) - win)
        and not (unormcsignal[i] < unormsignal[i] < -_wid * sd[i])
    ]
    infidx = list({idx for rng in _tmp for idx in rng})
    signal_buffer[gene["gene"]] = _t

    ratio = ones(len(signal))
    orgmedian = median(ratiofull)
    ratio[infidx] = ratiofull[infidx]

    try:
        approx, alev = HMM_long(
            ratio, sd, gene["gene"],
            booster=1, lev=_lev, LOGFILE=_LOGFILE,
            mask=genethere, minlev=_minlev, orgmedian=orgmedian,
        )
    except Exception:
        logreport(f"HMM problem: {gene['gene']}", logfile=_LOGFILE)
        return None

    detDEL = sum(approx == 0.5)
    detDUP = sum(approx == 1.5)

    if detDEL:
        del_ok = sum((approx == 0.5) * informative) / sum(approx == 0.5) > 0.1
        delreg_approx = (approx == 0.5) if del_ok else (approx == 0)
    else:
        delreg_approx = 0

    if detDUP:
        dup_ok = sum((approx == 1.5) * informative) / sum(approx == 1.5) > 0.1
        dupreg_approx = (approx == 1.5) if dup_ok else (approx == 0)
    else:
        dupreg_approx = 0

    approx = (approx - 1) * delreg_approx + (approx - 1) * dupreg_approx + 1
    mask   = (approx != 1) * 1
    Q = []

    if sum(mask):
        ovr = 15
        dmask    = diff(mask)
        poscalls = where(dmask != 0)[0].tolist()
        if dmask[poscalls[0]] != 1:
            poscalls = [0] + poscalls
        if dmask[poscalls[-1]] != -1:
            poscalls = poscalls + [len(mask) - 1]
        ncalls = len(poscalls)

        for i in range(0, ncalls, 2):
            mask_tmp = zeros(len(mask))
            val  = mask[int((poscalls[i] + poscalls[i + 1]) / 2)]
            ovrl = ovr if poscalls[i] > ovr else 0
            ovrr = ovr if poscalls[i + 1] + ovr < len(mask) else 0
            seg  = concatenate((
                val * ones(ovrl),
                mask[poscalls[i]: poscalls[i + 1]],
                ones(ovrr) * val,
            ))
            mask[poscalls[i] - ovrl: poscalls[i + 1] + ovrr] = seg
            mask_tmp[poscalls[i] - ovrl: poscalls[i + 1] + ovrr] = seg
            base = mask_tmp * genethere
            if sum(base):
                psignal = (
                    0.5 * (abs(signal - 1) >= 0.5)
                    + abs(signal - 1) * (abs(signal - 1) <= 0.5)
                )
                Q.append(float(
                    correlate(psignal, 0.5 * base) / correlate(0.5 * base, 0.5 * base)
                ))
            else:
                Q.append(0.0)
            logreport(f"Call quality: {gene['gene']} – {Q[-1]:.4f}", logfile=_LOGFILE)

    has_call = sum(mask) > 0 and sum(signal * mask * genethere) != 0
    if _FORCE or _write_wig or (has_call and Q and max(Q) > _callThreshold):
        if _write_wig:
            wig = wig_writer(gene["chr"], data_n)
            with open(f"{_output_px}_{gene['gene']}.wig", "w") as fw:
                fw.write("\n".join(wig))
        call = convertHMM(approx, data_n)
        LONG = "LONG" if any(int(c[-1]) > 1e5 for c in call) else ""
        logreport(f"{gene['gene']} Call: {call} {LONG}", logfile=_LOGFILE)
        return [
            gene, signal, csignal, enlight, stdm, stdM,
            approx, ref_exon_avg / ref_exon_avg.sum(), ratio, call, Q,
        ]
    else:
        signal_buffer.pop(gene["gene"], None)
        logreport(f"Removing {gene['gene']}", logfile=_LOGFILE)
        return None


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------

def _build_parser():
    parser = argparse.ArgumentParser(
        prog="coveragemaster",
        description="coverageMaster – nucleotide-resolution CNV caller (Wavelet + HMM).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single gene, one control (hg38 default)
  coveragemaster sample.cov sample.report.txt PGM1 \\
      -c control.cov -r total_ref_m_std -o out/sample

  # Gene list, multiple controls
  coveragemaster sample.cov sample.report.txt genes.txt \\
      -c controls.txt -r total_ref_m_std -o out/sample

  # BED regions, hg19
  coveragemaster sample.cov sample.report.txt regions.bed -b \\
      -c control.cov -r total_ref_m_std -o out/sample --genome hg19

  # Custom reference files
  coveragemaster sample.cov sample.report.txt PGM1 \\
      -c control.cov -r total_ref_m_std -o out/sample \\
      --refseq /data/REF/REFSEQ_hg38_HGMD3 \\
      --exon-bed /data/REF/hg38.exons.merged.bed
""",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    # Positional arguments
    parser.add_argument("cov_file",   help="samtools depth output (.cov)")
    parser.add_argument("stats_file", help="samtools flagstat output (.report.txt)")
    parser.add_argument(
        "targets",
        help=(
            "Gene name, comma-separated gene names, gene-list file, "
            "genomic region (chr:start-end), or BED file path (requires -b)"
        ),
    )

    # Required options
    req = parser.add_argument_group("required")
    req.add_argument("-c", "--control", required=True,
                     help="Control .cov file, or text file with one .cov path per line")
    req.add_argument("-r", "--ref", required=True,
                     help="Reference coverage file (mean+std, built with coveragemaster-build-ref)")
    req.add_argument("-o", "--out", dest="output_px", required=True,
                     help="Output file prefix")

    # Reference genome selection
    genome_grp = parser.add_argument_group("reference genome")
    genome_grp.add_argument(
        "--genome", choices=["hg19", "hg38"], default="hg38",
        help="Genome assembly for transcript/exon annotation (default: hg38)",
    )
    genome_grp.add_argument(
        "--refseq", default=None, metavar="FILE",
        help="Custom RefSeq transcript table — overrides --genome",
    )
    genome_grp.add_argument(
        "--exon-bed", default=None, metavar="FILE",
        help="Custom merged exon BED file — overrides --genome",
    )

    # Algorithm tuning
    algo = parser.add_argument_group("algorithm")
    algo.add_argument("-e", "--exons",    type=int,   default=5,   metavar="N",
                      help="Extra flanking exons on each side of the query window (default: 5)")
    algo.add_argument("-d", "--width",    dest="wid", type=float, default=1.0,
                      help="Significance band width in units of std (default: 1)")
    algo.add_argument("-l", "--level",    dest="lev", type=int,   default=5,
                      help="Maximum wavelet compression level (default: 5)")
    algo.add_argument("-m", "--minlevel", dest="minlev", type=int, default=0,
                      help="Minimum wavelet level for zooming (default: 0)")
    algo.add_argument("-x", "--offset",   type=float, default=0.0,
                      help="Offset added to the signal (default: 0)")
    algo.add_argument("-T", "--threshold", dest="cT", type=float, default=0.0,
                      help="Minimum call quality score to report (default: 0)")

    # Output options
    out_grp = parser.add_argument_group("output")
    out_grp.add_argument("-g", "--cgd", default="",
                         help="CGD inheritance file (gene → mode), adds annotation to calls")
    out_grp.add_argument("-f", "--force", action="store_true", default=False,
                         help="Force PDF output even when no CNV is called")
    out_grp.add_argument("-w", "--wig", action="store_true", default=False,
                         help="Write per-gene WIG files for UCSC Genome Browser")
    out_grp.add_argument("-k", "--dgv", default="",
                         help="DGV tabix file for population frequency annotation")

    # Execution
    exec_grp = parser.add_argument_group("execution")
    exec_grp.add_argument("-n", "--nproc", type=int, default=5,
                           help="Number of parallel worker processes (default: 5)")
    exec_grp.add_argument("-b", "--bed", action="store_true", default=False,
                           help="Treat <targets> as a BED file of regions")

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    global _cr, _cref, _ccont, _stats, _cstats
    global _normX, _normCX, _LOGFILE
    global _lev, _minlev, _wid, _offset, _FORCE
    global _output_px, _callThreshold, _write_wig

    parser = _build_parser()
    args   = parser.parse_args(argv)

    # Validate input
    if not os.path.exists(args.cov_file):
        parser.error(f"Input .cov file not found: {args.cov_file}")
    if not os.path.exists(args.stats_file):
        parser.error(f"Stats file not found: {args.stats_file}")

    # Set worker globals
    _output_px     = args.output_px
    _LOGFILE       = open(_output_px + ".CM.log", "a")
    _lev           = args.lev
    _minlev        = args.minlev
    _wid           = args.wid
    _offset        = args.offset
    _FORCE         = args.force
    _callThreshold = args.cT
    _write_wig     = args.wig

    logreport(f"{'–' * 80}\ncoverageMaster {__version__} warming up", logfile=_LOGFILE)
    logreport(f"Args: {vars(args)}", logfile=_LOGFILE)

    # Load reference genome annotations
    try:
        reference, exon_reference = load_references(
            genome=args.genome,
            refseq=args.refseq,
            exon_bed=args.exon_bed,
        )
        logreport(f"Genome: {args.genome}", logfile=_LOGFILE)
    except FileNotFoundError as exc:
        parser.error(str(exc))

    # Parse stats
    try:
        _stats = extract_tot_reads(open(args.stats_file).read())
    except Exception:
        parser.error(f"Cannot parse stats file: {args.stats_file}")

    # Load CGD
    cgd = {}
    if args.cgd:
        with open(args.cgd) as fh:
            for line in fh:
                try:
                    gname, inh = line.strip().split()
                    cgd[gname] = inh
                except ValueError:
                    pass

    # Build query region list
    nexons = args.exons
    XIST   = gene_reference(["XIST"], reference, exon_reference, _LOGFILE, nexons=50)
    qregions = []

    if args.bed:
        with open(args.targets) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                qregions.append({"chr": parts[0], "start": parts[1], "end": parts[2]})
        qregions += XIST

    elif ":" in args.targets:
        chr_, coords = args.targets.split(":")
        start, end = coords.split("-")
        qregions = [{"chr": chr_, "start": start, "end": end}] + XIST

    else:
        # Gene name(s) or gene-list file
        try:
            with open(args.targets) as fh:
                glist = fh.read().strip().split()
            if len(glist) < 2:
                glist = open(args.targets).read().strip().split("\n")[0].split()
        except (OSError, IOError):
            glist = args.targets.split(",")

        if ":" in glist[0]:
            for r in glist:
                chr_, coords = r.split(":")
                start, end = coords.split("-")
                qregions.append({"chr": chr_, "start": start, "end": end})
            qregions += XIST
        else:
            cache = args.targets + "_coderef" if os.path.isfile(args.targets) else None
            qregions = (
                gene_reference(
                    glist, reference, exon_reference, _LOGFILE,
                    nexons=nexons, cache_path=cache,
                )
                + XIST
            )

    if not qregions:
        parser.error("No valid query regions resolved. Check gene names or region format.")

    # DGV
    dgv_xplr = None
    if args.dgv:
        filt_columns = [
            "#chrom", "chromStart", "chromEnd", "name", "varType",
            "supportingVariants", "sampleSize", "observedGains", "observedLosses",
        ]
        dgv_xplr = DGVExplorer(
            file_path=args.dgv,
            filt_cols=filt_columns,
            chr_idx_path=args.dgv + "_chr_indices.npy",
        )

    # Index files
    _cr   = Regions(args.cov_file)
    logreport("Query index created", logfile=_LOGFILE)
    _cref = Regions(args.ref)
    logreport("Reference index created", logfile=_LOGFILE)

    # chrX normalisation (test sample)
    XIST_region = qregions.pop()
    chrX_signal, _, _, _, _ = _signal_processor(XIST_region, _cr, _cref, _stats)
    _normX = float(median(chrX_signal)) if chrX_signal is not None else 1.0

    # Parse control list
    if args.control.endswith(".cov"):
        cc = [args.control]
    else:
        cc = open(args.control).read().strip().split("\n")
    cc = [c for c in cc if os.path.exists(c)]
    if not cc:
        parser.error("No valid control .cov files found.")

    # Main loop over controls
    qregions_failed     = []
    qregions_failed_in_control = []
    repository          = []
    signal_buffer       = Manager().dict()

    for n, c in enumerate(cc):
        if c.startswith("#") or os.path.basename(c) == os.path.basename(args.cov_file):
            continue

        terminal    = (n == len(cc) - 1)
        cstats_path = c.replace(".cov", "") + ".report.txt"
        _ccont  = Regions(c)
        _cstats = extract_tot_reads(open(cstats_path).read())
        logreport(f"Control: {c}", logfile=_LOGFILE)

        ctrl_chrX, _, _, _, _ = _signal_processor(XIST_region, _ccont, _cref, _cstats)
        _normCX = float(median(ctrl_chrX)) if ctrl_chrX is not None else 1.0

        pfunc = partial(_process_coverage, terminal, signal_buffer=signal_buffer)
        with Pool(processes=args.nproc) as pool:
            repository = pool.map(pfunc, qregions)

        if repository:
            qregions_next = [b[0] for b in repository if b and not isinstance(b, str)]
            qregions_failed_in_control = [
                qregions[q]
                for q in locate(repository, lambda x: x == "NO COVERAGE IN CONTROL")
            ]
            qregions_failed += [
                qregions[q]
                for q in locate(repository, lambda x: x == "NO COVERAGE")
            ]
            if qregions_failed:
                logreport(f"No coverage: {qregions_failed}", logfile=_LOGFILE)
            qregions   = qregions_next + qregions_failed_in_control
            repository = [
                r for r in repository
                if r and r not in ("NO COVERAGE", "NO COVERAGE IN CONTROL")
            ]
        else:
            break

        logreport(f"N.calls: {len(repository)}", logfile=_LOGFILE)

        if len(qregions) < 400 and not qregions_failed_in_control:
            plotter(repository, _output_px, qregions_failed, cgd, dgv_xplr)

        if _FORCE or args.wig or len(qregions) == 0:
            break

    if repository:
        plotter(repository, _output_px, qregions_failed, cgd, dgv_xplr)

    if qregions_failed_in_control:
        logreport(
            f"WARNING: regions not covered in last control:\n{qregions_failed_in_control}",
            logfile=_LOGFILE,
        )

    logreport("coverageMaster is done", logfile=_LOGFILE)
    _LOGFILE.close()


if __name__ == "__main__":
    main()
