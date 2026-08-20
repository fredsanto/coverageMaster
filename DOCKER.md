# coverageMaster_hg38 (Docker)

Docker packaging for `coverageMaster_hg38.py` / `libCoverageMaster_hg38.py` — the
actual CNV-calling script used in production (see `SeqONE/README.md` in the
Medigenome pipeline repo) — not a rewrite. Same code, same output, just
containerized so it runs without the bare-metal CentOS 7 / Python 3.6 /
glibc-pinned-dependency setup the host install requires.

## What's in the image

- `coverageMaster_hg38.py`, `libCoverageMaster_hg38.py`, `HMM_CM.py`,
  `ReadCount.py`, `DGVExplorer.py`, `log.py` — copied in as-is.
- `REF/REFSEQ_hg38_c_HGMD3.gz` + `REF/hg38.exons.merged.bed` — the hg38 gene
  annotation + exon reference the script loads at import time (~11MB,
  baked in rather than mounted since it's a fixed dependency, not per-run data).
- Dependencies pinned to the versions validated against this script on the
  production host: `numpy==1.16.2 sympy==1.0 PyWavelets==1.0.3
  matplotlib==2.2.3 scipy==1.2.1 more-itertools==8.8.0 pandas==1.0.0`, on
  Python 3.7 (newer matplotlib drops the top-level `pylab` module this
  script imports; newer numpy/scipy can shift floating-point behavior in
  the HMM/wavelet calls — don't bump these without re-validating output).

## Build

```bash
docker build -t coveragemaster-hg38 .
```

## Run

Mount whatever directory holds your `.cov` / `.report.txt` / control /
reference files at `/data`, then pass container-side `/data/...` paths as
you would on the bare-metal command line:

```bash
docker run --rm \
    -v /path/to/your/data:/data \
    -w /data \
    coveragemaster-hg38 \
    /data/sample.cov /data/sample.report.txt GENE_OR_GENELIST \
    -c /data/control.cov -r /data/total_ref_m_std \
    -o /data/out_prefix -g /data/CGD_20_25.inh -n 4 -f
```

- `-w /data` matters if you're using a multi-sample control **list** file
  (`-c controls.txt`, one `.cov` path per line) with paths relative to your
  data root — those resolve relative to the container's working directory.
- Output files (`.CMreport`, `.CMcalls`, `.CM.log`, `.CMpositives.pdf`) land
  back on your host under `/path/to/your/data`, since `/data` is a live
  mount.

## Validated 2026-08-20

- Real production case (DME-8598, gene ADK, real `control.DM8xxx_hg38`
  control + `total_ref_m_std` reference): identical result to the
  bare-metal run (`ADK Call: []`).
- Today's error-handling fixes (clear `FATAL SETUP ERROR` messages, exit
  code 1 on bad input) confirmed working inside the container.

## Note

There is a separate, independent Python-package rewrite of coverageMaster
(different code, different CLI, `coverageMaster_v2/` in this repo) — that
one is *not* what this Dockerfile builds. If you're looking for the exact
script the SeqONE pipeline runs today, this is it.
