# coverageMaster (Docker)

CNV caller from `samtools depth` coverage data, packaged as a self-contained Docker
image. This replaces the legacy bare-metal install (which required CentOS 7 /
Python 3.6 and pinned numpy/scipy/matplotlib versions to match the OS's glibc) —
the container ships its own Python 3.10 environment, so it runs on any Docker host.

## Quickstart

```bash
# 1. Build the image (from this directory)
docker build -t coveragemaster .

# 2. Run it against your data
./run.sh <data_dir> <ref_dir> \
    /data/sample.cov /data/sample.report.txt GENE_OR_GENELIST \
    -c /data/control.cov -r /data/ref_file \
    -o /data/out_prefix --genome hg38
```

- `<data_dir>` — local folder with your `.cov` / `.report.txt` / control / gene-list files. Mounted read-write at `/data` inside the container; use `/data/...` paths in the arguments.
- `<ref_dir>` — local folder with the RefSeq annotation + exon BED reference files (the `REF/` folder in this repo). Mounted read-only at `/ref`; the tool finds `REFSEQ_*` / `*.exons.merged.bed` there automatically via `COVERAGEMASTER_REF_DIR` unless you pass `--refseq` / `--exon-bed` explicitly.

Output files (`<prefix>.CMreport`, `.CMcalls`, `.CM.log`, `.CMpositives.pdf`) land back on your host wherever `<data_dir>` points, since `/data` is a live mount, not copied out of the container afterward.

## Try it on the bundled demo case

```bash
./run.sh ./DEMO ./REF \
    /data/test.PGM1.cov /data/test.PGM1.report.txt PGM1 \
    -c /data/control.PGM1.cov -r /data/ref.PGM1 \
    -o /data/test_v2.PGM1 --genome hg19 -f
```

Expect a single `DEL` call in `PGM1` (~chr1:64,125,260-64,240,084) — this matches the
known-good result from the legacy tool on the same demo data. Verified against both
this hg19 demo case and real hg38 production data on 2026-08-20.

## docker-compose alternative

```bash
DATA_DIR=/path/to/your/data REF_DIR=./REF \
  docker compose run coveragemaster \
    /data/sample.cov /data/sample.report.txt GENE \
    -c /data/control.cov -r /data/ref_file -o /data/out
```

## Entry points (installed inside the image)

- `coveragemaster` — the main CNV caller (what `ENTRYPOINT` runs)
- `coveragemaster-build-ref` — builds the per-chromosome reference `.res` files from a folder of control `.cov`/`.report.txt` pairs
- `coveragemaster-compute-stats` — computes per-position mean/std from the concatenated reference

## Notes

- The `REF/` folder (RefSeq gene annotations + exon BED, ~100MB+) is deliberately **not** baked into the image — it's mounted at runtime via `-v ref_dir:/ref`, so the image stays small and you can swap hg19/hg38 references without rebuilding.
- `--genome hg19` / `--genome hg38` selects which RefSeq/exon reference set to look for by default; override with `--refseq` / `--exon-bed` for a custom build.
- This package (`coveragemaster/`, `pyproject.toml`, `Dockerfile`, etc.) is a separate, independent rewrite from the legacy `coverageMaster_hg38.py` / `libCoverageMaster_hg38.py` scripts documented in `README.md` — it does not share code with them, so fixes made to one side don't automatically apply to the other.
