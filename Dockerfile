FROM python:3.7-slim

WORKDIR /app

# Pinned to match the versions validated against this exact script on the
# production host (bare-metal CentOS 7 / Python 3.6). Newer matplotlib drops
# the top-level `pylab` module this script imports, and newer numpy/scipy
# can shift floating-point behavior in the HMM/wavelet calls - pin, don't
# float, unless you re-validate output against known-good results.
RUN pip install --no-cache-dir \
    numpy==1.16.2 \
    sympy==1.0 \
    PyWavelets==1.0.3 \
    matplotlib==2.2.3 \
    scipy==1.2.1 \
    more-itertools==8.8.0 \
    pandas==1.0.0

COPY coverageMaster_hg38.py libCoverageMaster_hg38.py HMM_CM.py ReadCount.py DGVExplorer.py log.py ./
# The hg38 gene annotation + exon reference are loaded at import time from
# ./REF relative to libCoverageMaster_hg38.py itself - small enough (~11MB)
# to bake in rather than require a volume mount for a fixed dependency.
COPY REF/REFSEQ_hg38_c_HGMD3.gz REF/hg38.exons.merged.bed REF/

ENTRYPOINT ["python3", "/app/coverageMaster_hg38.py"]
