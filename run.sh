#!/usr/bin/env bash
# Wrapper to run coverageMaster in Docker, writing output to your local machine.
#
# Usage:
#   ./run.sh <data_dir> <ref_dir> [coveragemaster arguments...]
#
# Example (DEMO):
#   ./run.sh ./DEMO ./REF \
#       /data/test.PGM1.cov /data/test.PGM1.report.txt PGM1 \
#       -c /data/control.PGM1.cov -r /data/ref.PGM1 \
#       -o /data/test_v2.PGM1 --genome hg19 -f
#
# All paths after the first two arguments must use /data/... (mounted data dir)
# or /ref/... (mounted REF dir).

set -e

DATA_DIR=$(realpath "${1:?Usage: ./run.sh <data_dir> <ref_dir> [args...]}")
REF_DIR=$(realpath "${2:?Usage: ./run.sh <data_dir> <ref_dir> [args...]}")
shift 2

docker run --rm \
    -v "${DATA_DIR}:/data" \
    -v "${REF_DIR}:/ref" \
    -e COVERAGEMASTER_REF_DIR=/ref \
    coveragemaster \
    "$@"
