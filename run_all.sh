#!/usr/bin/env bash
# run_all.sh -- reproduce the entire Week 4 deliverable from scratch.
#
#   1. generate the sample dataset (seeded, so identical every time)
#   2. run the test-suite
#   3. run the full pipeline on the sample data
#   4. run a second time to demonstrate the stage cache
#   5. run the pipeline on a real public dataset (COX1 across 12 vertebrates)
set -euo pipefail

PY=${PY:-python3}
echo "=== 1. generating sample data ==="
$PY scripts/make_sample_data.py

echo
echo "=== 2. running the test-suite ==="
$PY -m pytest tests -q

echo
echo "=== 3. full pipeline run on simulated data ==="
$PY -m genepipe run \
    --input data/sample.fasta \
    --out results/run_sample \
    --groups data/groups.json \
    --bootstrap 200 \
    --log-file results/run_sample/pipeline.log

echo
echo "=== 4. re-run (demonstrates the stage cache) ==="
$PY -m genepipe run \
    --input data/sample.fasta \
    --out results/run_sample \
    --groups data/groups.json \
    --bootstrap 200 | tail -3

echo
echo "=== 5. pipeline on real public data (COX1, 12 vertebrates) ==="
$PY -m genepipe run \
    --input data/cox1_proteins.fasta \
    --out results/run_real \
    --config config.real.yaml \
    --log-file results/run_real/pipeline.log

echo
echo "done -- reports are in results/run_sample/run_report.html and results/run_real/run_report.html"
