# genePipe — An Automated Pipeline for Genetic Data Analysis

**Genetics Assistant Internship — Week 4** | Rabi Narayan Patra | 29 September 2026

A working, tested, resumable command-line pipeline that takes a FASTA file through **cleaning → quality control → alignment → variant calling → statistical association testing → phylogenetic reconstruction** and writes a documented run report, in one command.

```bash
python3 -m genepipe run --input data/sample.fasta --out results/run1 --groups data/groups.json
```

## Why it exists

Genetic sequence analysis is a chain of repetitive, order-dependent, easy-to-get-wrong steps. genePipe automates that chain, and — more importantly — makes its own behaviour auditable: every exclusion, fallback and assumption is logged and written to a run manifest, and every stage is covered by tests.

## Results it produced

**On a simulated cohort with a known ground truth** (22 sequences, 1.5 kb, seeded and reproducible):

| Stage | Outcome |
|---|---|
| cleaning | 24 → 22 records (1 masked, 1 too short, 1 duplicate — all reported) |
| QC | 22/22 pass, mean 1,499.8 bp, mean GC 50.3% |
| alignment | 22 × 1,505 columns (FAMSA), gap fraction 0.35% |
| variants | 161 events (151 SNV, 4 ins, 6 del) |
| association | **141 sites tested → 1 significant** |
| phylogeny | NJ tree + 200 bootstrap replicates |

The one significant site is the **planted marker**: `g.612C>G`, present in **10/10 cases and 0/10 controls**, raw p = 1.08e-05, adjusted p = **0.0015**. Recovering a known answer exactly is the strongest evidence the preceding stages worked.

**On real public data** (curated COX1 proteins, 12 vertebrates, fetched from UniProt by a script in this repo) — the same pipeline, only the config changed:

| Clade recovered | Bootstrap support |
|---|---|
| *Homo sapiens* + *Pan troglodytes* | 100% |
| *Gallus gallus* + *Columba livia* | 100% |
| *Danio rerio* + *Salmo salar* | 99.8% |
| mammals (5 species) | 90.6% |

## Engineering

- **12 modules, 4 layers** — entry points → orchestrator → analysis → output, plus shared utilities and config
- **86 automated tests** (`pytest`), one file per module, including end-to-end and caching tests
- **Stage cache** — fingerprints each stage's inputs and config; a repeat run goes from **17.2 s to 0.0 s**
- **Fail-safe orchestration** — a stage that raises is recorded as failed, independent stages still run, and the report is always written
- **No third-party dependencies in the core** — runs on a bare Python 3.10+ interpreter; FAMSA and PyYAML are optional accelerators
- **Provenance** — pipeline version, interpreter, platform, SHA-256 input fingerprint and a configuration snapshot recorded per run
- **Deterministic output** — seeded randomness and fixed ordering, so results can be diffed and checksummed

## Defects found and fixed during development

Eleven substantive defects were found by testing — each one a case where the pipeline would have produced output that *looked* correct:

| # | Symptom | Root cause |
|---|---|---|
| 1 | Protein sequences reported 90% ambiguous; masking destroyed them | the DNA ambiguity set (N/R/Y/K/**M/S/W/B/D/H/V**) overlaps the amino-acid alphabet |
| 2 | Tree building crashed on `NoneType` | K2P distance is undefined at saturation → NaN → no pair ever selected |
| 3 | A resumed run failed every downstream stage | the cache skipped computation without restoring state |
| 4 | Controls counted 9 instead of 7; p-value wrong | the association denominator included samples dropped by cleaning |
| 5 | Every default run ended with a failed "report" stage | "report" was listed as a runnable stage but has no handler |
| 6 | `run --config x.yaml` refused to start | `--input`/`--out` were declared required |
| 7 | `--log-file` rejected after the subcommand | option defined only on the top-level parser |
| 8 | QC failed every protein record on GC content | GC is a nucleotide metric, meaningless for proteins |
| 9 | Ts/Tv reported for protein data | ts/tv is a nucleotide concept, computed on amino acids |
| 10 | Logging error: `%f format: a number is required, not NoneType` | a log statement formatted a value that is `None` when the GC check is off |
| 11 | Cleaning to zero records produced an empty "successful" run | no stage checked that any records survived |

Each is now fixed **and guarded by a regression test**.

## Repository contents

**Committed in this repository** (the shared layer, the IO layer, the cleaning/QC layer and the reporting layer):

```
├── genepipe/
│   ├── __init__.py        # public API
│   ├── utils.py           # logging, timing, checksums, versioning
│   ├── config.py          # layered configuration with validation
│   ├── io.py              # FASTA parsing/validation (gzip aware)
│   ├── clean.py           # primer trimming, masking, filtering, de-duplication
│   ├── qc.py              # quality metrics and thresholds (alphabet-aware)
│   ├── seqops.py          # translation, ORF finding, composition
│   └── report.py          # Markdown / HTML / JSON reporting
├── config.example.yaml    # documented configuration
├── run_all.sh             # reproduces everything end to end
├── requirements.txt       # core needs nothing; optional extras listed
└── README.md
```

**In the [v1.0 release](https://github.com/rabinarayanpatra350-biotech/genetics-week4-automation-pipeline/releases/tag/v1.0) reproducibility bundle** — the complete project:

```
├── genepipe/              # all 12 modules: adds cli.py, pipeline.py, align.py,
│                          #   variants.py, stats.py and phylo.py to those above
├── tests/                 # the 86-test pytest suite (one file per module)
├── scripts/               # sample-data generator, public-data fetcher, figures, report builder
├── data/                  # simulated cohort + real public sequences with provenance
├── results/run_sample/    # complete output of the simulated run (log, artefacts, reports)
├── results/run_real/      # complete output of the real-data run
├── config.real.yaml       # configuration for the real dataset
└── the report and all 8 figures
```

## Reproduce it

```bash
python3 scripts/make_sample_data.py          # seeded simulated cohort
python3 scripts/fetch_public_data.py         # real COX1 sequences from UniProt
python3 -m pytest tests -q                   # 86 passed
python3 -m genepipe run --input data/sample.fasta --out results/run_sample \
    --reference reference_SIMG1 --groups data/groups.json --bootstrap 200
python3 -m genepipe run --config config.real.yaml
```

The full report, figures and reproducibility bundle are attached to the [v1.0 release](https://github.com/rabinarayanpatra350-biotech/genetics-week4-automation-pipeline/releases/tag/v1.0).

## Author

**Rabi Narayan Patra** — Genetics Assistant Internship (2026). GitHub: [@rabinarayanpatra350-biotech](https://github.com/rabinarayanpatra350-biotech)
