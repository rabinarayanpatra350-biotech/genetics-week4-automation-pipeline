"""
genePipe -- an automated pipeline for genetic data analysis.

Genetics Assistant Internship, Week 4 | Rabi Narayan Patra

A single command takes a FASTA file (or a directory of them) through cleaning,
quality control, alignment, variant calling, statistical association testing and
phylogenetic reconstruction, and produces a documented run report.

Typical use::

    from genepipe import Pipeline, load_config

    cfg = load_config("config.yaml", {"input": "data/sample.fasta",
                                      "output_dir": "results/run1"})
    manifest = Pipeline(cfg).run()

Or from the shell::

    python -m genepipe run --input data/sample.fasta --out results/run1
"""
from .config import load_config, validate
from .pipeline import Pipeline, STAGE_ORDER
from .utils import __version__

__all__ = ["Pipeline", "STAGE_ORDER", "load_config", "validate", "__version__"]
