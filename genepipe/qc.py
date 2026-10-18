"""
genepipe.qc -- quality-control metrics for sequence data.

Quality control is the first analytical stage of the pipeline because every
downstream result inherits its inputs' problems. The module computes the
standard per-sequence metrics (length, GC, ambiguity, composition) and applies
configurable pass/fail thresholds, returning a machine-readable table plus a
summary that the report stage renders directly.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Sequence

from .io import SequenceRecord
from .utils import get_logger

log = get_logger("qc")

# IUPAC ambiguity codes. Note the overlap with the amino-acid alphabet: M, S, W,
# B, D, H and V are all valid residues, so the set MUST depend on the alphabet --
# applying the nucleotide set to a protein flags most of the sequence as
# ambiguous (a bug this pipeline hit and now guards against).
AMBIGUOUS_DNA = set("NRYKMSWBDHV")
AMBIGUOUS_PROTEIN = set("BJZXUO*")


@dataclass
class QCThresholds:
    """Configurable acceptance criteria; defaults are sane for barcode data.

    ``check_gc`` exists because GC content is a meaningful quality signal for
    nucleotides but meaningless for proteins (where G and C are just glycine and
    cysteine), so the check is disabled automatically for protein input.
    """
    min_length: int = 100
    max_length: int = 100_000
    max_ambiguous_fraction: float = 0.05
    min_gc: float = 0.20
    max_gc: float = 0.70
    check_gc: bool = True
    alphabet: str = "dna"           # dna | protein -- selects the ambiguity set


@dataclass
class QCRecord:
    id: str
    length: int
    gc_pct: float
    ambiguous: int
    ambiguous_pct: float
    n_counts: dict[str, int]
    passed: bool
    reasons: list[str] = field(default_factory=list)


def gc_content(seq: str) -> float:
    """GC as a percentage of A+C+G+T (ambiguities excluded from the denominator)."""
    acgt = sum(seq.count(b) for b in "ACGT")
    if acgt == 0:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C")) / acgt


def qc_record(rec: SequenceRecord, th: QCThresholds) -> QCRecord:
    """Compute metrics for one record and apply the thresholds."""
    amb_set = AMBIGUOUS_PROTEIN if th.alphabet == "protein" else AMBIGUOUS_DNA
    n_amb = sum(rec.seq.count(c) for c in amb_set)
    amb_pct = 100.0 * n_amb / rec.length if rec.length else 0.0
    gc = gc_content(rec.seq)
    reasons: list[str] = []
    if rec.length < th.min_length:
        reasons.append(f"length {rec.length} < {th.min_length}")
    if rec.length > th.max_length:
        reasons.append(f"length {rec.length} > {th.max_length}")
    if amb_pct > th.max_ambiguous_fraction * 100:
        reasons.append(f"ambiguous {amb_pct:.1f}% > {th.max_ambiguous_fraction*100:.0f}%")
    if th.check_gc and not (th.min_gc * 100 <= gc <= th.max_gc * 100):
        reasons.append(f"GC {gc:.1f}% outside [{th.min_gc*100:.0f}, {th.max_gc*100:.0f}]")
    return QCRecord(id=rec.id, length=rec.length, gc_pct=round(gc, 2),
                    ambiguous=n_amb, ambiguous_pct=round(amb_pct, 2),
                    n_counts=dict(Counter(rec.seq)), passed=not reasons, reasons=reasons)


def run_qc(records: Sequence[SequenceRecord],
           th: QCThresholds | None = None) -> tuple[list[QCRecord], dict]:
    """QC every record and return (per-record results, aggregate summary)."""
    th = th or QCThresholds()
    results = [qc_record(r, th) for r in records]
    passed = [r for r in results if r.passed]
    lengths = [r.length for r in results] or [0]
    summary = {
        "n_input": len(results),
        "n_passed": len(passed),
        "n_failed": len(results) - len(passed),
        "pass_rate_pct": round(100 * len(passed) / len(results), 1) if results else 0.0,
        "mean_length": round(sum(lengths) / len(lengths), 1),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_gc_pct": (round(sum(r.gc_pct for r in results) / len(results), 2)
                        if results and th.check_gc else None),
        "total_ambiguous": sum(r.ambiguous for r in results),
        "failure_reasons": dict(Counter(x for r in results for x in r.reasons)),
    }
    gc_text = (f"mean GC {summary['mean_gc_pct']:.1f}%"
               if summary["mean_gc_pct"] is not None else "GC check not applicable")
    log.info("QC: %d/%d passed (%.1f%%), mean length %.1f, %s",
             summary["n_passed"], summary["n_input"], summary["pass_rate_pct"],
             summary["mean_length"], gc_text)
    return results, summary
