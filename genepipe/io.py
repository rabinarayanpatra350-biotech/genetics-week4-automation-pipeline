"""
genepipe.io -- robust FASTA input/output.

Design decisions
----------------
* The reader is *strict but forgiving*: it validates that each record has a
  header and non-empty sequence, normalises sequence case, and reports the
  problems it found instead of raising on the first bad record. A pipeline that
  dies on one malformed record in a 10,000-record file is not automating
  anything.
* gzip input is transparently supported because public sequence files are
  almost always distributed compressed.
* Writing is deterministic (fixed line width, fixed record order) so that
  outputs can be checksummed and compared between runs.
"""
from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, field
from typing import Iterable, Iterator

from .utils import get_logger

log = get_logger("io")

VALID_DNA = set("ACGTNRYKMSWBDHVU-")
VALID_PROTEIN = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")


@dataclass
class SequenceRecord:
    """One FASTA record with a stable identifier and a description."""
    id: str
    seq: str
    description: str = ""

    @property
    def length(self) -> int:
        return len(self.seq)

    def __post_init__(self) -> None:
        self.seq = self.seq.upper().replace(" ", "").replace("\n", "")


@dataclass
class FastaValidation:
    """Problems found while reading a FASTA file (empty means clean)."""
    n_records: int = 0
    empty_sequences: list[str] = field(default_factory=list)
    duplicate_ids: list[str] = field(default_factory=list)
    invalid_characters: dict[str, str] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return not (self.empty_sequences or self.duplicate_ids
                    or self.invalid_characters)


def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def read_fasta(path: str, alphabet: str = "dna") -> tuple[list[SequenceRecord], FastaValidation]:
    """Read a (optionally gzipped) FASTA file.

    Returns the parsed records *and* a validation report. Records with an empty
    sequence are dropped from the return value but recorded in the report, so
    the caller can choose to warn, fail, or continue.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"input FASTA not found: {path}")
    allowed = VALID_DNA if alphabet == "dna" else VALID_PROTEIN

    records: list[SequenceRecord] = []
    report = FastaValidation()
    seen: set[str] = set()
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal header, chunks
        if header is None:
            return
        seq = "".join(chunks).upper().replace(" ", "")
        rid = header.split()[0]
        desc = header[len(rid):].strip()
        report.n_records += 1
        if rid in seen:
            report.duplicate_ids.append(rid)
        seen.add(rid)
        if not seq:
            report.empty_sequences.append(rid)
        else:
            bad = {c for c in seq if c not in allowed}
            if bad:
                report.invalid_characters[rid] = "".join(sorted(bad))
            records.append(SequenceRecord(id=rid, seq=seq, description=desc))
        header, chunks = None, []

    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
            else:
                if header is None:
                    raise ValueError(f"{path}: sequence data before the first header")
                chunks.append(line)
        flush()

    log.info("read %d records from %s (%s)", len(records), path,
             "clean" if report.ok else "with warnings")
    if not report.ok:
        log.warning("input issues: %d empty, %d duplicate id(s), %d with invalid characters",
                    len(report.empty_sequences), len(report.duplicate_ids),
                    len(report.invalid_characters))
    return records, report


def write_fasta(records: Iterable[SequenceRecord], path: str, width: int = 60) -> str:
    """Write records to FASTA deterministically; returns the path written."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for rec in records:
            head = rec.id if not rec.description else f"{rec.id} {rec.description}"
            fh.write(f">{head}\n")
            for i in range(0, len(rec.seq), width):
                fh.write(rec.seq[i:i + width] + "\n")
    log.info("wrote %d records -> %s", len(list(records)) if isinstance(records, list) else -1, path)
    return path


def iter_fasta(path: str) -> Iterator[SequenceRecord]:
    """Streaming reader for very large files (memory-safe)."""
    header: str | None = None
    chunks: list[str] = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    rid = header.split()[0]
                    yield SequenceRecord(rid, "".join(chunks), header[len(rid):].strip())
                header, chunks = line[1:].strip(), []
            elif line:
                chunks.append(line)
        if header is not None:
            rid = header.split()[0]
            yield SequenceRecord(rid, "".join(chunks), header[len(rid):].strip())
