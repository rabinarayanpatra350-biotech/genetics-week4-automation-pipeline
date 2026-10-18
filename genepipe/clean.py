"""
genepipe.clean -- data-cleaning operations.

The task explicitly asks for "data cleaning" to be automated, and in sequence
work that means: removing technical sequence (primers/adapters) that would
otherwise be mistaken for biology, masking or dropping ambiguous positions,
filtering records that fail QC, de-duplicating identical records, and
normalising identifiers so that downstream joins between files are reliable.

Every operation is pure (it returns new records and a report) so that it can be
tested in isolation and composed in any order by the orchestrator.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

from .io import SequenceRecord
from .utils import get_logger

log = get_logger("clean")

AMBIGUOUS_DNA = set("NRYKMSWBDHV")
AMBIGUOUS_PROTEIN = set("BJZXUO*")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class CleanReport:
    """What each cleaning step changed -- reported, never silent."""
    trimmed_primers: dict[str, int] = field(default_factory=dict)
    masked_records: list[str] = field(default_factory=list)
    dropped_short: list[str] = field(default_factory=list)
    dropped_qc: list[str] = field(default_factory=list)
    duplicates_removed: list[str] = field(default_factory=list)
    renamed: dict[str, str] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "n_trimmed": len(self.trimmed_primers),
            "n_masked": len(self.masked_records),
            "n_dropped_short": len(self.dropped_short),
            "n_dropped_qc": len(self.dropped_qc),
            "n_duplicates_removed": len(self.duplicates_removed),
            "n_renamed": len(self.renamed),
        }


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def trim_primers(records: Sequence[SequenceRecord], five_prime: str = "",
                 three_prime: str = "") -> tuple[list[SequenceRecord], dict[str, int]]:
    """Remove known primer/adapter sequences from the ends of each read.

    Matching is case-insensitive and tolerant of the common situation where the
    primer is present in reverse-complement orientation. Only the *anchored*
    case is trimmed (the primer must sit at the very end), because trimming an
    internal match would delete real sequence.
    """
    out: list[SequenceRecord] = []
    trimmed: dict[str, int] = {}
    fwd5, rev3 = five_prime.upper(), reverse_complement(three_prime.upper()) if three_prime else ""
    for rec in records:
        s = rec.seq
        removed = 0
        if fwd5 and s.startswith(fwd5):
            s, removed = s[len(fwd5):], removed + len(fwd5)
        # the 3' primer may be supplied in either orientation, so accept both
        for candidate in (rev3, three_prime.upper()):
            if candidate and s.endswith(candidate):
                s, removed = s[:-len(candidate)], removed + len(candidate)
                break
        if removed:
            trimmed[rec.id] = removed
            log.debug("trimmed %d nt of primer from %s", removed, rec.id)
        out.append(SequenceRecord(rec.id, s, rec.description))
    if trimmed:
        log.info("primer trimming removed sequence from %d record(s)", len(trimmed))
    return out, trimmed


def mask_ambiguous(records: Sequence[SequenceRecord], mask: str = "N",
                   alphabet: str = "dna") -> tuple[list[SequenceRecord], list[str]]:
    """Replace ambiguous codes with a single symbol so distance models stay valid.

    The ambiguity set is alphabet-specific: masking a protein with the nucleotide
    set would destroy legitimate M/S/W/B/D/H/V residues.
    """
    amb = AMBIGUOUS_PROTEIN if alphabet == "protein" else AMBIGUOUS_DNA
    out, masked = [], []
    for rec in records:
        if any(c in amb for c in rec.seq):
            masked.append(rec.id)
        out.append(SequenceRecord(rec.id, "".join(mask if c in amb else c for c in rec.seq),
                                  rec.description))
    if masked:
        log.info("masked ambiguous positions in %d record(s)", len(masked))
    return out, masked


def filter_by_length(records: Sequence[SequenceRecord], min_length: int
                     ) -> tuple[list[SequenceRecord], list[str]]:
    """Drop records shorter than min_length (too short to align or call variants)."""
    keep = [r for r in records if r.length >= min_length]
    dropped = [r.id for r in records if r.length < min_length]
    if dropped:
        log.warning("dropped %d record(s) shorter than %d bp", len(dropped), min_length)
    return keep, dropped


def deduplicate(records: Sequence[SequenceRecord]
                ) -> tuple[list[SequenceRecord], list[str]]:
    """Remove exact duplicate sequences, keeping the first occurrence.

    Duplicates are common when the same accession is retrieved twice, and they
    silently bias any downstream statistic that treats sequences as independent.
    """
    seen: dict[str, str] = {}
    keep, dupes = [], []
    for rec in records:
        if rec.seq in seen:
            dupes.append(rec.id)
        else:
            seen[rec.seq] = rec.id
            keep.append(rec)
    if dupes:
        log.warning("removed %d duplicate sequence(s): %s", len(dupes), ", ".join(dupes[:5]))
    return keep, dupes


def normalise_ids(records: Sequence[SequenceRecord]
                  ) -> tuple[list[SequenceRecord], dict[str, str]]:
    """Make identifiers filesystem- and join-safe (no spaces/pipes/colons)."""
    out, renamed = [], {}
    for rec in records:
        safe = "".join(ch if (ch.isalnum() or ch in "._-") else "_" for ch in rec.id)
        if safe != rec.id:
            renamed[rec.id] = safe
        out.append(SequenceRecord(safe, rec.seq, rec.description))
    if renamed:
        log.info("normalised %d identifier(s)", len(renamed))
    return out, renamed
