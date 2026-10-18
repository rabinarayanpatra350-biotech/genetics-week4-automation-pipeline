"""
genepipe.seqops -- sequence arithmetic: composition, translation, ORF finding.

These are the primitive operations the analytical stages are built from. They
are deliberately dependency-free (no Biopython) so the pipeline installs and
runs anywhere, and they are the most heavily unit-tested part of the codebase
because every downstream result depends on them.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator

STANDARD_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
STOPS = {"TAA", "TAG", "TGA"}


def gc_content(seq: str) -> float:
    acgt = sum(seq.count(b) for b in "ACGT")
    return 0.0 if acgt == 0 else 100.0 * (seq.count("G") + seq.count("C")) / acgt


def composition(seq: str) -> dict[str, float]:
    """Nucleotide composition as percentages of the full length."""
    n = len(seq) or 1
    return {b: round(100.0 * seq.count(b) / n, 2) for b in "ACGT"}


def translate(seq: str, frame: int = 0, stop_symbol: str = "*") -> str:
    """Translate a nucleotide sequence in the given reading frame (0, 1, 2).

    Codons containing ambiguous bases are rendered as 'X' rather than being
    silently dropped, so a translation always has length floor(len/3) and
    coordinates in the protein correspond to coordinates in the DNA.
    """
    if frame not in (0, 1, 2):
        raise ValueError(f"frame must be 0, 1 or 2 (got {frame})")
    seq = seq[frame:]
    out = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        codon = seq[i:i + 3]
        aa = STANDARD_CODE.get(codon, "X")
        out.append(stop_symbol if aa == "*" else aa)
    return "".join(out)


@dataclass
class ORF:
    start: int
    end: int
    frame: int
    protein: str

    @property
    def length_nt(self) -> int:
        return self.end - self.start


def find_orfs(seq: str, min_aa: int = 30, require_start: bool = True,
              require_stop: bool = True) -> list[ORF]:
    """Find open reading frames on the forward strand.

    An ORF runs from a start codon to the first in-frame stop codon. Frames
    whose ORF never closes are reported only when require_stop is False, because
    an unterminated ORF in a fragment is usually an artefact of the fragment
    boundary rather than a real gene.
    """
    orfs: list[ORF] = []
    for frame in range(3):
        i = frame
        while i + 3 <= len(seq):
            codon = seq[i:i + 3]
            if codon in STOPS:
                i += 3
                continue
            if require_start and codon != "ATG":
                i += 3
                continue
            start = i
            j = i
            while j + 3 <= len(seq) and seq[j:j + 3] not in STOPS:
                j += 3
            protein = translate(seq[start:j], 0)
            closed = j + 3 <= len(seq) and seq[j:j + 3] in STOPS
            if len(protein) >= min_aa and (closed or not require_stop):
                end = j + 3 if closed else j
                orfs.append(ORF(start=start, end=end, frame=frame, protein=protein))
            i = j + 3
    return sorted(orfs, key=lambda o: o.length_nt, reverse=True)


def hamming(a: str, b: str) -> int:
    """Positional differences between two equal-length strings."""
    if len(a) != len(b):
        raise ValueError("hamming() requires equal-length sequences")
    return sum(1 for x, y in zip(a, b) if x != y)


def pairwise_identity(a: str, b: str) -> float:
    """Percent identity over the overlap of two sequences (no gaps assumed)."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return 100.0 * sum(1 for x, y in zip(a[:n], b[:n]) if x == y) / n
