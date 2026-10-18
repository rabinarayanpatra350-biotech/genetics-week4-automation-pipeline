"""
genepipe.utils -- shared infrastructure: logging, checksums, timing, versioning.

Every stage of the pipeline uses these helpers so that behaviour (log format,
timing, provenance) is identical across the codebase. This is deliberate: the
task asks for an *automated* pipeline, and automation is only trustworthy when
its bookkeeping is uniform.
"""
from __future__ import annotations

import hashlib
import logging
import os
import platform
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

__version__ = "1.0.0"
LOG_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)-18s | %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def setup_logging(level: str = "INFO", logfile: str | None = None,
                  quiet: bool = False) -> logging.Logger:
    """Configure root logging to console and (optionally) a file.

    Idempotent: calling it twice does not duplicate handlers, which matters
    because the CLI and the test-suite both call it.
    """
    root = logging.getLogger("genepipe")
    root.setLevel(getattr(logging, level.upper(), logging.INFO))
    root.handlers.clear()

    if not quiet:
        con = logging.StreamHandler(sys.stdout)
        con.setFormatter(logging.Formatter(LOG_FORMAT, DATE_FORMAT))
        root.addHandler(con)
    if logfile:
        os.makedirs(os.path.dirname(os.path.abspath(logfile)), exist_ok=True)
        fh = logging.FileHandler(logfile, encoding="utf-8")
        fh.setFormatter(logging.Formatter(LOG_FORMAT, DATE_FORMAT))
        root.addHandler(fh)
    root.propagate = False
    return root


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"genepipe.{name}")


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    """SHA-256 of a file, streamed so that large inputs do not exhaust memory."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


@dataclass
class StageResult:
    """Outcome of one pipeline stage: name, status, timing and artefacts."""
    name: str
    status: str = "pending"           # ok | skipped | failed
    seconds: float = 0.0
    message: str = ""
    outputs: list[str] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)


@contextmanager
def timer(result: StageResult):
    """Time a block and record the duration on the StageResult.

    Exceptions are recorded as a failed stage and re-raised, so the orchestrator
    can decide whether to abort or continue.
    """
    start = time.perf_counter()
    try:
        yield
        result.status = "ok" if result.status == "pending" else result.status
    except Exception as exc:                       # noqa: BLE001 - re-raised
        result.status = "failed"
        result.message = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        result.seconds = round(time.perf_counter() - start, 3)


def environment_report() -> dict[str, str]:
    """Capture the runtime environment for the provenance record.

    Reproducibility is a stated requirement of the task, so every run records
    the interpreter, platform and pipeline version.
    """
    return {
        "pipeline_version": __version__,
        "python": sys.version.split()[0],
        "platform": f"{platform.system()} {platform.release()}",
        "machine": platform.machine(),
        "executable": sys.executable,
    }
