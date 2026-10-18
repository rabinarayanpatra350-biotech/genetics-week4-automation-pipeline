"""
genepipe.config -- configuration loading, defaults and validation.

Configuration is layered so that a run can be reproduced exactly from a file but
still overridden from the command line:

    built-in defaults  <  config file (YAML or JSON)  <  CLI overrides

Validation is explicit: an unknown key, a wrong type or a nonsensical value is
reported with the key name rather than silently ignored, because a pipeline that
ignores a typo in a config file is worse than one that crashes.
"""
from __future__ import annotations

import json
import os
from typing import Any

from .utils import get_logger

log = get_logger("config")

DEFAULTS: dict[str, Any] = {
    "input": None,
    "output_dir": "results",
    "alphabet": "dna",              # dna | protein
    "seed": 20260929,
    "resume": True,
    "stages": None,                 # None = all stages in STAGE_ORDER
    "clean": {
        "five_prime": "",
        "three_prime": "",
        "mask_ambiguous": True,
        "min_length": 100,
        "deduplicate": True,
    },
    "qc": {
        "min_length": 100,
        "max_length": 100000,
        "max_ambiguous_fraction": 0.05,
        "min_gc": 0.20,
        "max_gc": 0.70,
        "enforce": False,
    },
    "align": {
        "prefer_external": True,
    },
    "variants": {
        "reference_id": None,
        "distance_model": "k2p",
    },
    "stats": {"alpha": 0.05},
    "phylo": {"enabled": True, "bootstrap": 200, "distance_model": "k2p"},
    "groups": {},                   # sample id -> case/control
}

_KNOWN = set(DEFAULTS)


def _deep_merge(base: dict, override: dict) -> dict:
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def _load_file(path: str) -> dict:
    if not os.path.exists(path):
        raise FileNotFoundError(f"config file not found: {path}")
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if path.endswith((".yaml", ".yml")):
        try:
            import yaml                                  # type: ignore
        except ImportError as exc:                       # pragma: no cover
            raise ImportError(
                "PyYAML is required for YAML config files; "
                "use a .json config or `pip install pyyaml`") from exc
        return yaml.safe_load(text) or {}
    return json.loads(text)


def validate(cfg: dict) -> list[str]:
    """Return a list of human-readable problems (empty means valid)."""
    problems: list[str] = []
    unknown = set(cfg) - _KNOWN
    if unknown:
        problems.append(f"unknown configuration key(s): {', '.join(sorted(unknown))}")
    if not cfg.get("input"):
        problems.append("'input' is required")
    if cfg.get("alphabet") not in ("dna", "protein"):
        problems.append("'alphabet' must be 'dna' or 'protein'")
    for key in ("clean", "qc"):
        if not isinstance(cfg.get(key, {}), dict):
            problems.append(f"'{key}' must be a mapping")
    gc = cfg.get("qc", {})
    if gc.get("min_gc", 0) >= gc.get("max_gc", 1):
        problems.append("qc.min_gc must be less than qc.max_gc")
    if cfg.get("phylo", {}).get("bootstrap", 0) < 0:
        problems.append("phylo.bootstrap must be >= 0")
    return problems


def load_config(path: str | None = None, overrides: dict | None = None,
                strict: bool = True) -> dict:
    """Build the effective configuration from defaults, a file and overrides."""
    cfg = dict(DEFAULTS)
    if path:
        cfg = _deep_merge(cfg, _load_file(path))
        log.info("loaded configuration from %s", path)
    if overrides:
        cfg = _deep_merge(cfg, overrides)
    problems = validate(cfg)
    if problems and strict:
        raise ValueError("invalid configuration:\n  - " + "\n  - ".join(problems))
    for p in problems:
        log.warning("configuration: %s", p)
    return cfg


def dump_config(cfg: dict, path: str) -> str:
    """Persist the effective configuration so a run can be reproduced exactly."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=2, sort_keys=True)
    return path
