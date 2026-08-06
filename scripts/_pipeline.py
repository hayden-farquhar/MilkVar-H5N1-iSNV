"""
Loader giving the later analysis scripts importable access to the numbered
pipeline modules.

The pipeline scripts are numbered to document their run order, but a leading
digit is not a valid Python identifier, so `import 11_adaptation_analysis` is a
syntax error. Scripts 17 onward reuse the site definitions, codon classification,
and caller parsers from the earlier stages rather than restating them, so they
load those modules here by file path instead of by name.

Reusing rather than restating matters: a sensitivity analysis that re-derives its
own site coordinates or codon logic is not testing the primary analysis, it is
testing a different analysis that happens to resemble it.

Usage:
    from _pipeline import SITE_DEFINITIONS, classify_variant, parse_ivar
"""

import importlib.util
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent


def _load(stem: str, alias: str):
    """Import a numbered script as a module under an importable alias."""
    path = _SCRIPTS_DIR / f"{stem}.py"
    if not path.exists():
        raise FileNotFoundError(f"expected pipeline script not found: {path}")
    spec = importlib.util.spec_from_file_location(alias, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


adaptation_analysis = _load("11_adaptation_analysis", "adaptation_analysis")
concordance_filter = _load("04_concordance_filter", "concordance_filter")
merge_results = _load("05_merge_results", "merge_results")

# Re-exports used by scripts 17+
SITE_DEFINITIONS = adaptation_analysis.SITE_DEFINITIONS
classify_variant = adaptation_analysis.classify_variant
AF_THRESHOLD = adaptation_analysis.AF_THRESHOLD
DEPTH_THRESHOLD = adaptation_analysis.DEPTH_THRESHOLD
ALPHA = adaptation_analysis.ALPHA

parse_ivar = concordance_filter.parse_ivar
parse_lofreq = concordance_filter.parse_lofreq
strand_bias_filter = concordance_filter.strand_bias_filter

concordance = merge_results.concordance

__all__ = [
    "SITE_DEFINITIONS", "classify_variant", "AF_THRESHOLD", "DEPTH_THRESHOLD",
    "ALPHA", "parse_ivar", "parse_lofreq", "strand_bias_filter", "concordance",
    "adaptation_analysis", "concordance_filter", "merge_results",
]
