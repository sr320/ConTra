#!/usr/bin/env python3
"""DEPRECATED: Use `context_dependent_analysis.py` instead.

This file is retained only for backward compatibility. The functionality of
both the former full and subset analyses has been unified into
`context_dependent_analysis.py` with a --mode flag (full/subset) or
interactive prompt.

Usage (subset mode):
    python3 context_dependent_analysis.py --mode subset

Usage (full mode):
    python3 context_dependent_analysis.py --mode full

The old implementation has been removed to avoid divergence.
"""

import sys
from pathlib import Path

here = Path(__file__).resolve().parent
unified = here / 'context_dependent_analysis.py'

msg = (
    "This script is deprecated. Please run the unified analysis instead:\n\n"
    "  python3 context_dependent_analysis.py --mode subset  # or --mode full\n\n"
    f"Unified script path: {unified}\n"
)

print(msg)
# Optional: exit with non-error code
sys.exit(0)
