"""Per-baseline receiver-crosstalk (coherent delay line) QA.  See core.py for the algorithm and NOTES.md ("Crosstalk delay QA").

    from crosstalk import scan
    res = scan("obs.ms", metafits="obs.metafits"); res.flagged; res.tiles
"""

import os as _os
import sys as _sys

# Many tiny LAPACK calls: BLAS thread pools only oversubscribe (100x slower on a loaded box).
if "numpy" not in _sys.modules:
    for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        _os.environ.setdefault(_v, "1")

from .api import Result, analyse, scan  # noqa: E402
from .core import BandwidthError, Params  # noqa: E402

__all__ = ["scan", "analyse", "Result", "Params", "BandwidthError"]
