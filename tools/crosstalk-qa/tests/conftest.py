import os
import sys

for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = os.path.dirname(__file__)
sys.path[:0] = [os.path.dirname(HERE), HERE]
