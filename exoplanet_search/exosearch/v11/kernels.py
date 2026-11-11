"""Reproducible numerics for v1.1 (docs/preregistration/tie_break.md, "Arithmetic").

* one thread for every math library;
* wotan's numba kernels recompiled with fastmath=False, cache=False, parallel=False (fastmath permits reassociation,
  and kernels loaded from wotan's on-disk cache differed on arm64; docs/host_gate/README.md);
* a per-process empty NUMBA_CACHE_DIR, so nothing compiled elsewhere is ever loaded.
Call install() before the first detrend in every worker process (idempotent).
"""
from __future__ import annotations

import os
import tempfile

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"

_STATE = {"n": None}


def install():
    if _STATE["n"] is not None:
        return _STATE["n"]
    try:                                  # numpy may already be loaded (e.g. by a caller): cap its BLAS pool too
        from threadpoolctl import threadpool_limits
        _STATE["tp"] = threadpool_limits(1)
    except Exception:
        pass
    os.environ["NUMBA_CACHE_DIR"] = tempfile.mkdtemp(prefix="nbc-v11-")
    import importlib
    import pkgutil

    import numba
    import wotan
    from numba.core.registry import CPUDispatcher
    mods = [importlib.import_module(f"wotan.{m.name}") for m in pkgutil.iter_modules(wotan.__path__)]
    new = {}
    for m in mods:                                   # one replacement per original kernel, shared by every module
        for k, v in list(vars(m).items()):
            if isinstance(v, CPUDispatcher):
                if id(v) not in new:
                    opts = dict(v.targetoptions)
                    opts.update(fastmath=False, cache=False, parallel=False)
                    new[id(v)] = numba.jit(**opts)(v.py_func)
                setattr(m, k, new[id(v)])
    _STATE["n"] = len(new)
    return _STATE["n"]


def installed():
    return _STATE["n"] is not None
