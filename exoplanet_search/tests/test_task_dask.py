"""End-to-end: one real sector-70 star through process_star on a 2-worker LocalCluster."""
import os

import pytest

from exosearch import targets
from exosearch.task import process_star


@pytest.mark.skipif(os.environ.get("EXO_OFFLINE") == "1", reason="offline")
def test_process_star_on_localcluster():
    from dask.distributed import Client, LocalCluster
    sel = targets.select_targets(targets.fetch_manifest(70), 2, seed=7)
    with LocalCluster(n_workers=2, threads_per_worker=1, host="127.0.0.1", dashboard_address=None) as c, Client(c) as cl:
        rows = cl.gather([cl.submit(process_star, r.tic, 70, r.fname, r.url, "auto", pure=False)
                          for r in sel.itertuples()])
    for row in rows:
        assert row["status"] == "ok", row.get("error")
        assert row["backend"] == "cpu" and row["period"] > 0.5 and row["sde"] > 0
