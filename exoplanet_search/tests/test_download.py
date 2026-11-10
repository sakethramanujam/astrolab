import os
import tempfile
from pathlib import Path

import pytest

from exosearch import targets


def test_valid_fits_detection(tmp_path):
    p = tmp_path / "x_lc.fits"
    p.write_bytes(b"SIMPLE" + b" " * (2880 - 6))
    assert targets._valid_fits(p)
    p.write_bytes(b"garbage")
    assert not targets._valid_fits(p)


@pytest.mark.skipif(os.environ.get("EXO_OFFLINE") == "1", reason="offline")
def test_manifest_and_idempotent_download(tmp_path, monkeypatch):
    from exosearch import config
    monkeypatch.setattr(config, "DATA_ROOT", Path(os.environ.get("EXO_DATA_ROOT", "/media/nfs01/astro_data/tess")))
    man = targets.fetch_manifest(70)
    assert len(man) > 1000 and man.tic.is_unique
    sel = targets.select_targets(man, 5, seed=1)
    assert sel.equals(targets.select_targets(man, 5, seed=1))  # deterministic
    r = sel.iloc[0]
    p1 = targets.download_one(r.fname, r.url, tmp_path)
    m1 = p1.stat().st_mtime_ns
    p2 = targets.download_one(r.fname, r.url, tmp_path)  # cached: no re-download
    assert p1 == p2 and p2.stat().st_mtime_ns == m1
