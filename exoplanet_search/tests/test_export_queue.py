from types import SimpleNamespace as NS

from exosearch import export, sectorqueue


def _r(**k):
    base = dict(tic=1, passed=0, toi_period_match=0, ctoi_period_match=0, planet_period_match=0,
                known_toi=None, known_planet=None, known_ctoi=None)
    base.update(k)
    base.setdefault("passed_screens_v10", base["passed"] == 1)   # export.load_rows derives it from `passed`
    return NS(**base)


def test_candidate_status():
    assert export._status(_r(passed=1), {}) == "needs_review"
    assert export._status(_r(passed=1, toi_period_match=1, known_toi="465.01"), {}) == "known"
    assert export._status(_r(), {}) == "rejected"
    assert export._status(_r(passed=1), {"1": {"where": "PHT"}}) == "submitted"
    assert export._label(_r(toi_period_match=1, known_toi="465.01")) == "TOI-465.01"


def test_queue_order(monkeypatch):
    monkeypatch.setattr(sectorqueue, "status", lambda: [
        dict(sector=1, status="done"), dict(sector=5, status="pending"),
        dict(sector=70, status="in_progress"), dict(sector=90, status="pending")])
    assert sectorqueue.next_sectors(3) == [70, 90, 5]
