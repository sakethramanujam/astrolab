import numpy as np

from exosearch.v11 import beta_select as bs


def test_tost_gate_pass_and_fail():
    assert bs.tost([-0.3, 0.4])
    assert bs.tost([-0.5, 0.5])                       # boundary inside
    assert not bs.tost([-0.6, 0.1])                   # lower end outside
    assert not bs.tost([0.2, 0.7])                    # excludes 0 is irrelevant; upper end outside
    assert bs.tost([0.1, 0.3])                        # CI excludes 0 but is equivalent: passes


def _ind(p, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    return (rng.uniform(size=n) < p)


def test_ties_resolve_by_fixed_order():
    # identical indicators -> every option tied -> 'a' wins by order a, d, c, b
    base = _ind(0.6)
    ind = {o: base.copy() for o in "abcd"}
    comp = {o: {"vetted": float(ind[o].mean())} for o in "abcd"}
    o, _, tied = bs.pick(list("abcd"), comp, bs.paired_boot(ind, n_boot=300))
    assert o == "a" and set(tied) == set("abcd")
    # without a: d before c before b
    o, _, _ = bs.pick(list("bcd"), comp, bs.paired_boot(ind, n_boot=300))
    assert o == "d"
    o, _, _ = bs.pick(list("bc"), comp, bs.paired_boot(ind, n_boot=300))
    assert o == "c"


def test_clear_leader_wins_when_not_tied():
    base = _ind(0.5, seed=1)
    better = base | _ind(0.3, seed=2)                 # strictly more vetted, paired: CI excludes 0
    ind = {"a": base, "b": better, "c": base, "d": base}
    comp = {o: {"vetted": float(v.mean())} for o, v in ind.items()}
    o, _, tied = bs.pick(list("abcd"), comp, bs.paired_boot(ind, n_boot=300))
    assert o == "b" and tied == ["b"]


def test_no_option_passes_gives_b():
    o, why, _ = bs.pick([], {}, {})
    assert o == "b" and "not launched" in why
