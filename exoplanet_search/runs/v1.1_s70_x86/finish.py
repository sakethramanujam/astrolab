# S70 x86 production scoring (UNVALIDATED run). ONE verdict for every row, real / null / inject alike
# (exosearch.vet_final.final_verdict = vet_amend.amended() + common_period + checks_logged; no re-vetting
# pass and no light-curve rebuilds, so deep 'secondary_scan' rows get status 'needs_lc' and are bracketed:
# lower bound counts them as rejected, upper bound as exempted). Known-planet gate at 13.71 (frozen list),
# flipped-null check on known hosts, null counts above S*, candidates above S*, injection recovery scored
# with the same verdict. Writes report/ only; commits nothing.
import gzip, json, os, pandas as pd
from exosearch.v11 import run, known_planets as kp
from exosearch import vet, vet_final
P='/tmp/s70-v11-prod-x86-2026-11-12/'; R=P+'report/'; S=13.71
FROZEN='docs/preregistration/v1.1_known_planets_s70_2026-11-12_s13.71/frozen_eligibility_s70_s13.71.csv'


# ------------------------------------------------- the one verdict path (all modes)
def row_verdict(row, jobs=None, *, fail_common_period=None, scan=None, needs_lc='reject'):
    """final verdict for one row of ANY mode: vet_amend.amended() + common_period + checks_logged.
    jobs: manifest jobs by TIC (only the TIC density rho is used; no light curve is ever rebuilt)."""
    j = (jobs or {}).get(int(row['tic'])) or {}
    return vet_final.final_verdict(row, common_period=fail_common_period, checks_logged=row.get('checks_logged'),
                                   rho=j.get('rho'), scan=scan, needs_lc=needs_lc)


def real_path(row, **kw):
    """Verdict path for real rows (same function as every other mode)."""
    return row_verdict(row, **kw)


def null_path(row, **kw):
    """Verdict path for null rows (invert / scramble; same function as every other mode)."""
    return row_verdict(row, **kw)


def inject_path(row, **kw):
    """Verdict path for injection rows (same function as every other mode; no inject skip, no
    search-time passed_lc_checks any more)."""
    return row_verdict(row, **kw)


PATH_FOR_MODE = {'real': real_path, 'invert': null_path, 'scramble': null_path, 'inject': inject_path}


def path_for(row):
    return PATH_FOR_MODE.get(str(row.get('mode')), null_path)


def load_manifest_jobs(p):
    """Manifest jobs by TIC and the run cfg (rho per job feeds the amended vetter's a/R*)."""
    jobs = {}; cfg = None
    for m in ('known_h1', 'rest_h1'):
        with gzip.open(p + f'manifests/{m}.jsonl.gz', 'rt') as fh:
            cfg = json.loads(fh.readline())['header']['cfg']; jobs.update({j['tic']: j for j in map(json.loads, fh)})
    return jobs, cfg


def flag_common_periods_by_task(v):
    """common_period flag per (sector, task) group, exactly as before, for every mode."""
    out = []
    for _, g in v.groupby(['sector', 'task']):
        g = vet.flag_common_periods(g.assign(status='ok'))
        g['fail_common_period'] = g['flags'].fillna('').str.contains('common_period'); out.append(g)
    return pd.concat(out, ignore_index=True)


def main():
    os.makedirs(R, exist_ok=True)
    jobs, cfg = load_manifest_jobs(P)
    rows = pd.concat([run.load_rows(P + 'rows_known'), run.load_rows(P + 'rows_rest')], ignore_index=True)
    summ = dict(n_rows=int(len(rows)), n_ok=int((rows.status == 'ok').sum()), by_task=rows.task.value_counts().to_dict())
    ok = rows[rows.status == 'ok']
    # amended fields for every ok row (all modes; final_verdict_row never rebuilds a light curve).
    # The per-row verdict fields of this first pass are placeholders: the common-period term is only
    # known after the group flagging below, so they are dropped and re-attached there.
    v = pd.DataFrame([{k: x for k, x in vet_final.final_verdict_row(
                          r, checks_logged=r.get('checks_logged'), rho=(jobs.get(int(r['tic'])) or {}).get('rho')).items()
                       if k not in vet_final.VERDICT_KEYS} for r in ok.to_dict('records')])
    v = flag_common_periods_by_task(v)
    # final verdict columns: lower bound (needs_lc rejected) and upper bound (needs_lc exempted)
    lo, hi = [], []
    for r in v.to_dict('records'):
        kw = dict(jobs=jobs, fail_common_period=r['fail_common_period'])
        a = path_for(r)(r, **kw)
        b = path_for(r)(r, needs_lc='exempt', **kw)
        lo.append(dict(passed_final=a.passed, first_failing_rule=a.first_failing_rule,
                       all_failing=','.join(a.all_failing), status_final=a.status))
        hi.append(dict(passed_upper=b.passed, all_failing_upper=','.join(b.all_failing)))
    v = pd.concat([v, pd.DataFrame(lo), pd.DataFrame(hi)], axis=1)
    v['vetted_v2'] = v.passed_final.astype(bool) & (v.status_final == vet_final.STATUS_OK)
    v['vetted_v2_upper'] = v.passed_upper.astype(bool)
    v.to_parquet(R + 'rows_v2.parquet')
    real = v[v['mode'] == 'real'].drop(columns=['passed_lc_checks_v2']).rename(columns={'vetted_v2': 'passed_lc_checks_v2'})
    real.to_parquet(R + 'real_v2.parquet')
    df, ks = kp.run_known_planet_recovery(s_star=S, search_results_path=R + 'real_v2.parquet', out_dir=R + 'known', sector=70,
                                          frozen_eligibility_file=FROZEN)
    el_toi = df[(df.eligibility == 'eligible') & (df.type == 'toi')]
    ks['toi_eligible'] = int(len(el_toi)); ks['toi_recovered'] = int(el_toi.recovered.fillna(False).astype(bool).sum())
    summ['known_gate'] = ks
    fz = kp.read_frozen(FROZEN); el = fz[fz.eligibility == 'eligible']
    inv = v[v.task == 'invert'].set_index('tic')
    fl = []
    for r in el.itertuples():
        t = int(r.tic)
        if t not in inv.index: fl.append(dict(tic=t, name=r.name, invert_row=False)); continue
        x = inv.loc[t]; P_ = float(x.P_grid); tol = 0.01 * float(r.period)
        near = any(abs(P_ - k * float(r.period)) <= tol * k for k in (0.5, 1, 2))
        fl.append(dict(tic=t, name=r.name, required=bool(r.required), invert_S=float(x.S), invert_P=P_, invert_vetted=bool(x.vetted_v2),
                       at_catalog_period=bool(near), vanished=not (near and float(x.S) >= S and bool(x.vetted_v2))))
    fl = pd.DataFrame(fl); fl.to_csv(R + 'flipped_known_hosts.csv', index=False)
    summ['flipped_known_hosts'] = dict(n=int(len(fl)), not_vanished=int((~fl.vanished.fillna(True)).sum()) if 'vanished' in fl else None)
    summ['null_above_sstar'] = {t: dict(n=int((v.task == t).sum()), vetted_ge_S=int(((v.task == t) & v.vetted_v2 & (v.S >= S)).sum()),
        vetted_ge_S_upper=int(((v.task == t) & v.vetted_v2_upper & (v.S >= S)).sum()),
        all_ge_S=int(((v.task == t) & (v.S >= S)).sum())) for t in sorted(v.task.unique()) if t != 'real' and not t.startswith('inject')}
    known_tics = set(fz.tic.astype(int))
    cand = real[real.passed_lc_checks_v2 & (real.S >= S)].copy()
    cand['known_planet_host'] = cand.tic.astype(int).isin(known_tics)
    cand.sort_values('S', ascending=False)[['tic', 'S', 'P_grid', 'depth', 'duration', 'n_transits', 'rp_earth_v2', 'known_planet_host', 'warn_v2', 'holdout']].to_csv(R + 'candidates_above_sstar_INTERNAL.csv', index=False)
    summ['candidates_above_sstar'] = dict(total=int(len(cand)), on_known_planet_or_toi_hosts=int(cand.known_planet_host.sum()),
        other=int((~cand.known_planet_host).sum()), label='candidates above S* in an UNVALIDATED run; none is called a detection')
    inj_all = rows[rows['mode'] == 'inject']
    inj = v[v['mode'] == 'inject'].copy()          # ok injection rows, unified verdict columns attached
    inj['allowed'] = (inj.depth_inj < 0.05) & (inj.rp_inj < 25)
    inj['recovered'] = inj.matched.fillna(False).astype(bool) & (inj.S >= S)
    inj['vetted'] = inj['recovered'] & inj['vetted_v2']
    inj['vetted_upper'] = inj['recovered'] & inj['vetted_v2_upper']
    inj['in_bin_15_25'] = inj['allowed'] & (inj.snr_exp >= 15) & (inj.snr_exp < 25)
    inj['allowed_ge_S'] = inj['allowed'] & (inj.snr_exp >= S)
    def _inj_stat(g):
        n = int(len(g))
        return dict(n=n, matched=int(g.matched.fillna(False).astype(bool).sum()),
                    recovered_ge_S=int(g.recovered.sum()), vetted_ge_S=int(g.vetted.sum()),
                    vetted_ge_S_upper=int(g.vetted_upper.sum()), needs_lc=int((g.status_final == vet_final.STATUS_NEEDS_LC).sum()),
                    recovered_rate=round(float(g.recovered.mean()), 4) if n else None,
                    vetted_rate=round(float(g.vetted.mean()), 4) if n else None)
    summ['injection_recovery'] = dict(
        verdict='exosearch.vet_final.final_verdict (vet_amend.amended + common_period + checks_logged), same path as real and null rows',
        s_star=S, n_total=int(len(inj_all)), n_ok=int(len(inj)),
        n_unphysical=int((inj_all.status == 'unphysical').sum()),
        n_error=int((~inj_all.status.isin(['ok', 'unphysical'])).sum()),
        all_ok=_inj_stat(inj),
        allowed=_inj_stat(inj[inj.allowed]),
        allowed_snr_exp_ge_sstar=_inj_stat(inj[inj.allowed_ge_S]),
        allowed_snr_exp_15_25=_inj_stat(inj[inj.in_bin_15_25]),
        by_task={t: dict(n_total=int((inj_all.task == t).sum()), n_unphysical=int(((inj_all.task == t) & (inj_all.status == 'unphysical')).sum()),
                         **_inj_stat(inj[inj.task == t])) for t in sorted(inj_all.task.unique())},
        by_dens_class={c: _inj_stat(g) for c, g in inj.groupby('dens_class')},
        aliases={str(k): int(val) for k, val in inj.alias.dropna().value_counts().items()})
    if os.path.exists(P + 'thermal.json.log'):
        tl = [json.loads(l) for l in open(P + 'thermal.json.log') if l.strip()]
        summ['thermal_pauses'] = dict(
            pause_episodes=sum(1 for i, x in enumerate(tl) if x.get('action') == 'paused' and (i == 0 or tl[i - 1].get('action') != 'paused')),
            paused_ticks=sum(1 for x in tl if x.get('action') == 'paused'),
            resumed=sum(1 for x in tl if x.get('action') == 'resumed'),
            max_temp_c=max((x['temp_c'] for x in tl if x.get('temp_c') is not None), default=None))
    summ['status'] = 'UNVALIDATED (S72 gate failed 28/30)'
    json.dump(summ, open(R + 'summary.json', 'w'), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


if __name__ == '__main__':
    main()
