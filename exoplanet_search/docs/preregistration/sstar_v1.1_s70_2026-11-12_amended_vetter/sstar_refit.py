# S* refit on the x86 null rows: old vetter (beta_rows.with_vetting as is) and amended vetter (vet_amend rules a-d).
# Uses the repo's functions only: beta_rows reader/merge/with_vetting/beta_frames/sector_gate/decide/injection_completeness,
# vet_amend.revet. Leakage drop of the third registration (4 TICs, after pairing, alpha N unchanged) applied as recorded.
import os
for v in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS","NUMBA_NUM_THREADS","VECLIB_MAXIMUM_THREADS"): os.environ[v]="1"
import gzip, json, sys, time
import numpy as np, pandas as pd
from concurrent.futures import ProcessPoolExecutor
from exosearch.v11 import beta_rows as br, beta_select as bs, settings, vet_amend
from exosearch import quality
B='/tmp/s70-v11-nulls-x86-2026-11-12/'
OUT='/tmp/sstar-2026-11-12/'
LEAK=[243015498, 2052927695, 620213822, 630275176]
W=int(os.environ.get("W","10"))
t0=time.time()
c=settings.cfg(); rel=settings.exclusion_files(c)["70"]; ex=pd.read_csv(settings.ROOT/rel)
man=[B+'manifests/nulls_s70_all.jsonl.gz']
rows_df=br.read_dirs([B+'rows', B+'h3/tail/rows'])
x,rrep=br.merge_rows(rows_df, br.EXPECTED, None, 70, br.manifest_tasks(man), set(ex.tic.astype(int)), br.sha256(settings.ROOT/rel))
q=quality.load()
d=br.with_vetting(x)
# consistency: old vetted == row passed_lc_checks & no common period & checks logged
chk=(d.passed_lc_checks.astype(bool) & ~d.fail_common_period.astype(bool) & d.checks_logged.fillna(False).astype(bool))
print("old vetted vs passed_lc_checks&~cp&logged mismatches:", int((chk!=d.vetted).sum()), flush=True)
# amended verdict
with gzip.open(man[0],'rt') as fh:
    head=json.loads(fh.readline())['header']; jobs={int(j['tic']):j for j in map(json.loads,fh)}
cache=OUT+'revet_nulls.parquet'
if os.path.exists(cache):
    v=pd.read_parquet(cache)
else:
    v=vet_amend.revet(d.drop(columns=['vetted']), jobs, head['cfg'], workers=W)
    v.to_parquet(cache)
assert len(v)==len(d) and (v.tic.values==d.tic.values).all() and (v.task.values==d.task.values).all()
d2=d.copy()
for col in [c_ for c_ in v.columns if c_.endswith('_v2')]: d2[col]=v[col].values
d2['vetted']=d2.passed_lc_checks_v2.astype(bool) & ~d2.fail_common_period.astype(bool) & d2.checks_logged.fillna(False).astype(bool)
def run(dd, label):
    xf,xs,prep=br.beta_frames(dd,q,ex)
    prep['paired_before_leak_drop']=int(len(xf))
    keep=~xf.index.isin(LEAK); xf,xs=xf[keep],xs[keep]; prep['leakage_dropped']=LEAK; prep['paired']=int(len(xf))
    N=prep['n_alpha']
    with ProcessPoolExecutor(4) as exe:
        rr=list(exe.map(br._gate_job,[(xf,xs,o,N) for o in bs.OPTIONS]))
    res={(r['sector'],r['option']):r for r in rr}
    s_star={k:v_['s_star'] for k,v_ in res.items()}
    try:
        inj,comp,ind=br.injection_completeness(q,s_star); note="v1.0 interim sector-70 injections, scored at each option's S*"
    except FileNotFoundError as e:
        comp=ind=None; note=f"unavailable: {e}"
    ch=br.decide(res,comp,ind)
    return dict(label=label, per_sector_option=[res[(70,o)] for o in bs.OPTIONS], completeness_15_25=comp,
                completeness_note=note, decision=ch, population=prep,
                s_star_option_b=res[(70,'b')]['s_star'])
old=run(d,'old vetter (beta_rows.with_vetting)')
new=run(d2,'amended vetter (vet_amend rules a-d)')
chg=pd.DataFrame(dict(task=d.task, old=d.vetted, new=d2.vetted))
trans={t:{'fail_to_pass':int(((~g.old)&g.new).sum()),'pass_to_fail':int((g.old&~g.new).sum()),'vetted_old':int(g.old.sum()),'vetted_new':int(g.new.sum())} for t,g in chg.groupby('task')}
rep=dict(reader=rrep, verdict_changes=trans, old=old, new=new, n_revet_rebuilt=int(d.deep_flags.fillna('').str.contains('secondary_scan').sum()), runtime_s=time.time()-t0)
json.dump(rep, open(OUT+'sstar.json','w'), indent=1, default=br._js)
print(json.dumps(dict(changes=trans, old_b=old['s_star_option_b'], new_b=new['s_star_option_b'],
  table={r['label']:{p['option']:[round(p['s_star_flip'],3),round(p['s_star_scr'],3),round(p['s_star'],3),[round(z,3) for z in p['ci90_diff']],p['gate_pass']] for p in r['per_sector_option']} for r in (old,new)}), indent=1, default=br._js))
