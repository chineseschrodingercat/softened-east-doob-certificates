"""Recompute local rigorous certificates from frozen decimal candidates and seeds."""
import argparse,hashlib,json,time,platform,statistics,sys
from pathlib import Path
from fractions import Fraction as F
import mpmath
from candidate_inputs import candidate_path
from phase0 import certify,dec,rec
from certify_infinite import verify as center,seed
from verify_certificate import verify as graph
ROOT=Path(__file__).resolve().parents[1]
def dump(p,o):p.write_text(json.dumps(o,indent=2)+'\n',encoding='utf-8')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--dps',type=int,default=120);ap.add_argument('--repeats',type=int,default=1);ap.add_argument('--curve-only',action='store_true');ap.add_argument('--scope',choices=['all','curve','table1','example'],default='all');ap.add_argument('--generate-seeds',action='store_true');a=ap.parse_args()
 out=a.output.resolve();assert out.is_relative_to(ROOT)
 if out.exists() and any(out.iterdir()):raise FileExistsError(out)
 out.mkdir(parents=True,exist_ok=True);rows=[]
 cases=json.loads((ROOT/'cases.json').read_text(encoding='utf-8'))
 for case in cases:
  if (a.curve_only or a.scope=='curve') and not case['range_curve']:continue
  if a.scope=='table1' and not (case['N']==18 and (case['R'] in [7,9] or (case['kappa']=='0.1' and case['R']==8))):continue
  if a.scope=='example' and not (case['N']==18 and case['R']==7 and case['kappa']=='0.1'):continue
  path=candidate_path(case['candidate']);raw=path.read_bytes();sha=hashlib.sha256(raw).hexdigest();assert sha==case['sha256'];data=json.loads(raw)
  sf=ROOT/'seeds'/(path.stem+'_seed.json')
  if case['quadratic']:
   if not sf.exists():
    assert a.generate_seeds,'Missing frozen seed';se=seed(data);se['candidate_sha256']=sha;dump(sf,se)
   seeds=json.loads(sf.read_text());assert seeds['candidate_sha256']==sha
  reps=[]
  for rep in range(a.repeats):
   g=graph(path,a.dps,write_result=False);row=dict(case,graph=g)
   if case['quadratic']:
    qt=time.perf_counter();data=json.loads(path.read_bytes());seeds=json.loads(sf.read_bytes())
    c=certify(data,a.dps);assert c['status']=='INTERVAL_CONSTANTS_VERIFIED'
    c.update(candidate=path.name,candidate_sha256=sha)
    q=center(data,seeds,a.dps);quadratic_elapsed=time.perf_counter()-qt;q.update(candidate=path.name,candidate_sha256=sha,seed_sha256=hashlib.sha256(sf.read_bytes()).hexdigest())
    sl=F(q['Q_infinity']['lower']);su=F(q['Q_infinity']['upper'])+F(c['correction_upper']);interval=rec(sl,su)
    assert max(sl,F(g['lower_density']))<=min(su,F(g['upper_density']))
    q.update(correction_upper=c['correction_upper'],spectral_density_limit=interval,spectral_width_upper=dec(F(interval['upper'])-F(interval['lower']),True))
    row.update(constants=c,infinite=q,quadratic_seconds=quadratic_elapsed)
   reps.append(row)
   print(json.dumps({'case':path.stem,'dps':a.dps,'rep':rep+1,'graph_s':g['seconds'],'quadratic_s':quadratic_elapsed if case['quadratic'] else None,'quadratic_width':q['spectral_width_upper'] if case['quadratic'] else None}),flush=True)
  row=reps[0];row['timing']={'repetitions':a.repeats,'graph_seconds':[v['graph']['seconds'] for v in reps],'scope':'verification including reading/parsing the frozen inputs; excludes candidate/seed generation, imports and output I/O; no simultaneous numerical job'}
  # Both timers include reading/parsing the relevant frozen inputs.
  if case['quadratic']:row['timing']['quadratic_seconds']=[v['quadratic_seconds'] for v in reps]
  dump(out/(path.stem+'_certificate.json'),row);rows.append(row)
 dump(out/'summary.json',rows);dump(out/'runtime.json',dict(python=sys.version,platform=platform.platform(),mpmath=mpmath.__version__,dps=a.dps,repetitions=a.repeats,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))
if __name__=='__main__':main()
