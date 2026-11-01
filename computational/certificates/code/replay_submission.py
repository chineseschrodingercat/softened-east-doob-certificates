"""Reviewer entry point: verify shipped hashes, recompute certificates, compare exact endpoints."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,subprocess,sys,csv
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rounded(x,up=False,places=30):
 x=F(x)*10**places;k=-((-x.numerator)//x.denominator) if up else x.numerator//x.denominator
 return ('-' if k<0 else '')+str(abs(k)//10**places)+'.'+str(abs(k)%10**places).zfill(places)
def worked(row):
 c,q=row['constants'],row['infinite'];h={}
 for key,v,up in [('alpha_upper',c['alpha']['upper'],1),('beta_lower',c['beta']['lower'],0),('d_upper',c['d']['upper'],1),('Q_lower',q['Q_infinity']['lower'],0),('Q_upper',q['Q_infinity']['upper'],1)]:h[key]=rounded(v,up)
 a,b,d=map(F,[h['alpha_upper'],h['beta_lower'],h['d_upper']]);assert 0<=a<1 and b>0 and d>=0
 h['correction_upper']=rounded(d*d/(4*b*(1-a)**2),True);h['Lambda_lower']=h['Q_lower'];h['Lambda_upper']=rounded(F(h['Q_upper'])+F(h['correction_upper']),True);h['width']=str(F(h['Lambda_upper'])-F(h['Lambda_lower']))
 h.update(candidate=row['candidate'],candidate_sha256=row['sha256'],seed_sha256=q['seed_sha256'],places=30)
 return h
def compare(a,b):
 assert a['sha256']==b['sha256']
 for key in ['lower_density','upper_density']:assert a['graph'][key]==b['graph'][key],(a['candidate'],key)
 if a['quadratic']:
  for key in ['alpha','beta','d','correction_upper']:assert a['constants'][key]==b['constants'][key],(a['candidate'],key)
  for key in ['Q_infinity','spectral_density_limit','spectral_width_upper','seed_sha256']:assert a['infinite'][key]==b['infinite'][key],(a['candidate'],key)
  q=a['infinite'];interval=q['spectral_density_limit'];assert F(interval['upper'])-F(interval['lower'])==F(q['spectral_width_upper'])
  assert F(q['Q_infinity']['upper'])-F(q['Q_infinity']['lower'])==F(1,10**65)
def main():
 if not __debug__:raise RuntimeError('Do not use Python -O or PYTHONOPTIMIZE for certificate validation')
 ap=argparse.ArgumentParser();ap.add_argument('--scope',choices=['example','table1','curve','all'],default='example');ap.add_argument('--output',type=Path,required=True);ap.add_argument('--dps',type=int,default=120);args=ap.parse_args()
 out=(ROOT/args.output).resolve() if not args.output.is_absolute() else args.output.resolve()
 assert out.is_relative_to(ROOT),'Output must stay inside certificates/'
 if out.exists():raise FileExistsError('Choose a fresh output directory')
 package=ROOT.parent;manifest=read(package/'MANIFEST.json')
 for name,item in manifest.items():
  p=(package/name).resolve();assert p.is_relative_to(package.resolve());assert sha(p)==item['sha256'],name
 subprocess.run([sys.executable,'-B',str(ROOT/'code/run_final.py'),'--scope',args.scope,'--dps',str(args.dps),'--output',str(out)],check=True)
 actual=read(out/'summary.json');expected={x['candidate']:x for x in read(ROOT/'published/certificates.json')}
 for row in actual:compare(row,expected[row['candidate']])
 example=next((x for x in actual if x['kappa']=='0.1' and x['N']==18 and x['R']==7),None)
 if example is not None:
  h=worked(example);assert h==read(ROOT/'published/worked_certificate.json');(out/'worked_certificate.json').write_text(json.dumps(h,indent=2)+'\n')
 fields=['candidate','c','kappa','s','N','R','sha256','graph_lower','graph_upper','Q_lower','Q_upper','Lambda_lower','Lambda_upper','alpha_upper','beta_lower','d_upper','correction_upper','seed_sha256']
 with (out/'endpoints.csv').open('w',newline='',encoding='utf-8') as f:
  writer=csv.DictWriter(f,fields);writer.writeheader()
  for x in actual:
   row={k:x[k] for k in fields if k in x};row.update(graph_lower=x['graph']['lower_density'],graph_upper=x['graph']['upper_density'])
   if x['quadratic']:
    c,q=x['constants'],x['infinite'];row.update(Q_lower=q['Q_infinity']['lower'],Q_upper=q['Q_infinity']['upper'],Lambda_lower=q['spectral_density_limit']['lower'],Lambda_upper=q['spectral_density_limit']['upper'],alpha_upper=c['alpha']['upper'],beta_lower=c['beta']['lower'],d_upper=c['d']['upper'],correction_upper=c['correction_upper'],seed_sha256=q['seed_sha256'])
   writer.writerow(row)
 report=dict(status='PASS',scope=args.scope,graph_certificates=len(actual),quadratic_certificates=sum(x['quadratic'] for x in actual),manifest_files_checked=len(manifest),interval_dps=args.dps,endpoint_comparison='exact equality to all published decimal endpoints',scope_note='Graph intervals cover all n>=2R-1; quadratic limit intervals certify Lambda only.')
 (out/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
