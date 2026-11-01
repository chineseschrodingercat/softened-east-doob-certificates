"""Frozen-certificate replay; optimization and SDP solvers are not required.

All mathematical endpoints are reconstructed. Refuse optimized Python,
overwrites and escaping output paths. The release manifest is checked first.
"""
from pathlib import Path
import argparse,json,time,sys,hashlib
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];PACKAGE=ROOT.parent
sys.path.insert(0,str(PACKAGE/'certificates/code'))
import exact_sdp as ex
from benchmark import verify_anderson
from verify_certificate import verify as graph
from phase0 import certify,rec
from certify_infinite import verify as center
from finite_variational import enclose_Q

def main():
    if not __debug__:raise RuntimeError('Python -O forbidden')
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--case',default='');ap.add_argument('--skip-manifest',action='store_true',help='Development only; release replay checks manifest by default');a=ap.parse_args()
    out=(ROOT/a.output).resolve() if not a.output.is_absolute() else a.output.resolve()
    if not out.is_relative_to(ROOT) or out.exists():raise ValueError('fresh output inside external required')
    if not a.skip_manifest:
        manifest=json.loads((PACKAGE/'MANIFEST.json').read_text())
        for name,item in manifest.items():
            p=(PACKAGE/name).resolve()
            if not p.is_relative_to(PACKAGE) or hashlib.sha256(p.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('manifest mismatch: '+name)
    out.mkdir(parents=True);records=[];st=time.perf_counter()
    for f in sorted((ROOT/'results').glob('*/result.json')):
        if a.case and a.case not in f.parent.name:continue
        data=json.loads(f.read_text());method=data['method'];t=time.perf_counter()
        if method=='anderson':
            result=verify_anderson(json.loads((f.parent/'candidate.json').read_text()))
            assert result['upper_density']==data['upper_density']
        elif method in ['lti','eisert','kull_rg','kull_native','ffs_t0']:
            result=ex.verify_sdp(json.loads((f.parent/'candidate.json').read_text()))
            assert result['upper_density']==data['upper_density']
        elif method in ['js_block_plus_certificate','variance_plus_certificate']:
            trial=json.loads((f.parent/'trial_graph.json').read_text());se=json.loads((f.parent/'seed.json').read_text())
            g=graph(f.parent/'trial_graph.json',100,write_result=False);c=certify(trial,100);q=center(trial,se,100);finite,_=enclose_Q(trial,18,100)
            assert g['lower_density']==data['graph']['lower_density'] and g['upper_density']==data['graph']['upper_density']
            assert c['correction_upper']==data['constants']['correction_upper']
            assert q['Q_infinity']==data['infinite_center']['Q_infinity'];assert rec(*finite)==data['finite_center']
            for row,Q in zip(data['intervals'],[q['Q_infinity'],rec(*finite)]):
                lo=max(F(Q['lower']),F(g['lower_density']));hi=F(g['upper_density'])
                if c['correction_upper'] is not None:hi=min(hi,F(Q['upper'])+F(c['correction_upper']))
                assert row['lower_density']==ex.decimal(lo) and row['upper_density']==ex.decimal(hi,True)
            result={'status':'EXACT_ENDPOINT_MATCH'}
        else:raise ValueError(method)
        records.append({'case':f.parent.name,'status':'PASS','seconds':time.perf_counter()-t})
        print(f.parent.name,'PASS',flush=True)
    table=ROOT/'tables/all_paired_results.json';table_checks=0
    if table.exists():
        for row in json.loads(table.read_text()):
            if a.case and row['case']!=a.case:continue
            d=json.loads((ROOT/'results'/row['run']/'result.json').read_text())
            if d['method'] in ['js_block_plus_certificate','variance_plus_certificate']:
                assert F(json.loads((ROOT/'results'/row['run']/'trial.json').read_text())['s'])>=0
                raw=next(q for q in d['intervals'] if q['target']==row['target'])
                lo=F(raw['lower_density']);hi=F(raw['upper_density'])
            else:
                assert F(json.loads((ROOT/'results'/row['run']/'candidate.json').read_text())['s'])>=0
                common=json.loads((ROOT/'results'/f"js_{d['case']}_R4"/'result.json').read_text())
                Q=common['infinite_center']['Q_infinity'] if row['target']=='limit' else common['finite_center']
                lo=F(ex.decimal(F(Q['lower'])));hi=F(d['upper_density'])
            # s >= 0 in every frozen case, so L_s 1 <= 0 and lambda/n <= 0.
            hi=min(hi,F(0))
            assert F(row['lower'])==lo and F(row['upper'])==hi
            assert F(row['width_exact'])==hi-lo
            table_checks+=1
    from constant_lower import verify as constant_verify
    lower_checks=0
    for p in (ROOT/'lower_sources').glob('constant_*.json'):
        d=json.loads(p.read_text());assert constant_verify(d)==d['Q'];lower_checks+=1
    cost_checks=0
    from summarize import lower_options
    if (ROOT/'tables/all_cost_pairings.json').exists():
        pairs=json.loads((ROOT/'tables/all_cost_pairings.json').read_text())
        for r in pairs:
            if a.case and r['case']!=a.case:continue
            raw=json.loads((ROOT/'results'/r['run']/'result.json').read_text())
            if raw['method'] in ['js_block_plus_certificate','variance_plus_certificate']:
                q=next(q for q in raw['intervals'] if q['target']==r['target'])
                assert F(r['lower'])==F(q['lower_density']) and F(r['upper'])==min(F(q['upper_density']),F(0))
                assert F(r['width_exact'])==F(r['upper'])-F(r['lower']) and r['total_s']==raw['total_seconds']
                continue
            q=next(q for q in lower_options(r['case'],r['target']) if q['source']==r['lower_source'])
            assert F(r['lower'])==F(q['lower']) and F(r['upper'])==min(F(raw['upper_density']),F(0))
            assert F(r['width_exact'])==F(r['upper'])-F(r['lower'])
            reuse=q['candidate_s'] if r['method'] in ['kull_rg','kull_native'] and q['source']==f"js_{r['case']}_R2" else 0
            assert r.get('shared_lower_reuse_s',0)==reuse
            assert abs(r['total_s']-(raw['total_seconds']+q['total_s']-reuse))<1e-9
            cost_checks+=1
        for r in json.loads((ROOT/'tables/matched_tolerances.json').read_text()):
            if a.case and r['case']!=a.case:continue
            eligible=[p for p in pairs if p['case']==r['case'] and p['method']==r['method'] and p['target']==r['target'] and F(p['width_exact'])<=F(r['tolerance'])]
            if r['status']=='ATTAINED':
                best=min(eligible,key=lambda p:p['total_s']);assert r['run']==best['run'] and r['lower_source']==best['lower_source'] and r['total_s']==best['total_s']
            else:assert not eligible
    report={'status':'PASS','certificates':len(records),'paired_table_rows_checked':table_checks,'constant_lower_certificates':lower_checks,'cost_pairings_checked':cost_checks,'elapsed_seconds':time.perf_counter()-st,'checks':records}
    (out/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='checks'}))
if __name__=='__main__':main()
