"""Read-only audit of the frozen experiments, independent of paper assembly.

No floating candidate generation is repeated. All graph edges and constants
are verified at 120 digits; original records and timings are never replaced.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_CEILING
import json, hashlib, time, csv
from verify_certificate import verify
from constant_baseline import constant_certificate

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'editorial_v14'
OUT.mkdir(parents=True,exist_ok=True)
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(name,x):(OUT/name).write_text(json.dumps(x,indent=2),encoding='utf-8')
def ceil_sci(x,digits=6):
    with localcontext() as ctx:
        ctx.prec=100
        y=Decimal(x.numerator)/Decimal(x.denominator)
        place=Decimal(10)**(y.adjusted()-digits+1)
        return format(y.quantize(place,rounding=ROUND_CEILING),'.5E').lower()

start=time.perf_counter(); rows=[]; constants=[]; checks=[]
for folder in ['reference_sensitivity','extension']:
    d=ROOT/'results'/folder
    for p in sorted(d.glob('*_metrics.json')):
        m=read(p); candidate=d/m['candidate']; old=read(candidate.with_name(candidate.stem+'_verified.json'))
        new=verify(candidate,120,write_result=False)
        for key in ['candidate_sha256','lower_density','upper_density','all_integer_n_at_least']:
            assert old[key]==new[key],(p,key,old[key],new[key])
        lo,hi=F(new['lower_density']),F(new['upper_density'])
        assert lo==F(m['lower']) and hi==F(m['upper'])
        data=read(candidate)
        row=dict(path=candidate.relative_to(ROOT).as_posix(),method=m['method'],
                 c=str(F(m['c'])),kappa=str(F(m['kappa'])),s=str(F(m['s'])),
                 R=m['R'],N=m['N'],n_min=new['all_integer_n_at_least'],
                 lower=str(lo),upper=str(hi),width=str(hi-lo),sha256=new['candidate_sha256'],
                 seconds_120=new['seconds'],edges=new['edges'],
                 generation_seconds=m['generation_seconds'],graph_seconds=m['graph_seconds'],
                 verification_seconds_80=m['verification_seconds'],
                 generation_peak_mib=m['generation_peak_rss_mb'],
                 verification_process_peak_mib=m['worker_peak_rss_mb'],
                 memory_protocol='separate historical generation and verification processes' if
                 'historical generation' in m['time_provenance'] else 'generation and verification in one process')
        rows.append(row)
        if len(rows)%10==0: print('Verified graphs',len(rows),flush=True)
    for p in sorted(d.glob('constant*_verified.json')):
        old=read(p); new=constant_certificate(old['c'],old['kappa'],old['s'],120)
        assert all(old[k]==new[k] for k in ['lower_density','upper_density','exact_width'])
        constants.append(dict(path=p.relative_to(ROOT).as_posix(),c=str(F(old['c'])),
                              kappa=str(F(old['kappa'])),s=str(F(old['s'])),
                              lower=old['lower_density'],upper=old['upper_density'],width=old['exact_width']))
assert len(rows)==92 and len(constants)==6
save('graph_rechecks_120.json',rows);save('constant_rechecks_120.json',constants)
def get(k,s,R,method='perron',N=18):
    found=[r for r in rows if F(r['kappa'])==F(k) and F(r['s'])==F(s) and r['R']==R and
           r['method']==method and (method!='perron' or r['N']==N)]
    assert len(found)==1
    return found[0]
def base(k,s):return next(r for r in constants if F(r['kappa'])==F(k) and F(r['s'])==F(s))
range_rows=[]
for k,s in [('1','.02'),('.1','.01'),('.01','.001')]:
    crossing=[]
    for R in range(3,10):
        p,v=get(k,s,R),get(k,s,R,'variational')
        if F(p['width'])<F(v['width']):crossing.append(R)
        if R>=7:range_rows.append(dict(kappa=k,s=s,R=R,perron=float(F(p['width'])),
            variational=float(F(v['width'])),ratio=float(F(p['width'])/F(base(k,s)['width']))))
    checks.append(dict(claim='Perron narrower than variational',kappa=k,ranges=crossing))
assert next(x for x in checks if x['kappa']=='.1')['ranges']==[9]
p=get('.01','.001',7);b=base('.01','.001')
assert F(p['lower'])<F(b['lower']) and F(p['upper'])>F(b['upper'])
assert all(F(get('.01','.001',R)['width'])<F(b['width']) for R in (8,9))
ref_rows=[]
for N in [14,16,17,18,19,20]:
    r=get('.01','.001',7,N=N)
    ref_rows.append(dict(N=N,width=float(F(r['width'])),
                         delta_exact=str(F(r['width'])-F(p['width'])),
                         delta_micro=float((F(r['width'])-F(p['width']))*10**6)))
assert min(ref_rows,key=lambda r:r['width'])['N']==17
assert all(ref_rows[i]['width']<ref_rows[i+1]['width'] for i in range(2,5))
cost=[]
for R in (7,8,9):
    rr=[get(k,s,R,m) for k,s in [('1','.02'),('.1','.01'),('.01','.001')]
        for m in ('perron','variational')]
    extrema=lambda f:[min(f(x) for x in rr),max(f(x) for x in rr)]
    cost.append(dict(R=R,generation_graph_seconds=extrema(lambda x:x['generation_seconds']+x['graph_seconds']),
        verification_seconds=extrema(lambda x:x['verification_seconds_80']),
        generation_peak_mib=extrema(lambda x:x['generation_peak_mib']),
        verification_process_peak_mib=extrema(lambda x:x['verification_process_peak_mib']),
        protocol=sorted(set(x['memory_protocol'] for x in rr))))
tilts=[]
for s in ['-.02','-.01','.01','.02']:
    r=get('1',s,7)
    cand=ROOT/r['path']
    if s=='.02': lp=ROOT/'results'/'practical'/(cand.stem+'_locality.json')
    else:lp=cand.with_name(cand.stem+'_locality.json')
    loc=read(lp);x=F(loc['exact_vacuum_reference_coefficient_error']['upper_exact'])
    displayed=f'{float(x):.5e}'; upward=ceil_sci(x)
    tilts.append(dict(s=s,width=float(F(r['width'])),ratio=float(F(r['width'])/F(base('1',s)['width'])),
        bound_exact=str(x),old_display=displayed,old_display_below_upper=F(displayed)<x,
        upward_display=upward,source=lp.relative_to(ROOT).as_posix()))
save('table_checks.json',dict(range_rows=range_rows,reference_rows=ref_rows,cost=cost,tilts=tilts,claims=checks))
summary=dict(status='PASS',graphs=92,constants=6,dps=120,edges=sum(r['edges'] for r in rows),
             all_outward_endpoints_identical=True,seconds=time.perf_counter()-start,
             historical_issues_resolved_in_v13=['Table 5 mixes separate-process and combined-process memory peaks',
                     'Table 4 prints nearest-rounded values in a certified-upper-bound column'])
save('audit_summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
