"""Recheck frozen locality certificates, plotted bounds, and intersections.
Read-only with respect to released inputs; save fresh audit records separately.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
from practical_locality import verify_candidate,I
from verify_certificate import endpoint
from mpmath import iv
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results'/'editorial_v14'
def read(p):return json.loads(p.read_text(encoding='utf-8'))
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2),encoding='utf-8')
start=time.perf_counter();records=[]
audit=read(OUT/'graph_rechecks_120.json');byhash={r['sha256']:r for r in audit}
keys=['infinite_X0_error','infinite_Xa_error','flip_log_error','exact_vacuum_reference_coefficient_error']
strong=[r for r in audit if F(r['kappa'])==1];assert len(strong)==32
for i,r in enumerate(strong):
    p=ROOT/r['path'];oldpath=(ROOT/'results'/'practical' if 'reference_sensitivity' in r['path'] else p.parent)/(p.stem+'_locality.json')
    old=read(oldpath);new=verify_candidate(p,dps=120)
    assert new['status']==old['status']=='CERTIFIED_INFINITE_AND_ALL_RING_LOG_ERROR'
    assert new['candidate_sha256']==r['sha256']
    assert new['n_min']==r['n_min']
    for key in keys:assert F(new[key]['upper_exact'])<=F(old[key]['upper_exact']),(p,key)
    records.append(new)
    if (i+1)%8==0: print('Rechecked locality',i+1,flush=True)
save('locality_rechecks_120.json',records)
# Independently evaluate the expressions actually plotted, from the frozen
# exact upper endpoints, using more precision than figure assembly.
iv.dps=120;q=F(7,10);a=F(3,10);cap=F(1,50);gamma=F(13,10);D=F(18,7)
ell0=(I(q)+(I(q)*I(gamma)+I(q*q))*I(D)*(iv.exp(I(D*cap))-1))/I(gamma)
old=read(ROOT/'results'/'practical'/'perron_k1_s0.02_N18_R7_locality.json')
E0=F(old['infinite_X0_error']['upper_exact']);Ea=F(old['infinite_Xa_error']['upper_exact'])
tails=read(ROOT/'results'/'revision_analysis'/'spatial_tail_bounds.json');assert len(tails)==74
for x in tails:
    bound=min(I(E0)/I(q),I(Ea)*iv.exp(-I(a*x['K']))/I(q))
    assert F(x['flip_tail_upper']['upper'])>=endpoint(bound._mpi_[1])
reference_bounds=read(ROOT/'results'/'revision_analysis'/'reference_bounds.json');assert len(reference_bounds)==7
for x in reference_bounds:
    n=18;R=x['R'];loc=read(ROOT/'results'/('practical' if R<=7 else 'extension')/f'perron_k1_s0.02_N18_R{R}_locality.json')
    assert F(x['posterior']['upper'])>=F(loc['exact_vacuum_reference_coefficient_error']['upper_exact'])
    prior=(
      2*I(cap)/(I(q)*(1-ell0))*iv.exp(-I(a*((n+1)//2)))+
      2*I(cap)/(I(q*q)*(1-iv.exp(-I(a))))*iv.exp(-I(a*(n-R+1))))
    assert F(x['prior']['upper'])>=endpoint(prior._mpi_[1])
# Retained reference-vector diagnostics are not re-executed without their caches.
reference_records=[]
for directory in ['reference_certification','extension_reference_checks']:
    reference_records.extend((ROOT/'results'/directory).glob('*_reference_certified.json'))
assert len(reference_records)==17
for p in reference_records:assert read(p)['status']=='VERIFIED_REFERENCE_LOG_ERROR'

combs=[]
for p in sorted((ROOT/'results'/'revision_analysis').glob('combined_*.json')):
    c=read(p);ll=[];uu=[];nn=[]
    for member in c['members']:
        assert all(F(c[k])==F(member[k]) for k in ['c','kappa','s'])
        if 'candidate_sha256' in member:
            r=byhash[member['candidate_sha256']]
            assert F(r['lower'])==F(member['lower_density']) and F(r['upper'])==F(member['upper_density'])
            assert r['n_min']==member['all_integer_n_at_least']
        else:
            from constant_baseline import constant_certificate
            cc=constant_certificate(member['c'],member['kappa'],member['s'],120)
            assert cc['lower_density']==member['lower_density'] and cc['upper_density']==member['upper_density']
        ll.append(F(member['lower_density']));uu.append(F(member['upper_density']));nn.append(member['all_integer_n_at_least'])
    assert max(ll)==F(c['lower_density'])<=F(c['upper_density'])==min(uu)
    assert max(nn)==c['all_integer_n_at_least']
    assert min(uu)-max(ll)==F(c['exact_width'])<=min(u-l for u,l in zip(uu,ll))
    combs.append(dict(file=p.name,members=len(ll),n_min=max(nn),nonempty=True))
assert len(combs)==24
save('intersection_checks.json',combs)
summary=dict(status='PASS',locality=len(records),reference_vectors_rechecked=0,stored_reference_records_retained=len(reference_records),intersections=len(combs),
             tail_curve_points=len(tails),reference_curve_points=len(reference_bounds),precision=120,
             original_upper_bounds_enclose_refined_bounds=True,seconds=time.perf_counter()-start)
save('locality_audit_summary.json',summary);print(json.dumps(summary,indent=2))
