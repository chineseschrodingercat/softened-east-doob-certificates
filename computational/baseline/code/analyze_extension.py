"""Assemble exact verified records; recheck new candidates at 120 digits."""
from pathlib import Path
from fractions import Fraction as F
import csv,json,hashlib,time,subprocess,sys
from mpmath import iv
from verify_certificate import verify
from constant_baseline import constant_certificate
from combine_certificates import intersection
from practical_locality import I,record,upper
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results'/'reference_sensitivity'; NEW=ROOT/'results'/'extension'
OUT=ROOT/'results'/'revision_analysis';OUT.mkdir(exist_ok=True)
def save(name,x): (OUT/name).write_text(json.dumps(x,indent=2),encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def writecsv(name,rows):
    keys=sorted(set().union(*(x.keys() for x in rows)))
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)
summary=json.loads((NEW/'run_summary.json').read_text())
rows=[];members=[];rechecks=[]
for folder in (OLD,NEW):
    for p in sorted(folder.glob('*_metrics.json')):
        row=json.loads(p.read_text());candidate=folder/row['candidate']
        cert=json.loads(candidate.with_name(candidate.stem+'_verified.json').read_text())
        cert['exact_width']=str(F(cert['upper_density'])-F(cert['lower_density']))
        assert cert['status']=='VERIFIED_ALL_EDGES' and sha(candidate)==cert['candidate_sha256']==row['sha256']
        assert cert['lower_density']==row['lower'] and cert['upper_density']==row['upper']
        if folder==NEW:
            repeated=verify(candidate,120,write_result=False)
            assert (repeated['lower_density'],repeated['upper_density'])==(cert['lower_density'],cert['upper_density'])
            rechecks.append(dict(candidate=candidate.name,status='IDENTICAL_ENDPOINTS_AT_120_DIGITS',seconds=repeated['seconds']))
        row.update(exact_width=cert['exact_width'],source=str(candidate.relative_to(ROOT)).replace('\\','/'))
        rows.append(row)
        members.append(dict(cert,model='periodic_soft_east_activity',member_id=candidate.name,
                            file=row['source'],verification_sha256=sha(candidate.with_name(candidate.stem+'_verified.json'))))
save('precision_rechecks.json',rechecks)
writecsv('all_candidates.csv',rows)
combined=[]
for k,s in sorted({(F(x['kappa']),F(x['s'])) for x in rows}):
    const=constant_certificate('0.3',str(k),str(s),120);const.update(member_id='constant')
    selected=[x for x in members if F(x['kappa'])==k and F(x['s'])==s]
    for R in sorted({x['R'] for x in selected}):
        answer=intersection([const]+[x for x in selected if x['R']<=R])
        answer['maximum_range']=R;answer['baseline_not_worsened']=F(answer['exact_width'])<=F(const['exact_width'])
        assert answer['baseline_not_worsened']
        save(f'combined_k{float(k):g}_s{float(s):g}_through_R{R}.json',answer)
        combined.append({key:answer[key] for key in ['kappa','s','maximum_range','all_integer_n_at_least','lower_density','upper_density','exact_width','width','lower_source','upper_source']})
writecsv('combined.csv',combined)
weak=sorted([x for x in rows if F(x['kappa'])==F('.01') and int(x['R'])==7 and x['method']=='perron'],key=lambda x:int(x['N']))
w18=F(next(x['exact_width'] for x in weak if int(x['N'])==18))
for x in weak:x['delta_width_from_N18_exact']=str(F(x['exact_width'])-w18)
writecsv('reference_zoom.csv',weak)
iv.dps=100;bounds=[]
for row in sorted(rows,key=lambda x:int(x['R'])):
    if not (F(row['kappa'])==1 and F(row['s'])==F('.02') and row['method']=='perron' and int(row['N'])==18):continue
    p=ROOT/row['source'];loc=(NEW if p.parent==NEW else ROOT/'results'/'practical')/(p.stem+'_locality.json')
    d=json.loads(loc.read_text());assert d['status']=='CERTIFIED_INFINITE_AND_ALL_RING_LOG_ERROR'
    R=int(row['R']);q=F('.7');cap=F('.02');a=F('.3');ell0=F(d['local_ell0']['upper_exact'])
    prior=I(2*cap/(q*(1-ell0)))*iv.exp(-I(a*9))+I(2*cap/(q*q))*iv.exp(-I(a*(18-R+1)))/(1-iv.exp(-I(a)))
    bounds.append(dict(R=R,N=18,prior=record(prior),posterior=d['exact_vacuum_reference_coefficient_error'],
                      Ea=d['infinite_Xa_error'],E0=d['infinite_X0_error'],locality_file=str(loc.relative_to(ROOT)).replace('\\','/')))
save('reference_bounds.json',bounds)
d=json.loads((ROOT/'results'/'practical'/'perron_k1_s0.02_N18_R7_locality.json').read_text())
E0=F(d['infinite_X0_error']['upper_exact']);Ea=F(d['infinite_Xa_error']['upper_exact'])
tails=[]
for K in range(7,81):
    bound=min(E0,upper(I(Ea)*iv.exp(-I(F('.3')*K))))/F('.7')
    tails.append(dict(K=K,flip_tail_upper=record(I(bound))))
save('spatial_tail_bounds.json',tails)
save('summary.json',dict(status='VERIFIED_RECORDS_ASSEMBLED',graph_candidates=len(rows),new_candidates=len(rechecks),
     intersections=len(combined),new_graph_precision=120,prior_graph_precision=80,
     constant_count=len({(x['kappa'],x['s']) for x in combined}),
     candidate_failures=summary['failures'],N21_decision='Not added: N17 and N19 resolve the local shape around the certified N18-to-N20 deterioration. No monotonicity or asymptotic shape is inferred.'))
print('Assembled',len(rows),'graph certificates;',len(rechecks),'new 120-digit rechecks;',len(combined),'intersections.')
