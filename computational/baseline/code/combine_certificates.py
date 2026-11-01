"""Intersection of independently reverified certificates, with provenance."""
from pathlib import Path
from fractions import Fraction as F
import argparse,csv,hashlib,json,time
from constant_baseline import constant_certificate
from verify_certificate import verify
ROOT=Path(__file__).resolve().parents[1]

def intersection(members):
    if not members: raise ValueError('No members')
    identity=lambda x:(x['model'],F(x['c']),F(x['kappa']),F(x['s']))
    if any(identity(x)!=identity(members[0]) for x in members):
        raise ValueError('Different models or exact parameters')
    if any(x['status'] not in ('VERIFIED_ALL_EDGES','VERIFIED_CONSTANT_CYCLES') for x in members):
        raise ValueError('Every member must be independently verified')
    low=max(members,key=lambda x:F(x['lower_density']))
    high=min(members,key=lambda x:F(x['upper_density']))
    L,U=F(low['lower_density']),F(high['upper_density'])
    if L>U: raise ArithmeticError('EMPTY_INTERSECTION: inconsistent verified inputs; no interval returned')
    return dict(model=members[0]['model'],c=members[0]['c'],kappa=members[0]['kappa'],s=members[0]['s'],
        status='VERIFIED_INTERSECTION',all_integer_n_at_least=max(x['all_integer_n_at_least'] for x in members),
        lower_density=low['lower_density'],upper_density=high['upper_density'],
        exact_width=str(U-L),width=float(U-L),lower_source=low['member_id'],upper_source=high['member_id'],
        members=members,
        proof='Intersection of enclosures on the common range; no claim of a new convergence theorem.')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--study',type=Path,default=ROOT/'results'/'reference_sensitivity')
    ap.add_argument('--output',type=Path,default=ROOT/'results'/'combined');a=ap.parse_args()
    if a.output.exists() and any(a.output.iterdir()):raise FileExistsError('Use an empty output directory')
    a.output.mkdir(parents=True,exist_ok=True)
    rows=list(csv.DictReader((a.study/'study.csv').open()));started=time.perf_counter();combined=[]
    for k,s in [('1','0.02'),('0.1','0.01'),('0.01','0.001')]:
        const=constant_certificate('0.3',k,s,80)
        const['member_id']='constant';const['verification']='recomputed from the three exact activity cycle means'
        members=[const]
        for R in range(3,8):
            for row in rows:
                if F(row['kappa'])!=F(k) or F(row['s'])!=F(s) or int(row['R'])!=R:continue
                p=a.study/row['candidate'];cert=verify(p,80,write_result=False)
                if cert['candidate_sha256']!=row['sha256']:raise ValueError('Candidate changed since the sensitivity study')
                cert.update(model='periodic_soft_east_activity',member_id=p.name,
                            file='../reference_sensitivity/'+p.name,
                            original_verification_file='../reference_sensitivity/'+p.stem+'_verified.json',
                            original_verification_sha256=hashlib.sha256(p.with_name(p.stem+'_verified.json').read_bytes()).hexdigest())
                members.append(cert)
            answer=intersection(members.copy());answer['maximum_range']=R
            answer['baseline_not_worsened']=F(answer['exact_width'])<=F(const['exact_width'])
            assert answer['baseline_not_worsened']
            assert all(F(answer['exact_width'])<=F(x['upper_density'])-F(x['lower_density']) for x in members)
            path=a.output/f'combined_k{k}_s{s}_through_R{R}.json'
            path.write_text(json.dumps(answer,indent=2),encoding='utf-8')
            combined.append({key:answer[key] for key in ('c','kappa','s','maximum_range','all_integer_n_at_least',
                'lower_density','upper_density','width','lower_source','upper_source','baseline_not_worsened')})
    # A deliberately inconsistent synthetic pair must not produce an interval.
    bad=[dict(const,member_id='synthetic_left',lower_density='0',upper_density='1'),
         dict(const,member_id='synthetic_right',lower_density='2',upper_density='3')]
    try:intersection(bad)
    except ArithmeticError:empty_rejected=True
    else:raise AssertionError('Empty intersection accepted')
    bad[1]=dict(bad[0],kappa='0.2')
    try:intersection(bad)
    except ValueError:mismatch_rejected=True
    else:raise AssertionError('Mismatched model accepted')
    with (a.output/'combined.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(combined[0]));w.writeheader();w.writerows(combined)
    summary=dict(status='ALL_COMBINATIONS_VERIFIED',combinations=len(combined),graph_members_reverified=75,
        constant_members_recomputed=3,empty_intersection_rejected=empty_rejected,
        mismatched_parameters_rejected=mismatch_rejected,seconds=time.perf_counter()-started)
    (a.output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
