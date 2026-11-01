"""Build paired endpoints and cost ledger from every retained run."""
from pathlib import Path
from fractions import Fraction as F
import json,csv,sys
from functools import lru_cache
import exact_sdp as ex
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text(encoding='utf8'))
@lru_cache(maxsize=None)
def lower_options(case,target):
    ans=[]
    for p in sorted((ROOT/'results').glob('*/result.json')):
        r=load(p)
        if r['case']!=case or r['method'] not in ['js_block_plus_certificate','variance_plus_certificate']:continue
        Q=r['infinite_center']['Q_infinity'] if target=='limit' else r['finite_center']
        ans.append(dict(source=p.parent.name,lower=ex.decimal(F(Q['lower'])),candidate_s=r['candidate_seconds'],verification_s=r['transfer_verification_seconds'],
            total_s=r['candidate_seconds']+r['transfer_verification_seconds'],peak_MiB=r['process_peak_memory_bytes']/2**20,n_min=r['n_min']))
    p=ROOT/'lower_sources'/f'constant_{case}.json'
    if p.exists():
        r=load(p);ans.append(dict(source=f'constant_{case}',lower=r['Q']['lower'],candidate_s=r['candidate_seconds'],verification_s=r['verification_seconds'],total_s=r['total_seconds'],peak_MiB=r['process_peak_memory_bytes']/2**20,n_min=r['n_min']))
    return ans
def main():
    out=ROOT/'tables';out.mkdir(exist_ok=True);records=[];failures=[]
    for f in sorted((ROOT/'results').glob('*/failure.json')):failures.append(dict(run=f.parent.name,**load(f)))
    for f in sorted((ROOT/'results').glob('*/result.json')):
        r=load(f);case=r['case'];method=r['method'];common=load(ROOT/'results'/f'js_{case}_R4'/'result.json')
        if method in ['js_block_plus_certificate','variance_plus_certificate']:
            for row in r['intervals']:
                upper=ex.decimal(min(F(row['upper_density']),F(0)),True)
                records.append(dict(run=f.parent.name,method=method,case=case,level=f"R={r['R']}",target=row['target'],
                    lower=row['lower_density'],upper=upper,width_exact=str(F(upper)-F(row['lower_density'])),
                    width=float(F(upper)-F(row['lower_density'])),
                    candidate_s=r['candidate_seconds']+r['graph_dual_seconds'],verification_s=r['verification_seconds'],total_s=r['total_seconds'],
                    peak_MiB=r['process_peak_memory_bytes']/2**20,n_min=r['n_min'],solver_status='; '.join(x['message'] for x in r['optimizer_logs']),
                    common_lower_cost_s=0,raw_upper=row['upper_density'],common_lower_run='own trial'))
        else:
            for target,key in [('limit','infinite_center'),('n18','finite_center')]:
                Q=common[key]['Q_infinity'] if target=='limit' else common[key]
                lo=F(Q['lower']);hi=min(F(r['upper_density']),F(0))
                if lo>hi:raise ArithmeticError('empty external interval '+f.parent.name)
                prep=common['candidate_seconds'];ver=common['transfer_verification_seconds']
                level=(f"m={r['m']},t={r['overlap']}" if method=='eisert' else f"D=2,k={r['last_k']}"+(', optimized MPS' if r.get('mps_strategy')=='general_real_D2_optimized' else ', block MPS') if method in ['kull_rg','kull_native'] else f"m={r['m']}")
                records.append(dict(run=f.parent.name,method=method,case=case,level=level,target=target,
                    lower=ex.decimal(lo),upper=ex.decimal(hi,True),width_exact=str(hi-F(ex.decimal(lo))),width=float(hi-F(ex.decimal(lo))),
                    candidate_s=r['candidate_seconds']+prep,verification_s=r['verification_seconds']+ver,
                    total_s=r['total_seconds']+prep+ver,peak_MiB=max(r['process_peak_memory_bytes'],common['process_peak_memory_bytes'])/2**20,
                    n_min=max(r['n_min'],7),solver_status=r.get('solver_status','not applicable'),common_lower_cost_s=prep+ver,
                    raw_upper=r['upper_density'],common_lower_run=f'js_{case}_R4'))
    fields=list(records[0]);
    with (out/'all_paired_results.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(records)
    (out/'all_paired_results.json').write_text(json.dumps(records,indent=2))
    (out/'failures.json').write_text(json.dumps(failures,indent=2))
    selected=[]
    for case in ['strong','moderate','weak']:
        for method in ['anderson','eisert','lti','kull_rg','kull_native','ffs_t0','js_block_plus_certificate','variance_plus_certificate']:
            rows=[r for r in records if r['case']==case and r['method']==method and r['target']=='limit']
            if rows:selected.append(min(rows,key=lambda r:r['width']))
    (out/'best_by_method.json').write_text(json.dumps(selected,indent=2))
    lines=['# All-method comparison','', 'Each external upper bound is paired with the same fresh R=4 Jack-Sollich lower bound. Generation and verification of that lower bound are charged. For every method, the universal bound lambda/n <= 0 for s >= 0 is intersected with the raw endpoint. Selected rows minimize certified width within the retained grid; all runs, both targets and failures remain in all_paired_results.csv. Times are one fresh run; no asymptotic cost claim.','',
           '| Case | Method | Level | Strict lower (outward) | Strict upper (outward) | Width | Gen s | Verify s | Total s | Peak MiB |','|---|---|---|---|---|---|---|---|---|---|---|']
    for r in selected:
        lines.append(f"| {r['case']} | {r['method']} | {r['level']} | {ex.decimal(F(r['lower']),p=15)} | {ex.decimal(F(r['upper']),True,p=15)} | {r['width']:.6g} | {r['candidate_s']:.4g} | {r['verification_s']:.4g} | {r['total_s']:.4g} | {r['peak_MiB']:.1f} |")
    lines+=['',f'Failed runs: {len(failures)}. Paired rows: {len(records)}. Finite n=18 and limit rows are separate in the machine-readable ledger.','', 'For extremely small width, the printed rounded endpoint difference can exceed the width column; the width uses the full 35-decimal endpoints.']
    campaign=[]
    for case in ['strong','moderate','weak']:
        for method in ['anderson','eisert','lti','kull_rg','kull_native','ffs_t0','js_block_plus_certificate','variance_plus_certificate']:
            rs=[load(p) for p in (ROOT/'results').glob('*/result.json') if load(p)['case']==case and load(p)['method']==method]
            fs=[r for r in failures if r.get('case')==case and r.get('method')==method]
            # Timeout arguments still identify the method/case, when present.
            for r in failures:
                args=r.get('args',[])
                if '--case' in args and '--method' in args and args[args.index('--case')+1]==case and args[args.index('--method')+1]==method:fs.append(r)
            if not rs and not fs:continue
            shared=0
            if method not in ['js_block_plus_certificate','variance_plus_certificate']:
                shared=sum(r['total_s'] for r in lower_options(case,'limit'))
            if method in ['kull_rg','kull_native']:
                mp=ROOT/'maps'/f'optimized_{case}.json'
                if mp.exists():shared+=load(mp)['map_optimization_seconds']
            campaign.append(dict(case=case,method=method,completed=len(rs),failed=len(fs),
                inaccurate=sum(r.get('solver_status')=='optimal_inaccurate' for r in rs),
                native_nonzero_problem=sum(r.get('yalmip_problem',0)!=0 for r in rs),
                solver_warning_count=sum(r.get('solver_status')=='optimal_inaccurate' or r.get('yalmip_problem',0)!=0 for r in rs),
                successful_core_s=sum(r['total_seconds']-r.get('shared_mps_generation_seconds',0) for r in rs),
                failed_worker_s=sum(r.get('wall_seconds',r.get('seconds',0)) for r in fs),
                shared_preparation_s=shared,
                campaign_s=sum(r['total_seconds']-r.get('shared_mps_generation_seconds',0) for r in rs)+sum(r.get('wall_seconds',r.get('seconds',0)) for r in fs)+shared,
                max_completed_peak_MiB=max((r['process_peak_memory_bytes']/2**20 for r in rs),default=None),
                timing_scope='All retained successful upper-bound core times plus failed worker elapsed times; eligible lower-pool preparation once for external methods. Imports excluded. Search accounting, not a complexity bound.',
                failure_peak_memory='not recorded for failed workers'))
    (out/'campaign_costs.json').write_text(json.dumps(campaign,indent=2))
    with (out/'campaign_costs.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,list(campaign[0]));w.writeheader();w.writerows(campaign)
    extras=[]
    for p in sorted((ROOT/'results').glob('*/result.json')):
        r=load(p)
        if r['method'] not in ['js_block_plus_certificate','variance_plus_certificate']:continue
        base=r['candidate_seconds']+r['transfer_verification_seconds']
        extras.append(dict(run=p.parent.name,case=r['case'],R=r['R'],candidate_s=r['candidate_seconds'],
            center_validation_s=r['transfer_verification_seconds'],
            additional_upper_certificate_s=r['graph_dual_seconds']+r['graph_verification_seconds']+r['quadratic_verification_seconds'],
            extra_fraction_of_strict_lower_cost=(r['total_seconds']-base)/base,
            graph_only_width=r['intervals'][0]['graph_only_width'],
            final_width=r['intervals'][0]['width'],quadratic_available=r['intervals'][0]['quadratic_available']))
    (out/'trial_certificate_overhead.json').write_text(json.dumps(extras,indent=2))
    tolerance_rows=[]
    methods=['anderson','eisert','lti','kull_rg','kull_native','ffs_t0','js_block_plus_certificate','variance_plus_certificate']
    pairings=[]
    for r in records:
        if r['method'] in ['js_block_plus_certificate','variance_plus_certificate']:
            pairings.append(dict(r,lower_source=r['run']))
            continue
        raw=load(ROOT/'results'/r['run']/'result.json')
        for q in lower_options(r['case'],r['target']):
            lo=F(q['lower']);hi=F(r['upper'])
            if lo>hi:raise ArithmeticError('inconsistent lower/upper pairing')
            reuse=q['candidate_s'] if r['method'] in ['kull_rg','kull_native'] and q['source']==f"js_{r['case']}_R2" else 0
            pairings.append(dict(r,lower=q['lower'],width_exact=str(hi-lo),width=float(hi-lo),
                lower_source=q['source'],candidate_s=raw['candidate_seconds']+q['candidate_s']-reuse,shared_lower_reuse_s=reuse,
                peak_MiB=max(raw['process_peak_memory_bytes']/2**20,q['peak_MiB']),n_min=max(raw['n_min'],q['n_min']),
                verification_s=raw['verification_seconds']+q['verification_s'],total_s=raw['total_seconds']+q['total_s']-reuse))
    (out/'all_cost_pairings.json').write_text(json.dumps(pairings,indent=2))
    for case in ['strong','moderate','weak']:
        for target in ['limit','n18']:
            for tol in ['0.0001','0.000001','0.00000001']:
                for method in methods:
                    eligible=[r for r in pairings if r['case']==case and r['target']==target and r['method']==method and F(r['width_exact'])<=F(tol)]
                    r=min(eligible,key=lambda x:x['total_s']) if eligible else None
                    tolerance_rows.append(dict(case=case,target=target,tolerance=tol,method=method,
                        status='ATTAINED' if r else 'NOT_ATTAINED_IN_RETAINED_GRID',
                        run=r['run'] if r else None,total_s=r['total_s'] if r else None,
                        lower_source=r['lower_source'] if r else None,
                        verification_s=r['verification_s'] if r else None,width=r['width'] if r else None))
    (out/'matched_tolerances.json').write_text(json.dumps(tolerance_rows,indent=2))
    with (out/'matched_tolerances.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,list(tolerance_rows[0]));w.writeheader();w.writerows(tolerance_rows)
    lines+=['','## Cost at a common width tolerance','',
       'Cheapest eligible pair of measured upper/lower computations at each tolerance, selected after the grid for descriptive comparison. External methods may choose any fresh range-2/3/4 block, variance, or constant-trial lower endpoint; they are not forced to pay for range 4. This is not cumulative search cost; see the following ledger. NA means not attained within the retained grid. All three tolerances and both targets are in matched_tolerances.csv.','',
       '| Case | Method | Total s at width 1e-6 | Total s at width 1e-8 |', '|---|---|---|---|']
    for case in ['strong','moderate','weak']:
        for method in methods:
            found=[next(r for r in tolerance_rows if r['case']==case and r['method']==method and r['target']=='limit' and r['tolerance']==t) for t in ['0.000001','0.00000001']]
            values=[f"{r['total_s']:.4g}" if r['total_s'] is not None else 'NA' for r in found]
            lines.append(f"| {case} | {method} | {values[0]} | {values[1]} |")
    lines+=['','## Cumulative campaign costs','',
       'This ledger includes every retained solver level, repair and failure; the eligible lower pool is charged once per external method/case, including map optimization where applicable. Timings exclude interpreter/import startup. Failed-worker peak memory was not recorded and is not inferred from successful runs.','',
       '| Case | Method | Completed certificates | No certificate | Solver warnings | Campaign s |', '|---|---|---|---|---|---|']
    for r in campaign:lines.append(f"| {r['case']} | {r['method']} | {r['completed']} | {r['failed']} | {r['solver_warning_count']} | {r['campaign_s']:.4g} |")
    lines+=['','Solver warnings count CVXPY optimal_inaccurate flags and native nonzero YALMIP problem codes. Native problem 5 or -1 is not called optimal_inaccurate; each raw code/info is retained independently of exact certificate success.']
    (out/'COMPARISON.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    print('\n'.join(lines))
if __name__=='__main__':main()
