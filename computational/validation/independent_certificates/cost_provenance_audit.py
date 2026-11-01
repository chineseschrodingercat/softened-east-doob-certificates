"""Independent selected-pair / full-campaign accounting, optional source audit."""
from independent_audit import *
import zipfile,sys,platform

def main():
    configure_output()
    results={p.parent.name:load(p) for p in sorted((EXT/'results').glob('*/result.json'))}
    failures={p.parent.name:load(p) for p in sorted((EXT/'results').glob('*/failure.json'))}
    lower={}
    for name,r in results.items():
        if r['method'] in ['js_block_plus_certificate','variance_plus_certificate']:
            lower[name]={'generation':r['candidate_seconds'],'verification':r['transfer_verification_seconds'],
                         'n18':r['finite_center']['lower'],'limit':r['infinite_center']['Q_infinity']['lower']}
    for p in (EXT/'lower_sources').glob('constant_*.json'):
        r=load(p);c,k,s=map(F,[r['c'],r['kappa'],r['s']]);lo,hi=bounds(fraciv(2*c*(1-c)*(k+c))*(iv.exp(-fraciv(s))-1))
        assert F(r['Q']['lower'])<=lo<=hi<=F(r['Q']['upper'])
        lower[p.stem]={'generation':r['candidate_seconds'],'verification':r['verification_seconds'],'n18':r['Q']['lower'],'limit':r['Q']['lower']}
    pairs=load(EXT/'tables/all_cost_pairings.json');npairs=0
    for p in pairs:
        r=results[p['run']]
        assert p['n_min']<=18
        if r['method'] in ['js_block_plus_certificate','variance_plus_certificate']:continue
        q=lower[p['lower_source']];assert F(p['lower'])<=F(q[p['target']])<F(p['lower'])+F(1,10**35)
        assert F(p['upper'])==min(F(0),F(r['upper_density']))
        reuse=q['generation'] if r['method'] in ['kull_rg','kull_native'] and p['lower_source']==f"js_{r['case']}_R2" else 0
        expected=r['total_seconds']+q['generation']+q['verification']-reuse
        assert abs(expected-p['total_s'])<1e-10
        assert abs(p['candidate_s']+p['verification_s']-p['total_s'])<1e-10
        npairs+=1
    campaign=load(EXT/'tables/campaign_costs.json');failed_seen=[]
    for row in campaign:
        case,method=row['case'],row['method'];runs=[r for r in results.values() if r['case']==case and r['method']==method];f=[]
        for name,r in failures.items():
            args=r.get('args',[]);flags=dict(zip(args[::2],args[1::2]))
            if r.get('case',flags.get('--case'))==case and r.get('method',flags.get('--method'))==method:f.append((name,r));failed_seen.append(name)
        assert row['completed']==len(runs) and row['failed']==len(f)
        assert row['inaccurate']==sum(r.get('solver_status')=='optimal_inaccurate' for r in runs)
        assert row['native_nonzero_problem']==sum(r.get('yalmip_problem',0)!=0 for r in runs)
        assert row['solver_warning_count']==sum(r.get('solver_status')=='optimal_inaccurate' or r.get('yalmip_problem',0)!=0 for r in runs)
        core=sum(r['total_seconds']-r.get('shared_mps_generation_seconds',0) for r in runs)
        failuretime=sum(r.get('wall_seconds',r.get('seconds',0)) for _,r in f)
        shared=0
        if method not in ['js_block_plus_certificate','variance_plus_certificate']:
            shared=sum(q['generation']+q['verification'] for name,q in lower.items() if case in name)
        if method in ['kull_rg','kull_native']:shared+=load(EXT/'maps'/f'optimized_{case}.json')['map_optimization_seconds']
        for key,val in [('successful_core_s',core),('failed_worker_s',failuretime),('shared_preparation_s',shared),('campaign_s',core+failuretime+shared)]:assert abs(row[key]-val)<1e-9,(case,method,key)
    assert sorted(failed_seen)==sorted(failures)
    native=[]
    for name,r in results.items():
        if r['method']!='kull_native':continue
        raw=load(EXT/'results'/name/'matlab_output.json')
        assert r['solver_status']=='native_problem_'+str(raw['problem'])
        assert raw['author_function'].replace('\\','/').endswith('/author_code/RDM_Constraints_Renormalization/spinChains_1D/MPS_Relaxation/Relaxation_dual_SDP.m')
        assert r['matlab_startup_and_io_seconds']==max(0,r['native_process_wall_seconds']-raw['candidate_seconds'])
        native.append({'run':name,'core_seconds':r['total_seconds'],'MATLAB_process_wall_seconds':r['native_process_wall_seconds'],
                       'startup_and_io_seconds':r['matlab_startup_and_io_seconds'],'solver_problem':raw['problem']})
    provenance={'declared_commit':'2e9015fff5d9bc5b170cdc6cee98fbbb928decda','scope':'Optional author source byte check; external code is intentionally not redistributed.'}
    author=ROOT/'author_code/RDM_Constraints_Renormalization';archive=ROOT/'sources/RDM_Constraints_Renormalization.zip'
    if author.is_dir() and archive.is_file():
        expected=load(EXT/'AUTHOR_FUNCTION_SHA256.json');checks=[]
        with zipfile.ZipFile(archive) as z:
            for rel,sha in expected.items():
                b=(author/rel).read_bytes();assert hashlib.sha256(b).hexdigest()==sha
                names=[n for n in z.namelist() if n.endswith('/'+rel)];assert len(names)==1 and z.read(names[0])==b
                checks.append({'source':rel,'sha256':sha})
        provenance['status']='EXACT_MATCH_TO_PINNED_SOURCE_ARCHIVE';provenance['functions']=checks
    else:provenance['status']='NOT_RERUN_SOURCE_ARCHIVE_NOT_SHIPPED'
    save('cost_and_native_provenance.json',{'status':'PASS','external_pairings':npairs,'total_pairings':len(pairs),'campaign_rows':len(campaign),
        'all_failed_attempts_accounted_once':len(failed_seen),'native_runs':native,'author_sources':provenance})
    rows=load(EXT/'tables/all_paired_results.json')
    for row in rows:
        if row['method']=='kull_native':assert row['solver_status']==results[row['run']]['solver_status']
    before=HERE/'exact_external_audit.json';checked=0
    if before.exists():
        for prior in load(before):assert prior['upper']==results[prior['run']]['upper_density'];checked+=1
    before_centers=HERE/'n18_enumeration_audit.json';centers_checked=0
    if before_centers.exists():
        for prior in load(before_centers):
            actual=results[prior['run']]['finite_center']
            assert prior['saved_Q']==[actual['lower'],actual['upper']];centers_checked+=1
    save('f1_repair_recheck.json',{'status':'RESOLVED','native_labels_match_raw_codes':len(native),
         'raw_code_counts':{str(code):sum(r['solver_problem']==code for r in native) for code in [-1,0,5]},
         'main_table_native_status_rows_checked':sum(r['method']=='kull_native' for r in rows),
         'campaign_warning_and_native_code_counters_checked':len(campaign),
         'upper_endpoints_unchanged_against_pre_repair_independent_audit':checked,
         'finite_centers_unchanged_against_pre_repair_independent_audit':centers_checked,
         'evidence':'Comparison to raw MATLAB exports and independent pre-repair numerical outputs; no solver rerun.'})
    files=[EXT/'code'/n for n in ['exact_sdp.py','benchmark.py','summarize.py','run_native_kull.py','native_kull.m']]
    files.extend(EXT/'tables'/n for n in ['all_paired_results.json','all_cost_pairings.json','matched_tolerances.json','campaign_costs.json'])
    snapshots={str(p.relative_to(PACKAGE)).replace('\\','/'):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in files}
    save('input_hashes.json',snapshots)
    save('environment.json',{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),'numpy':np.__version__,
        'mpmath':__import__('mpmath').__version__,'scipy':__import__('scipy').__version__,'interval_dps':iv.dps,
        'warning_seen':'Anaconda emits Failed to find real location of C:\\ProgramData\\anaconda3\\python.exe; both completed audits exit 0 with all assertions evaluated.',
        'package_layout':str(PACKAGE)})
    print(json.dumps({'status':'PASS','external_pairings':npairs,'campaign_rows':len(campaign),'native_runs':len(native),'author_source_status':provenance['status']}))

if __name__=='__main__':main()
