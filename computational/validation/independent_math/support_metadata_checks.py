"""Focused checks of n_min guards and native-status metadata.

This focused regression check intentionally calls production verifier entry
points. The separate independent_checks.py still imports no verifier code.
"""
from pathlib import Path
import argparse, copy, hashlib, json, runpy, sys, time

ROOT=Path(__file__).resolve().parents[2]
PACKAGE=ROOT/'package' if (ROOT/'package/external/results').is_dir() else ROOT
OUT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True,help='New output directory; existing directories are refused')
args=parser.parse_args()
output=args.output.resolve()
if output.exists():raise FileExistsError('Output directory already exists: '+str(output))
output.mkdir(parents=True,exist_ok=False)
sys.path.insert(0,str(PACKAGE/'external/code'))
import exact_sdp
from benchmark import verify_anderson

def read(p):return json.loads(p.read_text(encoding='utf8'))
records=[]
for name,required in [('anderson_strong_m2',2),('lti_strong_m4_scs',4),('eisert_moderate_m4_t3',6),('kull_native_moderate_k5_optimized',7),('ffs_t0_weak_m4',4)]:
    p=PACKAGE/'external/results'/name
    spec=read(p/'candidate.json');result=read(p/'result.json')
    verify=verify_anderson if spec['method']=='anderson' else exact_sdp.verify_sdp
    current=verify(spec)
    assert current['upper_density']==result['upper_density']
    assert current['energy_lower_exact']==result['energy_lower_exact']
    rejected=[]
    for bad_value in [1,required-1,str(required),float(required)]:
        bad=copy.deepcopy(spec);bad['n_min']=bad_value
        try:verify(bad)
        except ValueError as error:rejected.append({'value':bad_value,'message':str(error)})
        else:raise AssertionError((name,bad_value,'guard accepted invalid value'))
    records.append({'run':name,'required':required,'valid_endpoint_match':True,'rejected':rejected})

structures=[]
for name,change in [('ffs_t0_weak_m4',{'m':2}),('eisert_moderate_m4_t3',{'overlap':0}),('eisert_moderate_m4_t3',{'overlap':5}),('kull_native_moderate_k5_optimized',{'last_k':2})]:
    spec=read(PACKAGE/'external/results'/name/'candidate.json');spec.update(change)
    try:exact_sdp.verify_sdp(spec)
    except ValueError as error:structures.append({'run':name,'changes':change,'rejected':True,'message':str(error)})
    else:raise AssertionError((name,change,'invalid physical support accepted'))

native=[]
for p in sorted((PACKAGE/'external/results').glob('kull_native_*/result.json')):
    r=read(p);raw=read(p.with_name('matlab_output.json'))
    assert r['yalmip_problem']==raw['problem']
    assert r['solver_status']=='native_problem_'+str(raw['problem'])
    assert r['solver_info']==raw['info']
    verified=exact_sdp.verify_sdp(read(p.with_name('candidate.json')))
    assert verified['upper_density']==r['upper_density']
    assert verified['energy_lower_exact']==r['energy_lower_exact']
    native.append({'run':p.parent.name,'problem':raw['problem'],'label':r['solver_status'],'endpoint_match':True,'startup_and_io_seconds':r['matlab_startup_and_io_seconds']})
assert len(native)==12

# Rerun every independent check, preserving the original recorded output.
independent=runpy.run_path(str(OUT/'independent_checks.py'))
start=time.perf_counter()
for name in ['metadata','selected_certificates','mapping_and_external_endpoints','posterior','monthus']:
    independent[name]()
report=independent['REPORT'];report['elapsed_seconds']=time.perf_counter()-start
report['script_sha256']=hashlib.sha256((OUT/'independent_checks.py').read_bytes()).hexdigest()
(output/'independent_recheck.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')

hashes={str(p.relative_to(PACKAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PACKAGE/'external/code/exact_sdp.py',PACKAGE/'external/code/benchmark.py',PACKAGE/'external/code/run_native_kull.py']}
payload={'status':'PASS','volume_guard_checks':records,'structural_support_checks':structures,'native_metadata_and_certificate_checks':native,'independent_checks_rerun':len(report['checks']),'current_code_sha256':hashes}
(output/'support_metadata_checks.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf8')
print(json.dumps({'status':'PASS','guard_methods':len(records),'native_records':len(native),'independent_checks':len(report['checks']),'startup_min':min(x['startup_and_io_seconds'] for x in native),'startup_max':max(x['startup_and_io_seconds'] for x in native)}),flush=True)
