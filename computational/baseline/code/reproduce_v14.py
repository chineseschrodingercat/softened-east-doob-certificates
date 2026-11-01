"""Replay frozen certificates into a new directory; never overwrite evidence."""
from pathlib import Path
import argparse,sys,os,shutil,subprocess,json,time
ROOT=Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT.parent/'code'))
from check_manifest import verify_manifest
def main():
 if not __debug__:raise RuntimeError("Do not disable verification assertions")
 manifest_count=verify_manifest(ROOT.parent)
 p=argparse.ArgumentParser()
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--quick',action='store_true',help='Recheck the highlighted kappa=1, s=0.02, R=9 graph only')
 a=p.parse_args();out=a.output.resolve()
 if not out.is_relative_to(ROOT):raise ValueError("Output must remain inside baseline/")
 if out.exists() or out==ROOT or ROOT.is_relative_to(out):
  raise ValueError('Choose a new output directory that does not contain the source.')
 out.mkdir(parents=True)
 if a.quick:
  from verify_certificate import verify
  candidate=ROOT/'results/extension/perron_k1_s0.02_N18_R9.json'
  old=json.loads(candidate.with_name(candidate.stem+'_verified.json').read_text())
  result=verify(candidate,120,write_result=False)
  for key in ['candidate_sha256','lower_density','upper_density','all_integer_n_at_least']:
   assert old[key]==result[key],key
  (out/'highlight_graph.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
  print(json.dumps({'status':'PASS','mode':'quick','n_min':result['all_integer_n_at_least'],'width':result['width'],'edges':result['edges']}))
  return
 for name in ['code','inputs','results']:
  shutil.copytree(ROOT/name,out/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONUTF8='1')
 started=time.perf_counter();steps=[]
 for name in ['audit_graphs_v14.py','audit_locality_v14.py','verify_constants.py','check_output_precision.py']:
  r=subprocess.run([sys.executable,str(out/'code'/name)],cwd=out,env=env,capture_output=True,text=True,encoding='utf-8')
  (out/(name+'.log')).write_text(r.stdout+r.stderr,encoding='utf-8')
  steps.append({'script':name,'returncode':r.returncode})
  if r.returncode:raise RuntimeError(f'{name} failed; inspect its log in the replay directory')
  print(name,'PASS',flush=True)
 report={'status':'PASS','mode':'frozen_certificate_replay','seconds':time.perf_counter()-started,'steps':steps,'manifest_files_checked':manifest_count,'scope':'Frozen graph/locality certificates and plotted bounds; no reference-vector or candidate-generation rerun.'}
 (out/'replay_summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print(json.dumps(report,indent=2))
if __name__=='__main__':main()
