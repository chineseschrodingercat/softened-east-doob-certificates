"""Run separately implemented exact matrix and tuple-word certificate checks."""
from pathlib import Path
import argparse,subprocess,sys,json
from package_integrity import MODULE,check_manifest,new_output
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='replay_independent')
    args=ap.parse_args();integrity=check_manifest();dest=new_output(args.output)
    code=Path(__file__).resolve().parent
    calls=[['independent_rational_certificate.py','--output',str(dest/'finite_ring.json')],
           ['independent_simple_witness.py','--output',str(dest/'simple_pair.json')]]
    for R,candidate,lo,hi in [(3,'uniform_R3_simple_candidate.json','-0.351663900893898','-0.349944830273774'),(4,'uniform_R4_candidate.json','-0.350871819626929','-0.350726876487370')]:
        calls.append(['independent_uniform_certificate.py','--candidate',str(MODULE/'results'/candidate),'--output',str(dest/f'uniform_R{R}.json'),'--claimed-lower='+lo,'--claimed-upper='+hi,'--local-record',str(dest/'simple_pair.json')])
    for i,call in enumerate(calls):
        p=subprocess.run([sys.executable,'-B',str(code/call[0]),*call[1:]],capture_output=True,text=True,encoding='utf-8')
        (dest/f'step_{i+1}.txt').write_text(p.stdout+'\n'+p.stderr,encoding='utf-8')
        if p.returncode:raise RuntimeError('Independent check failed: '+call[0]+'; inspect '+str(dest))
    reports={name:json.loads((dest/(name+'.json')).read_text(encoding='utf-8')) for name in ['finite_ring','simple_pair','uniform_R3','uniform_R4']}
    assert all(x['status']=='PASS' for x in reports.values())
    out={'status':'PASS','implementations':'integer matrix and tuple-word full-product checks; no primary verifier imports','checks':{name:reports[name]['checked_counts'] for name in ['uniform_R3','uniform_R4']},'package_integrity':integrity}
    (dest/'REPLAY_REPORT.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,indent=2))
