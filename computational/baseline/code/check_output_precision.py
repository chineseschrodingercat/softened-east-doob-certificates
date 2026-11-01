"""Check that endpoint output does not impose an arbitrary precision floor."""
from pathlib import Path
from fractions import Fraction as F
import json
from verify_certificate import verify
from certified_reference import log_error_bound
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'results'/'precision_check';out.mkdir(exist_ok=True)
candidate=out/'tiny_tilt_constant_graph.json'
data=dict(model='periodic_soft_east_activity',density_convention='U=sum(density); no implicit s factor',
    bit_order='leftmost_is_least_significant',c='0.3',kappa='0.01',s='1e-60',R=2,n_min=3,
    density=['0']*4,z_lower=['0']*4,z_upper=['0']*4)
candidate.write_text(json.dumps(data,indent=2),encoding='utf-8')
coarse=verify(candidate,120,write_result=False,places=35)
fine=verify(candidate,120,write_result=False,places=80)
target=F('1e-55')
assert F(fine['exact_width'])<target<F(coarse['exact_width'])
p=ROOT/'results'/'reference_sensitivity'/'perron_k0.01_s0.001_N18_R7.json'
old=json.loads(p.with_name(p.stem+'_verified.json').read_text())
new=verify(p,80,write_result=False)
assert old['lower_density']==new['lower_density'] and old['upper_density']==new['upper_density']
result=dict(status='PASS',default_35_places=coarse,adaptive_80_places=fine,
    target=str(target),primary_endpoints_unchanged=True,
    scope='Tiny-tilt graph check tests output precision, not the A1 convergence radius.')
refcoarse=log_error_bound(3,'0.3','1','1e-100',[1]*8,dps=140,places=65)
reffine=log_error_bound(3,'0.3','1','1e-100',[1]*8,dps=140,places=110)
assert F(reffine['log_error_upper_exact'])<=F(reffine['log_error_upper'])<F('1e-80')<F(refcoarse['log_error_upper'])
result.update(reference_default=refcoarse,reference_adaptive=reffine,
              reference_exact_bound_preserved=True)
(out/'validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(dict(status='PASS',coarse_width=coarse['width'],fine_width=fine['width'],
    target=float(target),primary_endpoints_unchanged=True),indent=2))
