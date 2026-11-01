"""Independent algebra and original-generator implementation checks.

No production certificate module is imported. This script reads frozen inputs
and writes only beside itself. Numerical diagnostics are distinguished from
exact Fraction / directed-interval checks in the output.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json, time, math, hashlib, itertools, argparse
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from scipy.io import loadmat
from mpmath import iv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PACKAGE = ROOT / 'package' if (ROOT / 'package' / 'external').is_dir() else ROOT
EXT = PACKAGE / 'external'
OUTPUT = HERE
CASES = {'strong': ('1','0.02'), 'moderate': ('0.1','0.01'), 'weak': ('0.01','0.001')}
iv.dps = 130

def load(p): return json.loads(p.read_text(encoding='utf-8'))
def save(name,x): (OUTPUT/name).write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def configure_output():
    global OUTPUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='Fresh directory below validation/independent_certificates; preserves frozen numerical evidence.')
    args=parser.parse_args()
    if args.output is not None:
        dest=args.output.resolve()
        if not dest.is_relative_to(HERE):raise ValueError('Diagnostic output must stay below validation/independent_certificates')
        dest.mkdir(parents=True,exist_ok=False);OUTPUT=dest
def fraciv(x):
    x=F(x); return iv.mpf(x.numerator)/x.denominator
def bounds(x):
    return tuple(F((-1 if s else 1)*m)*F(2)**e for s,m,e,b in x._mpi_)
def roundout(x,up=False,p=40):
    x=F(x)*10**p
    q=-((-x.numerator)//x.denominator) if up else x.numerator//x.denominator
    return ('-' if q<0 else '')+str(abs(q)//10**p)+'.'+str(abs(q)%10**p).zfill(p)
def eye(d): return np.array([[F(i==j) for j in range(d)] for i in range(d)],dtype=object)
def mat(x): return np.array([[F(v) for v in row] for row in x],dtype=object)

@lru_cache(None)
def local_coefficient(case):
    k,s=map(F,CASES[case]); c=F('0.3')
    return bounds(iv.exp(-fraciv(s))*iv.sqrt(fraciv(c*(1-c))))

def operator(a,where,n):
    """Independent embedding by checking spectator equality in matrix entries."""
    where=tuple(where); rest=tuple(i for i in range(n) if i not in where)
    local=[]; spectator=[]
    for x in range(2**n):
        bits=tuple((x>>(n-1-i))&1 for i in range(n))
        local.append(sum(bits[j]<<(len(where)-1-t) for t,j in enumerate(where)))
        spectator.append(tuple(bits[j] for j in rest))
    return np.array([[a[local[x],local[y]] if spectator[x]==spectator[y] else F(0)
                      for y in range(2**n)] for x in range(2**n)],dtype=object)

def open_h(h,n):
    out=np.zeros((2**n,2**n),dtype=object)
    for j in range(n-1): out += operator(h,(j,j+1),n)
    return out

def inversion(a):
    """Exact Gauss-Jordan, independent of production gauge implementation."""
    n=len(a); rows=[[F(a[i,j]) for j in range(n)]+[F(i==j) for j in range(n)] for i in range(n)]
    for j in range(n):
        p=next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[p]=rows[p],rows[j]
        z=rows[j][j]; rows[j]=[v/z for v in rows[j]]
        for i in range(n):
            if i==j: continue
            z=rows[i][j]; rows[i]=[a-z*b for a,b in zip(rows[i],rows[j])]
    out=np.array([r[n:] for r in rows],dtype=object)
    assert np.array_equal(a@out,eye(n)) and np.array_equal(out@a,eye(n))
    return out

def maps(spec,k):
    A=list(map(mat,spec['mps'])); d=len(A[0]); W=mat(spec['gauge']) if 'gauge' in spec else eye(d*d)
    Wi=inversion(W)
    # Kronecker vectorization identities: vec(A X)=(A kron I)vec(X)
    # and vec(X A)=(I kron A.T)vec(X), for row-major vectorization.
    left=np.column_stack([np.kron(A[x],eye(d)) for x in range(2)])
    right=np.zeros((d*d,d*d*2),dtype=object)
    for x in range(2): right[:,x::2]=np.kron(eye(d),A[x].T)
    V=np.empty((d*d,2**k),dtype=object)
    for j,w in enumerate(itertools.product(range(2),repeat=k)):
        p=eye(d)
        for x in w: p=p@A[x]
        V[:,j]=p.reshape(-1)
    return W@V, W@left@np.kron(eye(2),Wi), W@right@np.kron(Wi,eye(2))

def reconstruct(spec,result):
    assert spec['case']==result['case']
    assert (F(spec['kappa']),F(spec['s']))==tuple(map(F,CASES[spec['case']]))
    assert F(spec['c'])==F('0.3')
    c,k=F(spec['c']),F(spec['kappa']); a=F(result['coefficient']['a'])
    lo,hi=local_coefficient(spec['case']); err=(k+1)*max(abs(a-lo),abs(a-hi))
    supplied=F(result['coefficient']['local_operator_error'])
    assert supplied>=err
    h=np.kron(np.diag([k,k+1]),mat([[c,-a],[-a,1-c]]))
    method=spec['method']; Y=list(map(mat,spec.get('duals',[])))
    details={}; correction=supplied
    if method=='lti':
        m=spec['m']; sl=[operator(h,(0,1),m)+np.kron(eye(2),Y[0])-np.kron(Y[0],eye(2))]
        traces=[F(1)]; divisor=1; nmin=m
    elif method=='eisert':
        m,t=spec['m'],spec['overlap']; assert 1<=t<=m
        sl=[2*open_h(h,m)-operator(Y[0],range(m-t,m),m)-operator(Y[1],range(t),m),
            operator(h,(t-1,t),2*t)+np.kron(Y[0],eye(2**t))+np.kron(eye(2**t),Y[1])]
        traces=[F(1),F(1)]; divisor=2*m-1; nmin=max(m,2*t)
    elif method=='kull_rg':
        k0,kn=spec['k0'],spec['last_k']; assert k0==3 and len(Y)==1+2*(kn-k0+1)
        V,L,R=maps(spec,k0); alpha=Y[0]; pairs=list(zip(Y[1::2],Y[2::2]))
        Vleft,Vright=np.kron(V,eye(2)),np.kron(eye(2),V)
        sl=[operator(h,(0,1),4)+np.kron(eye(2),alpha)-np.kron(alpha,eye(2))+
            Vleft.T@pairs[0][0]@Vleft+Vright.T@pairs[0][1]@Vright]
        traces=[F(1)]; L2,R2=np.kron(L,eye(2)),np.kron(eye(2),R)
        for j,kk in enumerate(range(k0,kn+1)):
            a,b=pairs[j]; sk=-np.kron(eye(2),a)-np.kron(b,eye(2))
            if kk<kn:
                a,b=pairs[j+1]; sk+=L2.T@a@L2+R2.T@b@R2
            sl.append(sk)
            Vk,_,_=maps(spec,kk); traces.append(sum(F(x*x) for x in Vk.flat))
            Vnext,_,_=maps(spec,kk+1)
            assert np.array_equal(Vnext,L@np.kron(eye(2),Vk))
            assert np.array_equal(Vnext,R@np.kron(Vk,eye(2)))
        details['gauge_inverse_and_composition']='EXACT'
        divisor=1; nmin=kn+2
    elif method=='ffs_t0':
        m=spec['m']; H=open_h(h,m); base=[]
        for i,j in itertools.product(range(2),repeat=2):
            E=np.zeros((2,2),dtype=object); E[i,j]=F(1); base.append(operator(E,(1,),m))
        unsym=[[aa.T@(H@bb-bb@H) for bb in base] for aa in base]
        C=[[(z+z.T)/2 for z in row] for row in unsym]
        B=[(H@a-a@H+(H@a-a@H).T)/2 for a in base]
        B += [C[i][j]-C[j][i] for i in range(4) for j in range(i)]
        zf=mat(spec['zfactor']); Z=zf@zf.T
        yy=list(map(F,spec['stationary_duals'])); assert len(yy)==10
        a=operator(h,(0,1),m)+np.kron(eye(2),Y[0])-np.kron(Y[0],eye(2))
        for i,j in itertools.product(range(4),repeat=2): a-=Z[i,j]*C[i][j]
        for y,b in zip(yy,B): a+=y*b
        sl=[(a+a.T)/2]; traces=[F(1)]; divisor=1; nmin=m
        extra=supplied*(4*sum(abs(z) for z in Z.flat)+4*sum(abs(v) for v in yy[:4])+8*sum(abs(v) for v in yy[4:]))
        assert extra==F(spec['coefficient_correction']); correction+=extra
    elif method=='anderson':
        m=spec['m']; H=open_h(h,m); v=list(map(F,spec['vector']))
        assert len(v)==2**m and min(v)>0
        lower=min(sum(H[i,j]*v[j] for j in range(2**m))/v[i] for i in range(2**m))/(m-1)-supplied
        assert lower==F(result['energy_lower_exact']) and F(result['upper_density'])>=-lower
        return {'method':method,'upper':result['upper_density'],'verified':'EXACT_POSITIVE_VECTOR', 'n_min':m}
    else: raise ValueError(method)
    assert nmin==spec['n_min']==result['n_min']
    mins=[]
    for a,cert in zip(sl,spec['psd_certificates']):
        assert np.array_equal(a,a.T)
        mu=F(cert['mu']); scale=int(cert['factor_scale']); assert scale>0
        B=mat(cert['factor'])/scale
        # Direct Fraction Gram multiplication, unlike integer common-denominator
        # accumulation in production.
        E=a-mu*eye(len(a))-B@B.T
        radius=max(sum(abs(v) for v in row) for row in E)
        val=mu-radius; assert val==F(cert['lower_exact']); mins.append(val)
    assert len(mins)==len(sl)==len(spec['psd_certificates'])
    lower=(sum(mins) if method=='eisert' else mins[0]+sum(min(x,F(0))*t for x,t in zip(mins[1:],traces[1:])))/divisor-correction
    assert list(map(F,result['matrix_lower_exact']))==mins
    assert list(map(F,result['trace_bounds']))==traces
    assert F(result['coefficient_error_exact'])==correction
    assert F(result['energy_lower_exact'])==lower and F(result['upper_density'])>=-lower
    return {'method':method,'upper':result['upper_density'],'verified':'EXACT_RATIONAL_RESIDUAL',
            'slack_dimensions':[len(a) for a in sl], 'n_min':nmin,**details}

def original_cw(n,case):
    c=F('0.3'); k,s=map(F,CASES[case]); N=2**n
    states=np.arange(N,dtype=np.int64); rows=[];cols=[];vals=[];diag=np.zeros(N)
    for i in range(n):
        b=(states>>i)&1; prev=(states>>((i-1)%n))&1
        r=(float(k)+prev)*np.where(b,float(1-c),float(c));diag+=r
        rows.extend(states);cols.extend(states^(1<<i));vals.extend(-(float(k)+prev)*math.exp(-float(s))*math.sqrt(float(c*(1-c))))
    rows.extend(states);cols.extend(states);vals.extend(diag)
    H=coo_matrix((vals,(rows,cols)),shape=(N,N)).tocsr()
    val,vec=eigsh(H,k=1,which='SA',tol=1e-13,v0=np.ones(N))
    x=np.abs(vec[:,0]);ones=np.array([int(j).bit_count() for j in states]);sqrtpi=np.sqrt(float(c)**ones*float(1-c)**(n-ones))
    hv=x/sqrtpi; hv/=max(hv)
    v=[F(format(max(z,1e-30),'.17g')) for z in hv]
    el,eu=bounds(iv.exp(-fraciv(s)));lo=None;hi=None
    for xx in range(N):
        rowlo=F(0);rowhi=F(0)
        for i in range(n):
            b=(xx>>i)&1;prev=(xx>>((i-1)%n))&1;rate=(k+prev)*(1-c if b else c)
            ratio=v[xx^(1<<i)]/v[xx]
            rowlo+=rate*(el*ratio-1);rowhi+=rate*(eu*ratio-1)
        lo=rowlo/n if lo is None else min(lo,rowlo/n)
        hi=rowhi/n if hi is None else max(hi,rowhi/n)
    out={'case':case,'n':n,'dimension':N,'lower':roundout(lo),'upper':roundout(hi,True),
         'width':float(hi-lo),'ritz_diagnostic':float(-val[0]/n),'minimum_positive_rational_vector_entry':str(min(v)),
         'arithmetic':'Original nonsymmetric row generator; rational positive vector; exp at 130 dps outward; all remaining operations Fraction.'}
    save(f'cw_n{n}_{case}.json',dict(out,vector=[str(x) for x in v]))
    return out,lo,hi

def trial_n18(path):
    tr=load(path/'trial.json'); r=load(path/'result.json'); n=18;R=tr['R']; N=2**n;states=np.arange(N,dtype=np.uint32)
    counts=np.zeros((N,2**R),dtype=np.int16)
    for i in range(n):
        ids=sum(((states>>((i+j)%n))&1)<<j for j in range(R))
        counts[np.arange(N),ids]+=1
    one=np.array([int(x).bit_count() for x in states]); logpi=one*math.log(.3)+(n-one)*math.log(.7)
    k,s=float(tr['kappa']),float(tr['s']);v=np.array(tr['density'],dtype=float)
    def evaluate(v,gradient=False):
        U=counts@v; logw=logpi+2*U;w=np.exp(logw-logw.max());w/=w.sum()
        V=np.zeros(N);DV=np.zeros((N,len(v))) if gradient else None
        for i in range(n):
            y=states^(1<<i);rate=(k+((states>>((i-1)%n))&1))*np.where((states>>i)&1,.7,.3)
            tilted=rate*np.exp(U[y]-U-s);V+=tilted-rate
            if gradient:DV+=tilted[:,None]*(counts[y]-counts)
        Q=float(w@V/n);dev=V-float(w@V);var=float(w@(dev*dev));out=[Q,var/(n*s*s)]
        if gradient:
            # Different derivation: sum d(w_x*V_x^2)-2 E(V) dE(V).
            meanC=w@counts; dw=2*w[:,None]*(counts-meanC)
            dmean=dw.T@V+w@DV
            dsecond=dw.T@(V*V)+2*(w*V)@DV
            out.append((dsecond-2*(w@V)*dmean)/(n*s*s))
        return out
    Q,var=evaluate(v)
    qlo,qhi=map(float,(r['finite_center']['lower'],r['finite_center']['upper']))
    qerr=min(abs(Q-qlo),abs(Q-qhi)); assert qerr<2e-14
    out={'run':path.name,'n':n,'range':R,'enumerated_Q':Q,'saved_Q':[r['finite_center']['lower'],r['finite_center']['upper']],
         'double_diagnostic_distance':qerr,'normalized_variance':var,'scope':'Independent direct 2^18-state double enumeration; diagnostic, not an interval proof.'}
    if r['method']=='variance_plus_certificate':
        best=min(q['fun'] for q in r['optimizer_logs']); assert abs(best-var)<max(1e-12,abs(best)*1e-8)
        # Away from stationary point to avoid falsely easy zero-gradient check.
        probe=v+np.array([.003*math.sin(j+.4) for j in range(len(v))]);a=evaluate(probe,True)[2];fd=[]
        step=2e-6
        for j in range(len(v)):
            e=np.zeros(len(v));e[j]=step; fd.append((evaluate(probe+e)[1]-evaluate(probe-e)[1])/(2*step))
        diff=float(np.max(np.abs(a-fd)));reldiff=diff/max(1,float(np.max(np.abs(a))));assert reldiff<2e-6
        out.update(saved_optimizer_variance=best,gradient_probe_max_absolute_error=diff,gradient_probe_relative_error=reldiff,gradient_analytic=a.tolist(),gradient_central_difference=fd)
    return out

def ledger_audit():
    results={p.parent.name:load(p) for p in sorted((EXT/'results').glob('*/result.json'))}
    failures={p.parent.name:load(p) for p in sorted((EXT/'results').glob('*/failure.json'))}
    assert len(results)==114 and len(failures)==10 and not (results.keys()&failures.keys())
    rows=load(EXT/'tables/all_paired_results.json'); assert len(rows)==2*len(results)
    for q in rows:
        r=results[q['run']];assert q['target'] in ['limit','n18']
        assert F(q['width_exact'])==F(q['upper'])-F(q['lower'])>=0
        assert math.isclose(q['candidate_s']+q['verification_s'],q['total_s'],abs_tol=1e-10)
        if r['method'] in ['js_block_plus_certificate','variance_plus_certificate']:
            rr=next(x for x in r['intervals'] if x['target']==q['target']);assert F(q['lower'])==F(rr['lower_density']);assert F(q['upper'])==min(F(rr['upper_density']),F(0))
        else:
            common=results[f"js_{q['case']}_R4"];qc=common['finite_center'] if q['target']=='n18' else common['infinite_center']['Q_infinity']
            assert F(q['lower'])<=F(qc['lower']) and F(qc['lower'])-F(q['lower'])<F(1,10**35)
            assert F(q['upper'])==min(F(r['upper_density']),F(0))
            assert abs(q['total_s']-r['total_seconds']-common['candidate_seconds']-common['transfer_verification_seconds'])<1e-10
            assert q['n_min']==max(r['n_min'],7)
    pairings=load(EXT/'tables/all_cost_pairings.json');matched=load(EXT/'tables/matched_tolerances.json')
    for m in matched:
        eligible=[p for p in pairings if all(m[k]==p[k] for k in ['case','method','target']) and F(p['width_exact'])<=F(m['tolerance'])]
        assert bool(eligible)==(m['status']=='ATTAINED')
        if eligible:
            p=min(eligible,key=lambda p:p['total_s']);assert p['run']==m['run'] and p['lower_source']==m['lower_source'] and p['total_s']==m['total_s']
    for q in load(EXT/'tables/best_by_method.json'):
        rs=[r for r in rows if all(r[k]==q[k] for k in ['case','method','target'])];assert q['width']==min(r['width'] for r in rs)
    failtab=load(EXT/'tables/failures.json');assert {q['run'] for q in failtab}==set(failures)
    counts={}
    for r in results.values():counts[r['method']]=counts.get(r['method'],0)+1
    return {'successful_runs':len(results),'failed_runs':len(failures),'paired_rows':len(rows),'cost_pairings':len(pairings),'matched_tolerance_rows':len(matched),'methods':counts},results

def main():
    configure_output()
    t=time.perf_counter();summary,results=ledger_audit();save('ledger_audit.json',summary);print('ledger',summary,flush=True)
    exact=[]
    for run,r in results.items():
        p=EXT/'results'/run/'candidate.json'
        if not p.exists(): continue
        q=reconstruct(load(p),r);q['run']=run;exact.append(q)
        if r['method']=='kull_native':
            raw=load(p.with_name('matlab_output.json'));sp=load(p)
            assert len(raw['duals'])==len(sp['duals'])
            for a,b in zip(raw['duals'],sp['duals']):
                x=np.asarray(a,float); expected=np.rint((x+x.T)/2*10**14)/10**14
                assert np.max(np.abs(expected-np.asarray(mat(b),float)))<1e-15
            inp=loadmat(p.with_name('input.mat')) if p.with_name('input.mat').exists() else None
            q['native_export_dual_count']=len(raw['duals']);q['native_solver_problem']=raw['problem']
            if inp is not None:
                V,L,R=maps(sp,3)
                for key,a in [('V',V),('L',L),('R',R)]:assert np.array_equal(np.asarray(a,float),inp[key])
        print('exact',run,flush=True)
    save('exact_external_audit.json',exact)
    checks=[]
    for case in CASES:
        cw,lo,hi=original_cw(10,case);valid=[q for q in exact if results[q['run']]['case']==case and q['n_min']<=10]
        for q in valid:assert F(q['upper'])>=hi, (q['run'],q['upper'],cw)
        cw['external_uppers_checked']=len(valid);cw['minimum_upper_above_cw_upper']=float(min(F(q['upper'])-hi for q in valid));checks.append(cw);print('CW',cw,flush=True)
    save('original_generator_crosschecks.json',checks)
    trials=[]
    for run,r in results.items():
        if r['method'] not in ['js_block_plus_certificate','variance_plus_certificate']:continue
        q=trial_n18(EXT/'results'/run);trials.append(q);save('n18_enumeration_audit.json',trials);print('n18',run,q['double_diagnostic_distance'],flush=True)
    summary['exact_external_certificates']=len(exact);summary['independent_n18_trials']=len(trials);summary['elapsed_seconds']=time.perf_counter()-t
    save('audit_summary.json',summary);print('COMPLETE',summary,flush=True)

if __name__=='__main__':main()
