"""Generate external candidates and independently verify deterministic endpoints."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1');os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import sys,time,json,argparse,platform,hashlib
ROOT=Path(__file__).resolve().parents[1];PACKAGE=ROOT.parent
sys.path.insert(0,str(PACKAGE/'certificates/code'))
import numpy as np
from scipy.optimize import minimize,linprog
from scipy.sparse import coo_matrix
from scipy.linalg import eig,eigh
from fractions import Fraction as F
import psutil
import exact_sdp as ex
from phase0 import certify,transfer,local_tables,rec
from certify_infinite import seed,verify as verify_center
from finite_variational import enclose_Q
from verify_certificate import verify as verify_graph

CASES={'strong':('1','0.02'),'moderate':('0.1','0.01'),'weak':('0.01','0.001')}
def write(p,d):p.write_text(json.dumps(d,indent=2)+'\n',encoding='utf8')
def data_for(case,R):
    k,s=CASES[case]
    return dict(model='periodic_soft_east_activity',density_convention='U=sum(density); no implicit s factor',
                bit_order='leftmost_is_least_significant',c='0.3',kappa=k,s=s,R=R,n_min=2*R-1,reference_N=18)
def graph_duals(data):
    R=data['R'];V=1<<(2*R-2);words=np.arange(2*V);tail=words&(V-1);head=words>>1
    phi=np.asarray(local_tables(data)[0]);mat=coo_matrix((np.r_[np.ones(2*V),-np.ones(2*V),-np.ones(2*V)],(np.r_[words,words,words],np.r_[tail,head,np.full(2*V,V)])),shape=(2*V,V+1)).tocsr()
    ans=[]
    for sign in [1,-1]:
        obj=np.zeros(V+1);obj[-1]=1
        r=linprog(obj,A_ub=mat,b_ub=-sign*phi,bounds=[(0,0)]+[(None,None)]*V,method='highs',options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
        if not r.success:raise RuntimeError(r.message)
        ans.append(sign*r.x[:V])
    data['z_upper']=[format(float(x),'.17g') for x in ans[0]]
    data['z_lower']=[format(float(x),'.17g') for x in ans[1]]
    return data

def js_trial(case,R,out,objective='energy'):
    t0=time.perf_counter();data=data_for(case,R);optlogs=[]
    # All 2^R block weights, one constant fixed; redundancy is harmless.
    def evalq(z):
        data['density']=[0.0]+list(z)
        try:return -transfer(data)
        except (AssertionError,FloatingPointError):return 1e9
    start=np.zeros((1<<R)-1);jac=False
    if objective=='variance':
        n=18;states=np.arange(1<<n,dtype=np.uint32);idx=np.arange(1<<n);m=2*R-1
        CU=np.zeros((1<<n,1<<R));CV=np.zeros((1<<n,1<<m));ones=np.zeros(1<<n)
        for i in range(n):
            ones+=(states>>i)&1
            for length,C in [(R,CU),(m,CV)]:
                word=np.zeros(len(states),dtype=np.uint32)
                for j in range(length):word|=((states>>((i+j)%n))&1)<<j
                C[idx,word]+=1
        logpi=ones*np.log(.3)+(n-ones)*np.log(.7)
        words=np.arange(1<<m);diff=np.zeros((1<<m,1<<R));k,s=map(float,CASES[case])
        for j in range(R):
            w=(words>>j)&((1<<R)-1);diff[np.arange(1<<m),w^(1<<(R-1-j))]+=1;diff[np.arange(1<<m),w]-=1
        bit=(words>>(R-1))&1;prev=(words>>(R-2))&1;rate=(k+prev)*np.where(bit,.7,.3)
        def evalq(z):
            v=np.r_[0,z];logw=logpi+2*(CU@v);w=np.exp(logw-logw.max());w/=w.sum()
            r=rate*np.exp(diff@v-s);phi=r-rate;V=CV@phi;dV=CV@(r[:,None]*diff)
            mean=w@V;dev=V-mean;variance=w@(dev*dev)
            grad=2*((w*dev)@dV+(w*(dev*dev-variance))@CU)
            return variance/(n*s*s),grad[1:]/(n*s*s)
        jac=True
    sol=minimize(evalq,start,jac=jac,method='BFGS',options={'gtol':1e-10,'maxiter':150 if jac else 600})
    optlogs.append(dict(status=int(sol.status),message=str(sol.message),nit=int(sol.nit),fun=float(sol.fun)))
    # A deterministic second start is used at every parameter, not just adverse ones.
    alt=minimize(evalq,np.full_like(start,-0.04),jac=jac,method='BFGS',options={'gtol':1e-10,'maxiter':150 if jac else 600})
    optlogs.append(dict(status=int(alt.status),message=str(alt.message),nit=int(alt.nit),fun=float(alt.fun)))
    best=min([sol,alt],key=lambda r:r.fun)
    data['density']=['0']+[format(float(x),'.17g') for x in best.x]
    # Normalize additive gauge to avoid needless transfer growth.
    density=np.asarray(data['density'],float);density-=density.mean()
    data['density']=[format(float(x),'.17g') for x in density]
    se=seed(data,digits=85,steps=600)
    candidate_s=time.perf_counter()-t0
    write(out/'trial.json',data);write(out/'seed.json',se)
    t=time.perf_counter();q=verify_center(data,se,100);finite,fd=enclose_Q(data,18,100)
    center_s=time.perf_counter()-t
    t=time.perf_counter();data=graph_duals(data);dual_s=time.perf_counter()-t
    write(out/'trial_graph.json',data)
    t=time.perf_counter();g=verify_graph(out/'trial_graph.json',100,write_result=False);gtime=time.perf_counter()-t
    t=time.perf_counter();const=certify(data,100);ctime=time.perf_counter()-t
    rows=[]
    for target,Q in [('limit',q['Q_infinity']),('n18',rec(*finite))]:
        lo=F(Q['lower']);ghi=F(g['upper_density']);glo=F(g['lower_density'])
        hi=ghi
        if const['correction_upper'] is not None:hi=min(hi,F(Q['upper'])+F(const['correction_upper']))
        if lo>hi:raise ArithmeticError('empty certified intersection')
        rows.append({'target':target,'lower_density':ex.decimal(max(lo,glo)),'upper_density':ex.decimal(hi,True),
                     'width':float(hi-max(lo,glo)),'JS_variational_lower':Q['lower'],
                     'graph_only_width':float(ghi-glo),'graph_upper':g['upper_density'],
                     'quadratic_available':const['correction_upper'] is not None})
    return dict(method='variance_plus_certificate' if objective=='variance' else 'js_block_plus_certificate',case=case,R=R,n_min=2*R-1,candidate_seconds=candidate_s,
                transfer_verification_seconds=center_s,graph_dual_seconds=dual_s,graph_verification_seconds=gtime,
                quadratic_verification_seconds=ctime,verification_seconds=center_s+gtime+ctime,
                total_seconds=candidate_s+center_s+dual_s+gtime+ctime,optimizer_logs=optlogs,
                graph=g,constants=const,infinite_center=q,finite_center=rec(*finite),finite_diagnostics=fd,intervals=rows)

def anderson(case,m,out):
    k,s=CASES[case];t0=time.perf_counter();h,error,meta=ex.local_h('0.3',k,s);H=ex.open_h(h,m);hf=np.asarray(H,float)
    # First spin is fixed under the open Hamiltonian; solve both blocks.
    d=1<<m;v=np.zeros(d)
    for a in range(2):
        inds=slice(a*d//2,(a+1)*d//2);vals,vec=eigh(hf[inds,inds],subset_by_index=[0,0]);v[inds]=np.maximum(np.abs(vec[:,0]),1e-14)
    vi=[F(format(float(x),'.17g')) for x in v];candidate_s=time.perf_counter()-t0
    spec=dict(method='anderson',case=case,c='0.3',kappa=k,s=s,m=m,n_min=m,vector=[str(x) for x in vi])
    write(out/'candidate.json',spec);t=time.perf_counter();verified=verify_anderson(spec);vs=time.perf_counter()-t
    return dict(**spec,**verified,candidate_seconds=candidate_s,verification_seconds=vs,total_seconds=candidate_s+vs,coefficient=meta)

def verify_anderson(spec):
    m=spec['m'];h,e,_=ex.local_h(spec['c'],spec['kappa'],spec['s']);H=ex.open_h(h,m);v=[F(x) for x in spec['vector']]
    if not isinstance(spec['n_min'],int) or spec['n_min']<m:raise ValueError('claimed n_min is below block support')
    if len(v)!=1<<m or min(v)<=0:raise ValueError('positive vector required')
    bound=min(sum(F(H[i,j])*v[j] for j in np.flatnonzero(H[i,:]))/v[i] for i in range(1<<m))
    lower=bound/(m-1)-e
    return dict(status='EXACT_POSITIVE_VECTOR_BLOCK_CERTIFICATE',upper_density=ex.decimal(-lower,True),energy_lower_exact=str(lower))

def solve_sdp(case,method,m,out,overlap=1,last_k=4,solver='CLARABEL',mps_map=None):
    import cvxpy as cp
    k,s=CASES[case];t0=time.perf_counter();h,error,meta=ex.local_h('0.3',k,s);hf=np.asarray(h,float)
    spec=dict(method=method,case=case,c='0.3',kappa=k,s=s,m=m,overlap=overlap,solver=solver)
    equalities=[];constraints=[];mps_generation=0
    if method=='lti':
        rho=cp.Variable((1<<m,1<<m),symmetric=True)
        eq=ex.ptr(rho,m,range(1,m))-ex.ptr(rho,m,range(m-1))==0
        equalities=[eq];constraints=[rho>>0,cp.trace(rho)==1,eq];objective=cp.trace(np.asarray(ex.embed(h,[0,1],m),float)@rho);spec['n_min']=m
    elif method=='eisert':
        t=overlap;rho=cp.Variable((1<<m,1<<m),symmetric=True);sigma=cp.Variable((1<<(2*t),1<<(2*t)),symmetric=True)
        e1=ex.ptr(sigma,2*t,range(t))-ex.ptr(rho,m,range(m-t,m))==0
        e2=ex.ptr(sigma,2*t,range(t,2*t))-ex.ptr(rho,m,range(t))==0
        equalities=[e1,e2];constraints=[rho>>0,sigma>>0,cp.trace(rho)==1,cp.trace(sigma)==1,e1,e2]
        objective=cp.trace(2*np.asarray(ex.open_h(h,m),float)@rho)+cp.trace(np.asarray(ex.embed(h,[t-1,t],2*t),float)@sigma);spec['n_min']=max(m,2*t)
    elif method=='kull_rg':
        k0=3;spec.update(k0=k0,last_k=last_k,n_min=last_k+2)
        # Rational MPS products: precisely the unorthogonalized author-code variant.
        source=ROOT/'results'/f'js_{case}_R2'/'trial.json';tr=json.loads(source.read_text());tr_result=json.loads(source.with_name('result.json').read_text())
        mps_generation=tr_result['candidate_seconds']
        v=list(map(float,tr['density']));c=.3
        # Stored trial words use little endian; matrix physical order is big endian.
        A=[]
        for x in range(2):
            a=np.zeros((2,2),dtype=object)
            for y in range(2):a[x,y]=F(format(np.sqrt(c if x else 1-c)*np.exp(v[x+2*y]),'.12g'))
            A.append(a)
        if mps_map is not None:
            maprecord=json.loads(Path(mps_map).read_text())
            if maprecord['case']!=case:raise ValueError('MPS case mismatch')
            A=[ex.loadmat(a) for a in maprecord['mps']]
            mps_generation=maprecord['total_generation_seconds']
            spec['mps_map_provenance']=maprecord
        spec['mps']=[ex.dumpmat(a) for a in A]
        V,_,_=ex.mps_maps(A,k0);U,sv,_=np.linalg.svd(np.asarray(V,float),full_matrices=False)
        W=np.diag(1/sv)@U.T
        spec['gauge']=ex.dumpmat(np.array([[F(format(float(x),'.14g')) for x in row] for row in W],dtype=object))
        spec['ungauged_singular_values']=sv.tolist()
        V,L,R=ex.gauged_maps(A,k0,spec['gauge']);VxI=np.asarray(np.kron(V,ex.eye(2)),float);IxV=np.asarray(np.kron(ex.eye(2),V),float)
        LxI=np.asarray(np.kron(L,ex.eye(2)),float);IxR=np.asarray(np.kron(ex.eye(2),R),float)
        rho=cp.Variable((16,16),symmetric=True);omega=[cp.Variable((16,16),symmetric=True) for _ in range(k0,last_k+1)]
        constraints=[rho>>0,cp.trace(rho)==1]+[w>>0 for w in omega]
        eq=ex.ptr(rho,4,range(1,4))-ex.ptr(rho,4,range(3))==0;equalities.append(eq)
        equalities.extend([VxI@rho@VxI.T-cp.partial_trace(omega[0],(2,4,2),axis=0)==0,IxV@rho@IxV.T-cp.partial_trace(omega[0],(2,4,2),axis=2)==0])
        for j in range(len(omega)-1):equalities.extend([LxI@omega[j]@LxI.T-cp.partial_trace(omega[j+1],(2,4,2),axis=0)==0,IxR@omega[j]@IxR.T-cp.partial_trace(omega[j+1],(2,4,2),axis=2)==0])
        constraints+=equalities;objective=cp.trace(np.asarray(ex.embed(h,[0,1],4),float)@rho)
    elif method=='ffs_t0':
        if m<3:raise ValueError('FFS requires both neighbors of the perturbed site')
        rho=cp.Variable((1<<m,1<<m),symmetric=True);eq=ex.ptr(rho,m,range(1,m))-ex.ptr(rho,m,range(m-1))==0;equalities=[eq]
        ops=ex.ffs_operators(h,m);C=cp.bmat([[cp.trace(np.asarray(a,float)@rho) for a in row] for row in ops]);positive=C>>0
        station=[cp.trace(np.asarray(a,float)@rho)==0 for a in ex.ffs_stationary(h,m)]
        constraints=[rho>>0,cp.trace(rho)==1,eq,positive]+station
        objective=cp.trace(np.asarray(ex.embed(h,[0,1],m),float)@rho);spec['n_min']=m
    else:raise ValueError(method)
    prob=cp.Problem(cp.Minimize(objective),constraints)
    if solver=='CLARABEL':prob.solve(solver=cp.CLARABEL,tol_gap_abs=1e-10,tol_feas=1e-10,tol_gap_rel=1e-10,max_iter=300)
    else:prob.solve(solver=cp.SCS,eps=1e-10,max_iters=250000)
    if prob.status not in ['optimal','optimal_inaccurate']:raise RuntimeError(str(prob.status))
    spec['duals']=[ex.dumpmat(ex.rat(c.dual_value)) for c in equalities]
    if method=='ffs_t0':
        Z=np.asarray(positive.dual_value);val,U=np.linalg.eigh((Z+Z.T)/2);fac=U@np.diag(np.sqrt(np.maximum(val,0)))
        zf=np.array([[F(format(float(x),'.14g')) for x in row] for row in fac],dtype=object);spec['zfactor']=ex.dumpmat(zf)
        spec['stationary_duals']=[str(F(format(float(q.dual_value),'.14g'))) for q in station]
        z=zf@zf.T;y=list(map(F,spec['stationary_duals']))
        spec['coefficient_correction']=str(error*(4*sum(abs(v) for v in z.flat)+4*sum(abs(v) for v in y[:4])+8*sum(abs(v) for v in y[4:])))
    mats,_,_,_=ex.slacks(spec);spec['psd_certificates']=[ex.make_psd_certificate(a) for a in mats]
    generation=time.perf_counter()-t0;spec['solver_status']=prob.status;spec['raw_primal_energy']=float(prob.value/(2*m-1) if method=='eisert' else prob.value)
    spec['solver_seconds']=prob.solver_stats.solve_time;spec['iterations']=prob.solver_stats.num_iters
    write(out/'candidate.json',spec);t=time.perf_counter();result=ex.verify_sdp(spec);verification=time.perf_counter()-t
    return dict(method=method,case=case,m=m,overlap=overlap,last_k=last_k,**result,
                candidate_seconds=generation+mps_generation,shared_mps_generation_seconds=mps_generation,
                verification_seconds=verification,total_seconds=generation+mps_generation+verification,
                raw_primal_energy=spec['raw_primal_energy'],solver_status=prob.status,
                mps_strategy='general_real_D2_optimized' if mps_map is not None else 'JS_range_two' if method=='kull_rg' else None,
                mps_preparation_peak_bytes=spec.get('mps_map_provenance',{}).get('process_peak_memory_bytes',0),
                solver_seconds=spec['solver_seconds'],roundoff_correction=float(-F(result['upper_density'])-F(str(spec['raw_primal_energy']))),coefficient=meta)

def main():
    if not __debug__:raise RuntimeError('Python -O is not permitted')
    ap=argparse.ArgumentParser();ap.add_argument('--method',required=True);ap.add_argument('--case',choices=CASES,required=True);ap.add_argument('--m',type=int,default=4);ap.add_argument('--overlap',type=int,default=1);ap.add_argument('--last-k',type=int,default=5);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--solver',default='CLARABEL');ap.add_argument('--mps-map',type=Path);a=ap.parse_args()
    out=a.output.resolve()
    if not out.is_relative_to(ROOT.resolve()):raise ValueError('output outside external')
    out.mkdir(parents=True,exist_ok=False);t0=time.perf_counter()
    try:
        if a.method in ['js','variance']:r=js_trial(a.case,a.m,out,objective='variance' if a.method=='variance' else 'energy')
        elif a.method=='anderson':r=anderson(a.case,a.m,out)
        else:r=solve_sdp(a.case,a.method,a.m,out,a.overlap,a.last_k,a.solver,a.mps_map)
        r.update(status=r.get('status','VERIFIED'),process_peak_memory_bytes=psutil.Process().memory_info().peak_wset,
                 wall_seconds=time.perf_counter()-t0,python=sys.version,platform=platform.platform(),
                 memory_scope='Whole fresh worker process peak working set, including imports; excludes parent',
                 timing_scope='One fresh run; core total excludes interpreter/import/output; generation includes solver and certificate-factor construction')
        r['process_peak_memory_bytes']=max(r['process_peak_memory_bytes'],r.get('mps_preparation_peak_bytes',0))
        write(out/'result.json',r);print(json.dumps({k:v for k,v in r.items() if k in ['method','case','R','m','upper_density','total_seconds','intervals','process_peak_memory_bytes','status']}),flush=True)
    except Exception as e:
        import traceback
        r=dict(status='FAILED',method=a.method,case=a.case,m=a.m,overlap=a.overlap,last_k=a.last_k,error=repr(e),traceback=traceback.format_exc(),wall_seconds=time.perf_counter()-t0)
        write(out/'failure.json',r);print(r['traceback']);raise
if __name__=='__main__':main()
