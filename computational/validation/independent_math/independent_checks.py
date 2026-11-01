"""Independent implementation checks; imports no production verifier code.

Rational checks reconstruct selected slacks and use Fraction residual sums.
Coefficient enclosures use alternating Taylor series and integer square roots.
Dense spectra below are diagnostics, not the source of any rigorous endpoint.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt, exp, sqrt
import argparse, json, hashlib, time
import numpy as np
from scipy.linalg import eigh

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
PACKAGE = ROOT / 'package' if (ROOT / 'package/external/results').is_dir() else ROOT
RESULTS = PACKAGE / 'external/results'
REPORT = {'description': __doc__, 'checks': []}

def record(label, **kw):
    REPORT['checks'].append(dict(check=label, **kw))
    print(label, 'OK', flush=True)

def read(p): return json.loads(p.read_text(encoding='utf8'))
def mat(x): return np.array([[Q(v) for v in row] for row in x], dtype=object)
def ident(n): return np.eye(n, dtype=int).astype(object)
def tr(a): return sum(a[i, i] for i in range(len(a)))
def expect(a, rho): return sum(a[i,j]*rho[j,i] for i in range(len(a)) for j in range(len(a)))
def bits(i,n): return [(i >> (n-1-j)) & 1 for j in range(n)]
def number(b): return sum(v << (len(b)-1-j) for j,v in enumerate(b))

def insert(op, positions, n):
    # Build from matrix units |output><input| and the identity on the rest.
    positions=list(positions); free=[j for j in range(n) if j not in positions]
    a=np.full((2**n,2**n),Q(0),dtype=object)
    for env in range(2**len(free)):
        eb=bits(env,len(free))
        ids=[]
        for local in range(len(op)):
            b=[0]*n
            for j,x in zip(free,eb):b[j]=x
            for j,x in zip(positions,bits(local,len(positions))):b[j]=x
            ids.append(number(b))
        for j,jj in enumerate(ids):
            for k,kk in enumerate(ids):a[jj,kk]=op[j,k]
    return a

def partial(a, dims, axis):
    rest=dims[:axis]+dims[axis+1:]; r=int(np.prod(rest))
    out=np.full((r,r),Q(0),dtype=object)
    for i in range(r):
        bi=list(np.unravel_index(i,rest))
        for j in range(r):
            bj=list(np.unravel_index(j,rest))
            for x in range(dims[axis]):
                ii=bi[:axis]+[x]+bi[axis:];jj=bj[:axis]+[x]+bj[axis:]
                out[i,j]+=a[np.ravel_multi_index(ii,dims),np.ravel_multi_index(jj,dims)]
    return out

def rational_coefficient_interval(c,s,places=110):
    scale=10**places
    root=isqrt((c*(1-c)).numerator*scale**2//(c*(1-c)).denominator)
    rl,rh=Q(root,scale),Q(root+1,scale)
    term=Q(1); acc=term; elo=ehi=None
    for j in range(1,82):
        term*= -s/j;acc+=term
        if j==80:ehi=acc
        if j==81:elo=acc
    return elo*rl,ehi*rh

def ham(spec,result):
    c,k,s=map(Q,(spec['c'],spec['kappa'],spec['s']))
    a=Q(result['coefficient']['a'])
    al,ah=rational_coefficient_interval(c,s)
    err=(k+1)*max(abs(a-al),abs(a-ah))
    stored=Q(result['coefficient']['local_operator_error'])
    assert stored>=err, 'Stored coefficient enclosure is insufficient'
    h=np.kron(np.diag([k,k+1]),mat([[c,-a],[-a,1-c]]))
    return h,stored

def factor_lower(a,cert):
    assert np.array_equal(a,a.T)
    mu=Q(cert['mu']); scale=Q(cert['factor_scale'])
    b=mat(cert['factor'])/scale
    # Direct rational products, independent of production integer-LCM arithmetic.
    defect=a-mu*ident(len(a))-b@b.T
    gamma=mu-max(sum(abs(x) for x in row) for row in defect)
    assert gamma==Q(cert['lower_exact'])
    return gamma

def inv(a):
    n=len(a); rows=[list(a[i])+[Q(i==j) for j in range(n)] for i in range(n)]
    for col in range(n):
        p=next(j for j in range(col,n) if rows[j][col])
        rows[col],rows[p]=rows[p],rows[col]
        z=rows[col][col];rows[col]=[x/z for x in rows[col]]
        for j in range(n):
            if j==col:continue
            z=rows[j][col];rows[j]=[x-z*y for x,y in zip(rows[j],rows[col])]
    return np.array([row[n:] for row in rows],dtype=object)

def rg_maps(spec,k):
    A=list(map(mat,spec['mps'])); W=mat(spec['gauge']); Wi=inv(W); d=len(A[0]); dd=d*d
    assert np.array_equal(W@Wi,ident(dd))
    columns=[]
    for word in range(2**k):
        p=ident(d)
        for x in bits(word,k):p=p@A[x]
        columns.append(p.reshape(dd))
    V=W@np.array(columns,dtype=object).T
    lc=[];rc=[]
    for x in range(2):
        for j in range(dd):
            E=Wi[:,j].reshape(d,d)
            lc.append(W@(A[x]@E).reshape(dd))
    for j in range(dd):
        for x in range(2):
            E=Wi[:,j].reshape(d,d)
            rc.append(W@(E@A[x]).reshape(dd))
    return V,np.array(lc,dtype=object).T,np.array(rc,dtype=object).T

def reconstruct(spec,h,err):
    kind=spec['method']; dual=list(map(mat,spec['duals']))
    if kind=='eisert':
        m,t=spec['m'],spec['overlap'];yl,yr=dual
        hm=sum((insert(h,[j,j+1],m) for j in range(m-1)))
        sw=2*hm-insert(yl,range(m-t,m),m)-insert(yr,range(t),m)
        ss=insert(h,[t-1,t],2*t)+insert(yl,range(t),2*t)+insert(yr,range(t,2*t),2*t)
        return [sw,ss],[Q(1),Q(1)],Q(2*m-1),err
    if kind=='lti':
        m=spec['m']; y=dual[0]
        return [insert(h,[0,1],m)+insert(y,range(1,m),m)-insert(y,range(m-1),m)],[Q(1)],Q(1),err
    if kind=='ffs_t0':
        m=spec['m'];y=dual[0]
        hh=insert(h,[0,1],m)+insert(h,[1,2],m)
        units=[]
        for i in range(2):
            for j in range(2):
                x=np.full((2,2),Q(0),dtype=object);x[i,j]=Q(1);units.append(insert(x,[1],m))
        commute=[hh@x-x@hh for x in units]
        op=[[x.T@cb for cb in commute] for x in units]
        op=[[(x+x.T)/2 for x in row] for row in op]
        station=[(x+x.T)/2 for x in commute]
        station += [op[i][j]-op[j][i] for i in range(4) for j in range(i)]
        zf=mat(spec['zfactor']);z=zf@zf.T; ys=list(map(Q,spec['stationary_duals']))
        slack=insert(h,[0,1],m)+insert(y,range(1,m),m)-insert(y,range(m-1),m)
        for i in range(4):
            for j in range(4):slack-=z[i,j]*op[i][j]
        for weight,term in zip(ys,station):slack+=weight*term
        correction=err*(4*sum(abs(x) for x in z.flat)+4*sum(abs(x) for x in ys[:4])+8*sum(abs(x) for x in ys[4:]))
        assert correction==Q(spec['coefficient_correction'])
        return [slack],[Q(1)],Q(1),err+correction
    if kind=='kull_rg':
        k0,last=spec['k0'],spec['last_k'];V,L,R=rg_maps(spec,k0)
        y=dual[0];pairs=[dual[j:j+2] for j in range(1,len(dual),2)]
        VI=np.kron(V,ident(2));IV=np.kron(ident(2),V)
        LI=np.kron(L,ident(2));IR=np.kron(ident(2),R)
        s0=insert(h,[0,1],4)+insert(y,range(1,4),4)-insert(y,range(3),4)
        s0+=VI.T@pairs[0][0]@VI+IV.T@pairs[0][1]@IV
        slacks=[s0];bounds=[Q(1)]
        for k in range(k0,last+1):
            yl,yr=pairs[k-k0];sk=-np.kron(ident(2),yl)-np.kron(yr,ident(2))
            if k<last:
                zl,zr=pairs[k-k0+1];sk+=LI.T@zl@LI+IR.T@zr@IR
            Vk,_,_=rg_maps(spec,k)
            Vnext,_,_=rg_maps(spec,k+1)
            assert np.array_equal(Vnext,L@np.kron(ident(2),Vk))
            assert np.array_equal(Vnext,R@np.kron(Vk,ident(2)))
            slacks.append(sk);bounds.append(sum(x*x for x in Vk.flat))
        return slacks,bounds,Q(1),err
    raise ValueError(kind)

def ti_marginal(n,m):
    # Nonproduct, strictly positive, exactly cyclically invariant pure state.
    amp=[]
    for w in range(2**n):
        b=bits(w,n)
        amp.append(1+sum(b)+2*sum(b[j]*b[(j+1)%n] for j in range(n))+3*sum(b[j]*b[(j+2)%n] for j in range(n)))
    arr=np.array(amp,dtype=object).reshape(2**m,2**(n-m)); norm=sum(a*a for a in amp)
    return arr@arr.T/Q(norm)

def check_rg_feasibility(spec,h,slacks,bounds):
    last=spec['last_k']; n=last+2;rho=ti_marginal(n,4);V,L,R=rg_maps(spec,3)
    omega=[]
    for k in range(3,last+1):
        Vk,_,_=rg_maps(spec,k);compression=np.kron(np.kron(ident(2),Vk),ident(2))
        w=compression@ti_marginal(n,k+2)@compression.T
        assert tr(w)>=0 and tr(w)<=bounds[k-2]
        omega.append(w)
    assert np.array_equal(partial(rho,[2,2,2,2],0),partial(rho,[2,2,2,2],3))
    assert np.array_equal(np.kron(V,ident(2))@rho@np.kron(V,ident(2)).T,partial(omega[0],[2,4,2],0))
    assert np.array_equal(np.kron(ident(2),V)@rho@np.kron(ident(2),V).T,partial(omega[0],[2,4,2],2))
    for j in range(len(omega)-1):
        assert np.array_equal(np.kron(L,ident(2))@omega[j]@np.kron(L,ident(2)).T,partial(omega[j+1],[2,4,2],0))
        assert np.array_equal(np.kron(ident(2),R)@omega[j]@np.kron(ident(2),R).T,partial(omega[j+1],[2,4,2],2))
    lhs=expect(insert(h,[0,1],4),rho)
    rhs=expect(slacks[0],rho)+sum(expect(a,w) for a,w in zip(slacks[1:],omega))
    assert lhs==rhs
    return {'physical_ring':n,'compressed_traces':[float(tr(w)) for w in omega], 'trace_bounds':[float(x) for x in bounds[1:]],'duality_identity_exact':True}

def selected_certificates():
    names=['eisert_moderate_m4_t3','lti_strong_m4_scs','ffs_t0_weak_m4','ffs_t0_moderate_m3','kull_rg_weak_k5_optimized_mps','kull_native_moderate_k5_optimized','kull_native_strong_k3_block']
    for name in names:
        t=time.perf_counter();p=RESULTS/name;spec=read(p/'candidate.json');r=read(p/'result.json');h,err=ham(spec,r)
        slacks,bounds,div,corr=reconstruct(spec,h,err)
        gamma=[factor_lower(a,c) for a,c in zip(slacks,spec['psd_certificates'])]
        assert gamma==list(map(Q,r['matrix_lower_exact']))
        assert bounds==list(map(Q,r['trace_bounds']))
        lower=gamma[0]+sum(gamma[j] if spec['method']=='eisert' else min(Q(0),gamma[j])*bounds[j] for j in range(1,len(gamma)))
        lower=lower/div-corr
        assert lower==Q(r['energy_lower_exact'])
        assert Q(r['upper_density'])>=-lower
        extra=check_rg_feasibility(spec,h,slacks,bounds) if spec['method']=='kull_rg' else {}
        record('exact_slack_'+name,matrices=[len(a) for a in slacks],gamma=[float(x) for x in gamma],coefficient_correction=float(corr),upper=r['upper_density'],seconds=time.perf_counter()-t,**extra)

def dense_generator(c,k,s,n):
    size=2**n;L=np.zeros((size,size));p=np.zeros(size)
    for x in range(size):
        b=bits(x,n);p[x]=np.prod([c if a else 1-c for a in b])
        for i in range(n):
            rate=(k+b[(i-1)%n])*(1-c if b[i] else c)
            L[x,x^ (1<<(n-1-i))]+=exp(-s)*rate;L[x,x]-=rate
    H= -np.sqrt(p)[:,None]*L/np.sqrt(p)[None,:]
    return L,H,p

def mapping_and_external_endpoints():
    pars={'strong':(1,.02),'moderate':(.1,.01),'weak':(.01,.001)}
    runs=[(p.parent.name,read(p)) for p in RESULTS.glob('*/result.json')]
    failures=list(RESULTS.glob('*/failure.json'))
    assert len(runs)==114 and len(failures)==10
    for case,(k,s) in pars.items():
        n=9;c=.3;L,H,p=dense_generator(c,k,s,n)
        local=np.kron(np.diag([k,k+1]),[[c,-exp(-s)*sqrt(c*(1-c))],[-exp(-s)*sqrt(c*(1-c)),1-c]])
        assembled=np.asarray(sum(insert(local,[i,(i+1)%n],n) for i in range(n)),float)
        mapping_error=float(np.max(np.abs(H-assembled)));assert mapping_error<1e-14
        eig,vec=eigh(H,subset_by_index=[0,0]);lam=-float(eig[0])/n; margins=[]
        for name,r in runs:
            if r['case']!=case or 'upper_density' not in r:continue
            if r['n_min']>n:continue
            margin=float(r['upper_density'])-lam
            assert margin>=-2e-13,(name,margin)
            margins.append(margin)
        record('independent_generator_'+case,ring=n,mapping_max_error=mapping_error,spectral_density=lam,external_upper_count=len(margins),minimum_margin=min(margins))
    record('run_inventory',success=len(runs),failures=len(failures),native=sum(r['method']=='kull_native' for _,r in runs))

def posterior():
    for case in ['strong','moderate','weak']:
        path=RESULTS/f'js_{case}_R3';spec=read(path/'trial.json');reported=read(path/'result.json')['constants']
        c,k,s=map(float,(spec['c'],spec['kappa'],spec['s']));R=spec['R'];v=np.array(list(map(float,spec['density'])))
        for n in [5,7]:
            L,H,pi=dense_generator(c,k,s,n);sz=2**n;u=np.zeros(sz)
            for x in range(sz):
                b=bits(x,n)
                u[x]=sum(v[sum(b[(j+t)%n]<<t for t in range(R))] for j in range(n))
            h=np.exp(u);rates=np.zeros((n,sz));V=L@h/h;nu=pi*h*h;nu/=nu.sum()
            pc=np.zeros((n,sz))
            for i in range(n):
                mask=1<<(n-1-i)
                for x in range(sz):
                    y=x^mask;rates[i,x]=L[x,y]*h[y]/h[x]
                    lo=x&~mask;hi=lo|mask;pc[i,x]=nu[hi]/(nu[lo]+nu[hi])
            C=np.zeros((n,n))
            for i in range(n):
                for j in range(n):
                    if i!=j:C[i,j]=max(abs(pc[i,x]-pc[i,x^(1<<(n-1-j))]) for x in range(sz))
            alpha=float(C.sum(axis=0).max());beta=min(rates[i,x]+rates[i,x^(1<<(n-1-i))] for i in range(n) for x in range(sz))
            d=max(abs(V[x]-V[x^(1<<(n-1-i))]) for i in range(n) for x in range(sz))
            assert alpha<=float(reported['alpha']['upper'])+3e-14
            assert beta>=float(reported['beta']['lower'])-3e-14
            assert d<=float(reported['d']['upper'])+3e-14
            Qn=float(nu@V/n);ev=float(eigh(H,subset_by_index=[0,0],eigvals_only=True)[0]);lam=-ev/n
            eps=d*d/(4*beta*(1-alpha)**2)
            assert -2e-14<=lam-Qn<=eps+2e-14
            psi=np.sqrt(pi)*h;ray=-float(psi@H@psi/(psi@psi))/n
            assert abs(ray-Qn)<2e-14
            record(f'posterior_{case}_n{n}',alpha=alpha,beta=beta,d=d,Q_n=Qn,spectral_density=lam,actual_gap=lam-Qn,finite_D_bound=eps,stored_uniform_bound=reported['correction_upper'],rayleigh_identity_error=abs(ray-Qn))

def monthus():
    p=list(map(Q,['1.24410','0.80380','1.16128','0.86112','1.23520','0.80959','1.16854','0.85577']))
    z=list(map(Q,['0','-.00220','.29340','.30046','.02592','.03187','.23263','.22102','.00181','-.00068','.28753','.29023','.02373','.02885','.23713','.22918']))
    vals=[]
    for w in range(32):
        flip=w^4;ratio=Q(1)
        for j in range(3):ratio*=p[(flip>>j)&7]/p[(w>>j)&7]
        phi=Q(1+((w>>1)&1),2)*(ratio/2-1)
        vals.append(phi+z[w&15]-z[w>>1])
    maximum=max(vals)
    assert maximum==Q(-141078963529,403146300000)
    assert maximum<Q('-0.3499')
    B=[[0]*8 for _ in range(8)]
    for x in range(8):
        bitsx=[(x>>i)&1 for i in range(3)]
        B[x][x]=12-2*sum(1+a for a in bitsx)
        for i in range(3):B[x][x^(1<<i)]=1+bitsx[(i-1)%3]
    v=[1]*8
    for _ in range(100):v=[sum(a*b for a,b in zip(row,v)) for row in B]
    bv=[sum(a*b for a,b in zip(row,v)) for row in B]
    bounds=[(Q(a,b)-12)/12 for a,b in zip(bv,v)]
    assert min(bounds)>=Q('-0.350706886144410010') and max(bounds)<=Q('-0.350706886144409929')
    a,b,d=map(Q,['.110','.187','.516'])
    assert a+2*b+d==1 and Q('.1434')**2<a*b and Q('.3106')**2<b*d
    assert Q('.1434')+Q('.3106')/2-a-3*b/2-d/2==Q('-.3498')
    record('Monthus_exact_reconstruction',edges=32,max_edge=str(maximum),max_word=[(vals.index(maximum)>>i)&1 for i in range(5)],finite_ring_lower=float(min(bounds)),finite_ring_upper=float(max(bounds)),separation_proven=True)

def metadata():
    candidates=list(RESULTS.glob('*/candidate.json'));counts={}
    for p in candidates:
        x=read(p);method=x['method'];counts[method]=counts.get(method,0)+1
        expected={'anderson':lambda:x['m'],'lti':lambda:x['m'],'eisert':lambda:max(x['m'],2*x['overlap']),'ffs_t0':lambda:x['m'],'kull_rg':lambda:x['last_k']+2}[method]()
        assert x['n_min']>=expected,(str(p),expected,x['n_min'])
        assert Q(x['c'])==Q(3,10)
    record('candidate_scope_metadata',counts=counts,n_min_checked=len(candidates))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True,help='New output directory; existing directories are refused')
    args=parser.parse_args()
    output=args.output.resolve()
    if output.exists():raise FileExistsError('Output directory already exists: '+str(output))
    output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter()
    metadata();selected_certificates();mapping_and_external_endpoints();posterior();monthus()
    REPORT['elapsed_seconds']=time.perf_counter()-start
    REPORT['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (output/'independent_checks.json').write_text(json.dumps(REPORT,indent=2)+'\n',encoding='utf8')
    print('COMPLETED',len(REPORT['checks']),'checks',flush=True)
