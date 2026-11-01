"""Separate algebraic and floating diagnostics for the new manuscript claims.

This does not validate the infinite-volume proof by numerical experimentation.
The exact SymPy checks below are finite algebra; the CTMC checks use float64.
"""
from pathlib import Path
import argparse
import itertools
import json
import platform
import sys
import numpy as np
import sympy as sp


def exact_algebra():
    # Nontrivial rational Bernoulli parameter with rational sqrt(c(1-c)).
    c, kappa = sp.Rational(1, 5), sp.Rational(1, 10)
    w = c * (1-c)
    N = sp.diag(0, 1)
    X = sp.Matrix([[c, sp.sqrt(w)], [sp.sqrt(w), 1-c]])
    I = sp.eye(2)
    sx = sp.Matrix([[0, 1], [1, 0]])
    bond = sp.kronecker_product(kappa*I + X, N)
    mid = kappa+sp.Rational(1, 2)
    assert bond == mid*sp.kronecker_product(I, N) + sp.kronecker_product(X-sp.Rational(1, 2)*I, N)
    assert bond == sp.kronecker_product((kappa+c)*I+(1-2*c)*N, N) + sp.sqrt(w)*sp.kronecker_product(sx,N)
    R = -2/(1+kappa)*sp.kronecker_product(I-X,N)
    prev = sp.kronecker_product(N,I)
    comm = R*prev-prev*R
    comm_norm_sq = max((comm.T*comm).eigenvals())
    assert sp.simplify(comm_norm_sq-4*w/(1+kappa)**2) == 0
    ratios = []
    for length in (1, 2, 3, 4):
        def site(A, i):
            return sp.kronecker_product(*[A if j == i else I for j in range(length)])
        HB = kappa*sum((site(N,i) for i in range(length)), sp.zeros(2**length))
        for i in range(1,length):
            HB += site(X,i-1)*site(N,i)
        vac=sp.zeros(2**length,1); vac[0]=1
        first=sp.zeros(2**length,1); first[2**(length-1)]=1
        assert HB*vac == sp.zeros(2**length,1)
        assert HB*first == kappa*first
        PB=vac*vac.T
        assert HB*PB == PB*HB == sp.zeros(2**length)
        ratios.append({'block_length':length,'gap_lower_and_attained':str(kappa)})
    z,q=sp.symbols('z q', positive=True)
    assert sp.simplify(q*z*((z+1/z)/2-1)-q*(z-1)**2/2) == 0
    return {'status':'PASS','method':'exact rational/symbolic identities',
            'c':str(c),'kappa':str(kappa),'commutator_norm_squared':str(comm_norm_sq),
            'block_checks':ratios}


def prefix_graph_check():
    # Elementary excitation transitions toggle j only if j+1 is excited.
    # The exterior right neighbor is allowed to remain permanently excited,
    # which gives the easiest possible emptying, hence an adverse test.
    cases=0
    minima=[]
    for length in (3,4,5,6,7):
        core=range(1,length-1)
        for initial in range(1<<length):
            if not any(initial & (1<<j) for j in core):
                continue
            distance={initial:0}
            queue=[initial]
            for state in queue:
                if state == 0: break
                for j in range(length):
                    if j == length-1 or state & (1<<(j+1)):
                        target=state^(1<<j)
                        if target not in distance:
                            distance[target]=distance[state]+1
                            queue.append(target)
            # The rightmost initially excited core gives the weakest chain count.
            rightmost=max(j for j in core if initial & (1<<j))
            lower=length-rightmost
            assert distance[0]>=lower
            cases+=1
            minima.append(distance[0]-lower)
    return {'status':'PASS','finite_cases':cases,'minimum_margin':min(minima),
            'scope':'finite combinatorial diagnostic, not proof for all block sizes'}


def model(n,c,kappa,s,clock):
    words=list(itertools.product((0,1),repeat=n))
    lookup={x:i for i,x in enumerate(words)}
    S=len(words)
    pi=np.array([np.prod([c if bit else 1-c for bit in x]) for x in words])
    off=np.zeros((S,S)); a0=np.zeros(S)
    flip=np.empty((S,n),dtype=int)
    for row,x in enumerate(words):
        for i in range(n):
            y=list(x);y[i]=1-y[i]
            col=lookup[tuple(y)];flip[row,i]=col
            raw=clock*(kappa+x[(i-1)%n])*(1-c if x[i] else c)
            a0[row]-=raw;off[row,col]=np.exp(-s)*raw
    L=off+np.diag(a0)
    rootpi=np.sqrt(pi)
    H=rootpi[:,None]*L/rootpi[None,:]
    ev,vec=np.linalg.eigh(H)
    psi=vec[:,-1]*np.sign(vec[:,-1].sum())
    assert np.min(psi)>0
    logh=np.log(psi/rootpi)
    return words,pi,off,a0,flip,ev[-1],logh


def check_case(n,c,kappa,s,clock,kind):
    words,pi,off,a0,flip,lam,logh=model(n,c,kappa,s,clock)
    U=np.array([.13*sum(x)-.08*sum(x[i]*x[(i+1)%n] for i in range(n)) for x in words])
    if kind=='constant':U*=0
    if kind=='exact':U=logh.copy()
    nu=pi*np.exp(2*(U-U.max()));nu/=nu.sum()
    mu=pi*np.exp(2*(logh-logh.max()));mu/=mu.sum()
    B=off*np.exp(U[None,:]-U[:,None]);Q=off*np.exp(logh[None,:]-logh[:,None])
    G=B-np.diag(B.sum(1));V=a0+B.sum(1)
    log_ratio=(U-logh)[None,:]-(U-logh)[:,None]
    rate=(B*log_ratio-B+Q).sum(1)
    gap=(lam-nu@V)/n
    kl=nu@rate/n
    quadratic=nu@(Q*np.expm1(log_ratio)**2).sum(1)/(2*n)
    cancellation=nu@(B*log_ratio).sum(1)
    err=max(abs(gap-kl),abs(gap-quadratic),abs(cancellation),np.max(np.abs(nu@G)))
    assert err<2e-10,(n,c,kappa,s,clock,kind,err)
    qmax=Q.max()
    eta=max(abs(log_ratio[row,col]) for row in range(len(pi)) for col in flip[row])
    local_bound=qmax*np.expm1(eta)**2/2
    graph_bound=(V.max()-V.min())/n
    assert gap<=min(local_bound,graph_bound)+2e-10
    # Build the conditional table directly from the stationary probability vector.
    p=np.zeros((len(pi),n))
    for row,x in enumerate(words):
        for i in range(n):
            r0=row if x[i]==0 else flip[row,i]
            r1=flip[r0,i]
            p[row,i]=nu[r1]/(nu[r0]+nu[r1])
    C=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i!=j:C[i,j]=max(abs(p[:,i]-p[flip[:,j],i]))
    alpha=C.sum(axis=0).max()
    beta=min(B[row,col]+B[col,row] for row in range(len(pi)) for col in flip[row])
    d=max(abs(V[row]-V[col]) for row in range(len(pi)) for col in flip[row])
    correction=d*d/(4*beta*(1-alpha)**2) if alpha<1 else None
    if correction is not None:assert gap<=correction+2e-10
    # Nonstationary tests integrate the finite semigroup, including endpoint terms.
    rootnu=np.sqrt(nu)
    e,W=np.linalg.eigh(rootnu[:,None]*G/rootnu[None,:])
    T=.73
    integral=np.empty_like(e)
    small=np.abs(e)<1e-11
    integral[small]=T
    integral[~small]=np.expm1(T*e[~small])/e[~small]
    P=(W*np.exp(T*e))@W.T*rootnu[None,:]/rootnu[:,None]
    integ=(W*integral)@W.T*rootnu[None,:]/rootnu[:,None]
    rho=np.zeros_like(nu);rho[0]=1
    F=U-logh
    nonstation_direct=rho@integ@rate
    nonstation_identity=rho@P@F-rho@F+T*lam-rho@integ@V
    nonstation_error=abs(nonstation_direct-nonstation_identity)
    assert nonstation_error<2e-10
    initial_kl=nu@np.log(nu/mu)
    stationary_total=T*n*kl
    exact_initial_total=initial_kl+stationary_total
    assert exact_initial_total>=-2e-10
    # Reverse orientation is recorded, without incorrectly equating the two KLs.
    reverse_rate=mu@(-Q*log_ratio-Q+B).sum(1)/n
    return {'n':n,'c':c,'kappa':kappa,'s':s,'clock':clock,'candidate':kind,
            'entropy_density':float(kl),'gap':float(gap),'identity_error':float(err),
            'nonstationary_identity_error':float(nonstation_error),
            'nonstationary_stationary_formula_difference':float(nonstation_direct-stationary_total),
            'reverse_stationary_density':float(reverse_rate),'alpha':float(alpha),
            'posterior_bound':None if correction is None else float(correction),
            'initial_kl':float(initial_kl)}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise SystemExit('Refusing to overwrite diagnostic results')
    exact=exact_algebra()
    prefix=prefix_graph_check()
    rows=[]
    for n,c,(kappa,s),sign,clock,kind in itertools.product(
        (3,4,5),(.3,.5),((1.,.02),(.1,.01),(.01,.001)),(-1,1),(1.,3.7),
        ('constant','pair','exact')):
        rows.append(check_case(n,c,kappa,sign*s,clock,kind))
    # Clock scaling is separately checked against matching parameter rows.
    bykey={tuple(r[k] for k in ('n','c','kappa','s','candidate','clock')):r for r in rows}
    scale_errors=[]
    for r in rows:
        if r['clock']==1.:
            other=bykey[tuple(r[k] for k in ('n','c','kappa','s','candidate'))+(3.7,)]
            scale_errors.append(abs(other['entropy_density']-3.7*r['entropy_density']))
    assert max(scale_errors)<2e-10
    assert max(abs(r['nonstationary_stationary_formula_difference']) for r in rows)>1e-5
    result={'status':'PASS','exact_algebra':exact,'prefix_diagnostic':prefix,
        'float_cases':len(rows),'method':'float64 finite-state diagnostics; not interval certificates',
        'max_identity_error':max(r['identity_error'] for r in rows),
        'max_nonstationary_error':max(r['nonstationary_identity_error'] for r in rows),
        'max_clock_scaling_error':max(scale_errors),'python':sys.version,
        'executable':sys.executable,'numpy':np.__version__,'sympy':sp.__version__,
        'platform':platform.platform(),'results':rows}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('results','python','platform')},indent=2))
