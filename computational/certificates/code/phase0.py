"""Independent Phase 0 for quadratic variational corrections.

Inputs are exact stored decimals. Interval certification covers alpha, beta,
physical TOTAL-residual oscillation, and the resulting correction only.
The transfer-matrix centers in this first experiment are floating diagnostics.
No archived module is imported and no archived result is overwritten.
"""
import argparse
import csv
import hashlib
import json
import platform
import time
from fractions import Fraction as F
from pathlib import Path

import mpmath
from mpmath import iv
import numpy as np
from scipy.linalg import eig
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[1]
from candidate_inputs import candidate_paths


def endpoint(t):
    sign, mantissa, exponent, bits = t
    if bits < 0:
        raise ArithmeticError('Nonfinite interval')
    x = F((-1 if sign else 1) * mantissa)
    return x * 2**exponent if exponent >= 0 else x / 2**(-exponent)


def I(x):
    x = F(x)
    return iv.mpf(x.numerator) / x.denominator


def bounds(x):
    return endpoint(x._mpi_[0]), endpoint(x._mpi_[1])


def dec(x, up=False, places=65):
    t = F(x) * 10**places
    q = -((-t.numerator)//t.denominator) if up else t.numerator//t.denominator
    sign = '-' if q < 0 else ''
    q = abs(q)
    return f'{sign}{q//10**places}.{q%10**places:0{places}d}'


def rec(lo, hi=None):
    hi = lo if hi is None else hi
    return {'lower':dec(lo), 'upper':dec(hi,True)}


def local_tables(data, interval=False):
    R = data['R']; m = 2*R-1
    if interval:
        c,k,s = [F(data[t]) for t in ('c','kappa','s')]
        v = list(map(F,data['density']))
    else:
        c,k,s = [float(F(data[t])) for t in ('c','kappa','s')]
        v = list(map(float,data['density']))
    phi=[]; b=[]; p=[]
    for word in range(1<<m):
        delta=sum(v[((word>>j)&((1<<R)-1))^(1<<(R-1-j))]-v[(word>>j)&((1<<R)-1)] for j in range(R))
        bit=(word>>(R-1))&1; prev=(word>>(R-2))&1
        rate=(k+prev)*(c if bit==0 else 1-c)
        if interval:
            z=I(rate)*iv.exp(I(delta-s))
            phi.append(bounds(z-I(rate))); b.append(bounds(z))
            # Conditional log odds use U(x_i=1)-U(x_i=0).
            odds=I(c)*iv.exp(I(2*delta*(1-2*bit)))
            p.append(bounds(odds/(I(1-c)+odds)))
        else:
            z=rate*np.exp(delta-s)
            phi.append(rate*np.expm1(delta-s)); b.append(z)
            p.append(1/(1+((1-c)/c)*np.exp(-2*delta*(1-2*bit))))
    return phi,b,p


def maxplus(phi, m, interval=False):
    """Exact-rational DP if interval=True; unrestricted long-window maximum.

    Valid as an upper bound on every ring n>=m. Equal to the physical maximum
    for n>=2m-1. Short rings have extra identifications between noncentral bits.
    """
    V=1<<(m-1); mask=V-1
    if interval:
        dlo=[F(0)]*V; dhi=[F(0)]*V
        for t in range(m):
            flip=1<<(m-1-t)
            nlo=[]; nhi=[]
            for head in range(V):
                w0=head<<1; w1=w0|1
                a=w0&mask; b=w1&mask
                nlo.append(max(dlo[a]+phi[w0^flip][0]-phi[w0][1], dlo[b]+phi[w1^flip][0]-phi[w1][1]))
                nhi.append(max(dhi[a]+phi[w0^flip][1]-phi[w0][0], dhi[b]+phi[w1^flip][1]-phi[w1][0]))
            dlo,dhi=nlo,nhi
        return max(F(0),max(dlo)),max(dhi)
    phi=np.asarray(phi); words=np.arange(1<<m); tails=words&mask
    dp=np.zeros(V)
    for t in range(m):
        costs=phi[words^(1<<(m-1-t))]-phi
        dp=(dp[tails]+costs).reshape(V,2).max(axis=1)
    return float(dp.max())


def certify(data, dps):
    t0=time.perf_counter(); iv.dps=dps
    assert data['model']=='periodic_soft_east_activity'
    assert data['bit_order']=='leftmost_is_least_significant'
    assert data['density_convention']=='U=sum(density); no implicit s factor'
    R=data['R']; m=2*R-1; center=R-1
    assert R>=2 and data['n_min']==m and len(data['density'])==1<<R
    phi,b,p=local_tables(data,True)
    dl,du=maxplus(phi,m,True)
    cl=[]; cu=[]
    for j in range(m):
        if j==center:continue
        low=F(0); high=F(0)
        for word in range(1<<m):
            if (word>>center)&1 or (word>>j)&1:continue
            w2=word^(1<<j)
            low=max(low,p[word][0]-p[w2][1],p[w2][0]-p[word][1])
            high=max(high,p[word][1]-p[w2][0],p[w2][1]-p[word][0])
        cl.append(low);cu.append(high)
    al,au=sum(cl,F(0)),sum(cu,F(0))
    words=[w for w in range(1<<m) if not((w>>center)&1)]
    bl=min(b[w][0]+b[w^(1<<center)][0] for w in words)
    bu=min(b[w][1]+b[w^(1<<center)][1] for w in words)
    correction = du*du/(4*bl*(1-au)**2) if au<1 and bl>0 else None
    return {'status':'INTERVAL_CONSTANTS_VERIFIED' if correction is not None else 'CONDITIONS_INCONCLUSIVE',
            'dps':dps,'n_min':m,'exact_long_window_threshold':2*m-1,
            'alpha':rec(al,au),'beta':rec(bl,bu),'d':rec(dl,du),
            'correction_upper':dec(correction,True) if correction is not None else None,
            'B_upper':dec(max(x[1] for x in b),True),
            'influence_upper_by_offset':{str(j-center):dec(cu[j if j<center else j-1],True) for j in range(m) if j!=center},
            'states':1<<(m-1),'layers':m,'seconds':time.perf_counter()-t0,
            'scope':'Uniform correction to Q_n for all n>=n_min, conditional on the proved spectral lemma. This file does not enclose Q_n or Q_infinity.'}


def transfer(data, n=None):
    R=data['R']; m=2*R-1; D=1<<(R-1)
    c=float(F(data['c'])); v=np.asarray(data['density'],float)
    e=np.arange(1<<R);tail=e&(D-1);head=e>>1
    te=np.where(e>>(R-1),c,1-c)*np.exp(2*v)
    T=np.zeros((D,D));T[tail,head]=te
    words=np.arange(1<<m)
    path=np.ones(len(words))
    for j in range(R):path*=te[(words>>j)&((1<<R)-1)]
    start=words&(D-1);end=words>>R
    phi=np.asarray(local_tables(data)[0])
    if n is None:
        vals,left,right=eig(T,left=True,right=True)
        j=int(np.argmax(vals.real));rho=float(vals[j].real)
        l=left[:,j].real;r=right[:,j].real
        if l.sum()<0:l=-l
        if r.sum()<0:r=-r
        assert l.min()>0 and r.min()>0
        prob=l[start]*path*r[end]/(rho**R*(l@r))
    else:
        assert n>=m
        A=np.linalg.matrix_power(T,n-R)
        den=np.trace(np.linalg.matrix_power(T,n))
        prob=path*A[end,start]/den
    assert abs(prob.sum()-1)<2e-10
    return float(prob@phi)


def direct_ring(data,n,eigen=False):
    R=data['R'];c,k,s=[float(F(data[x])) for x in ('c','kappa','s')]
    v=np.asarray(data['density'],float);states=np.arange(1<<n)
    bits=(states[:,None]>>np.arange(n))&1
    U=np.zeros(len(states))
    for i in range(n):
        w=sum(bits[:,(i+j)%n]<<j for j in range(R))
        U+=v[w]
    rates=(k+np.roll(bits,1,axis=1))*(c+(1-2*c)*bits)
    delta=np.column_stack([U[states^(1<<i)]-U for i in range(n)])
    b=rates*np.exp(-s+delta)
    V=(rates*np.expm1(-s+delta)).sum(axis=1)
    pi=np.prod(np.where(bits,c,1-c),axis=1)
    nu=pi*np.exp(2*(U-U.max()));nu/=nu.sum()
    Q=float(nu@V/n)
    d=max(float(np.max(np.abs(V[states^(1<<i)]-V))) for i in range(n))
    ans={'Q':Q,'d':d,'n':n}
    if eigen:
        row=np.repeat(states,n);col=(states[:,None]^(1<<np.arange(n))).ravel()
        off=(np.exp(-s)*np.sqrt(c*(1-c))*(k+np.roll(bits,1,axis=1))).ravel()
        A=coo_matrix((off,(row,col)),shape=(len(states),len(states))).tocsr()-diags(rates.sum(axis=1))
        val,vec=eigsh(A,k=1,which='LA',v0=np.sqrt(pi),tol=1e-13)
        h=vec[:,0]/np.sqrt(pi)
        if h.sum()<0:h=-h
        assert h.min()>0
        logh=np.log(h)
        mismatch=delta-np.column_stack([logh[states^(1<<i)]-logh for i in range(n)])
        # 2*sinh(z/2)^2 avoids cancellation at small z.
        cosh_identity=float(nu@(b*(2*np.sinh(mismatch/2)**2)).sum(axis=1)/n)
        actual=float(val[0]/n-Q)
        true_rates=rates*np.exp(-s+delta-mismatch)
        kl=float(nu@(b*mismatch-b+true_rates).sum(axis=1)/n)
        ans.update(lambda_density=float(val[0]/n),actual_gap=actual,cosh_identity=cosh_identity,
                   identity_error=abs(actual-cosh_identity),kl_rate=kl,kl_error=abs(kl-cosh_identity),
                   eigen_residual=float(np.max(np.abs(A@vec[:,0]-val[0]*vec[:,0]))))
    return ans


def verification_tests():
    tests=[]
    for R in [2,3]:
        m=2*R-1
        # Deliberately non-tiny, nontrivial candidate; both tilt signs.
        for s in ['-0.03','0.04']:
            d={'R':R,'c':'0.3','kappa':'0.1','s':s,
               'density':[str(((i*7+3)%11-5)/80) for i in range(1<<R)]}
            phi=np.asarray(local_tables(d)[0]); dp=maxplus(phi,m)
            exhaustive=max(sum(phi[((w>>j)&((1<<m)-1))^(1<<(m-1-j))]-phi[(w>>j)&((1<<m)-1)] for j in range(m)) for w in range(1<<(2*m-1)))
            assert abs(dp-exhaustive)<1e-12
            for n in range(m,2*m):
                direct=direct_ring(d,n,eigen=True);q=transfer(d,n)
                assert abs(q-direct['Q'])<2e-12
                assert direct['d']<=dp+2e-12
                if n>=2*m-1:assert abs(direct['d']-dp)<2e-12
                assert direct['identity_error']<3e-12 and direct['kl_error']<3e-12
                assert direct['actual_gap']>=-3e-12
                tests.append({'R':R,'s':s,'dp':dp,'exhaustive_window':exhaustive,'transfer_error':abs(q-direct['Q']),**direct})
    # A nonconstant one-site product candidate checks the sharp quadratic order.
    asymptotic=[]
    c=0.3;k=0.1
    for eta in [0.02,0.01,0.005]:
        z=1-c+c*np.exp(2*eta);p=c*np.exp(2*eta)/z
        gap=(k+p)*c*(1-c)*(np.expm1(eta)**2)/z
        a=(1-c+c*np.exp(eta))/np.sqrt(z)
        asymptotic.append({'eta':eta,'gap':gap,'gap_over_eta_squared':gap/eta**2,'overlap_base':a})
        assert 0<a<1
    return {'status':'PASS_FLOATING_INDEPENDENT_CROSSCHECKS','ring_checks':tests,'product_sharpness':asymptotic,
            'note':'These crosschecks are diagnostics, not interval spectral certificates.'}


def bridge():
    data=json.loads((ROOT.parent/'baseline/results/practical/perron_k1_s0.02_N18_R7_locality.json').read_text())
    eps=F(data['flip_log_error']['upper']);q=F(7,10);cap=F(data['local_cap']);a=F(data['local_a'])
    iv.dps=100
    B=I(F(7,5))*iv.exp(I(F(1,50)+cap/q))
    error=B*(iv.exp(I(eps))+iv.exp(-I(eps))-2)/2
    alpha=I(cap)/(2*I(q*q*a)*iv.exp(I(1)))
    # tanh(eps/2)=(exp(eps)-1)/(exp(eps)+1).
    obs=(iv.exp(I(eps))-1)/(iv.exp(I(eps))+1)/(1-alpha)
    return {'status':'PROP_L_CONDITIONAL_CORRECTION','epsilon_D_upper':str(eps),
            'B':rec(*bounds(B)),'quadratic_correction':rec(*bounds(error)),
            'norm_based_alpha':rec(*bounds(alpha)),'single_site_observable_error':rec(*bounds(obs)),
            'candidate_sha256':data['candidate_sha256'],
            'scope':'Uses the previously certified Proposition L flip error. Does not yet enclose the variational center.'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dps',type=int,default=80)
    ap.add_argument('--case',default=None);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();out=args.output.resolve()
    if not out.is_relative_to(ROOT):raise ValueError('Output must stay inside extension')
    if out.exists() and any(out.iterdir()):raise FileExistsError('Choose a new run directory; archive old results explicitly')
    out.mkdir(parents=True,exist_ok=True)
    test=verification_tests();(out/'independent_tests.json').write_text(json.dumps(test,indent=2)+'\n')
    (out/'proposition_L_bridge.json').write_text(json.dumps(bridge(),indent=2)+'\n')
    print('Independent crosschecks passed:',len(test['ring_checks']),flush=True)
    rows=[]
    for path in candidate_paths():
        data=json.loads(path.read_text())
        if 'density' not in data:continue
        if args.case and args.case not in path.stem:continue
        t0=time.perf_counter();cert=certify(data,args.dps)
        cert.update(candidate=path.name,candidate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    c=data['c'],kappa=data['kappa'],s=data['s'],R=data['R'],reference_N=data['reference_N'])
        (out/(path.stem+'_constants.json')).write_text(json.dumps(cert,indent=2)+'\n')
        qinf=transfer(data);qn=transfer(data,data['reference_N'])
        old=json.loads(path.with_name(path.stem+'_verified.json').read_text())
        width=F(old['upper_density'])-F(old['lower_density'])
        correction=F(cert['correction_upper']) if cert['correction_upper'] else None
        row={'candidate':path.name,'kappa':data['kappa'],'s':data['s'],'N':data['reference_N'],'R':data['R'],
             'alpha_upper':float(F(cert['alpha']['upper'])),'beta_lower':float(F(cert['beta']['lower'])),
             'd_upper':float(F(cert['d']['upper'])),'correction_upper':float(correction) if correction else None,
             'old_graph_width':float(width),'old_width_over_correction':float(width/correction) if correction else None,
             'Q_infinity_FLOAT_ONLY':qinf,'Q_N_FLOAT_ONLY':qn,'Q_N_minus_Q_infinity_FLOAT_ONLY':qn-qinf,
             'seconds':time.perf_counter()-t0}
        rows.append(row);print(json.dumps(row),flush=True)
    (out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    if rows:
        with (out/'summary.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'runtime.json').write_text(json.dumps({'python':platform.python_version(),'numpy':np.__version__,
        'mpmath':mpmath.__version__,'interval_dps':args.dps,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')


if __name__=='__main__':main()
