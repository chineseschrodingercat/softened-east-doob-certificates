"""Finite-state diagnostics of the new coupling lemma; not a locality proof."""
from pathlib import Path
import argparse, hashlib, json, math, platform
import numpy as np
import scipy
from scipy.linalg import eigh, expm
import mpmath as mp

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--output',type=Path,default=Path('new_results_replay'))
OUT=(ROOT/parser.parse_args().output).resolve()
if not OUT.is_relative_to(ROOT.resolve()) or OUT==ROOT.resolve():
    raise ValueError('Use a new output subdirectory of the script directory')
OUT.mkdir(parents=True,exist_ok=False)
records=[]

def state_data(n):
    states=np.arange(2**n)
    bits=(states[:,None]>>np.arange(n))&1
    flip=states[:,None]^(1<<np.arange(n))
    return states,bits,flip

def generator(rates,flip):
    G=np.zeros((len(rates),len(rates)))
    G[np.arange(len(rates))[:,None],flip]=rates
    G[np.diag_indices(len(rates))]=-rates.sum(axis=1)
    return G

def diagnose(q,b,flip,bits,label):
    n=q.shape[1]
    influence=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i!=j: influence[i,j]=np.max(np.abs(q[flip[:,j],i]-q[:,i]))
    gamma=float(influence.sum(axis=1).max())
    A=float(q.max())
    eta=float(np.max(np.abs(np.log(b/q))))
    Gq,Gb=generator(q,flip),generator(b,flip)
    fs=[bits[:,0].astype(float),bits[:,0]*bits[:,1]+.3*bits[:,1],np.ones(len(bits))]
    for t in [.025,.25,.8]:
        Pq,Pb=expm(t*Gq),expm(t*Gb)
        rows=[]
        for f in fs:
            osc=float(np.abs(f[flip]-f[:,None]).max(axis=0).sum())
            gt=Pq@f
            evolved=float(np.abs(gt[flip]-gt[:,None]).max(axis=0).sum())
            err=float(np.max(np.abs((Pb-Pq)@f)))
            J=math.expm1(t*gamma)/gamma if gamma else t
            bound=A*math.expm1(eta)*J*osc
            assert evolved<=math.exp(t*gamma)*osc+2e-12,(label,t,'oscillation')
            assert err<=bound+2e-12,(label,t,'Duhamel')
            rows.append({'error':err,'bound':bound,'oscillation':osc,'evolved_oscillation':evolved})
        records.append({'case':label,'n':n,'time':t,'gamma_row':gamma,
                        'gamma_column':float(influence.sum(axis=0).max()),'eta':eta,'functions':rows})

for n in [3,4,5,6]:
    states,bits,flip=state_data(n)
    for c,kappa,s in [(.3,1.,.02),(.3,.1,.01),(.3,.01,.001),(.5,.1,-.03)]:
        pi=(c**bits.sum(axis=1))*((1-c)**(n-bits.sum(axis=1)))
        rates=(kappa+np.roll(bits,1,axis=1))*(c*(1-bits)+(1-c)*bits)
        L=generator(rates*np.exp(-s),flip)
        L[np.diag_indices(len(states))]=-rates.sum(axis=1)
        S=np.sqrt(pi)[:,None]*L/np.sqrt(pi)[None,:]
        assert np.max(np.abs(S-S.T))<1e-13
        vals,vecs=eigh(S)
        h=vecs[:,-1]/np.sqrt(pi)
        if h[0]<0:h=-h
        assert h.min()>0
        q=np.exp(-s)*rates*h[flip]/h[:,None]
        U=.025*bits.sum(axis=1)-.013*(bits*np.roll(bits,1,axis=1)).sum(axis=1)
        b=np.exp(-s)*rates*np.exp(U[flip]-U[:,None])
        diagnose(q,b,flip,bits,f'East c={c} kappa={kappa} s={s}')
    # Asymmetric influences guard the row/column orientation in the new proof.
    strengths=np.arange(1,n+1)/n
    q=(.2+strengths[None,:])*(1+.2*bits+.9*bits[:,[0]])
    F=(bits*np.linspace(-.04,.07,n)[None,:]).sum(axis=1)+.02*bits[:,0]*bits[:,-1]
    b=q*np.exp(F[flip]-F[:,None])
    diagnose(q,b,flip,bits,'nontranslation-invariant positive rates')

mp.mp.dps=90
scale=[]
eps=mp.mpf('1e-32')
for ctext in ['0.1','0.3','0.5']:
    c=mp.mpf(ctext);m=min(c,1-c);w=c*(1-c);C0=2/m**2
    K=8
    def k_ok(K):
        return ((mp.log(64*(K+1))+2)/(8*w)<=K/2 and
                32*K*mp.log(mp.e/16)<=mp.mpf('2.5')*mp.log(eps/(2*mp.e)) and
                mp.log(128*(K+2))+1-mp.mpf(K)/2<=126*mp.log(eps/(mp.e*C0)))
    while not k_ok(K):K*=2
    for ktext in ['1','0.1','0.01','0.0001','0.000001']:
        k=mp.mpf(ktext);delta=k/(1+k);Ld=mp.log(mp.e/delta)
        T=int(mp.ceil(K*Ld/delta));ell=32*T
        tau=mp.log(64*T-1)/(8*delta*w);t=T-tau
        logp=ell*mp.log(mp.e/16)
        alpha=1-delta
        logzeta=mp.log(128*(1+t)/delta**2)+max(ell/2*mp.log(alpha),-2*(1-mp.sqrt(alpha))*t)
        assert t>=1 and tau<=mp.mpf(T)/2
        assert logp<=mp.mpf('2.5')*mp.log(eps/(2*mp.e))
        assert logzeta<=126*mp.log(eps/(mp.e*C0))
        bstar=mp.exp(-1)*(eps/2)**97
        beta=min(bstar/(128*T*T),delta/(1024*T*(1+256*T*mp.log(1/m))))
        lower=min(bstar/128,1/(1024*(1+256*mp.log(1/m))))*delta/(T*T)
        upper=delta/(262144*mp.log(1/m)*T*T)
        assert lower<=beta<=upper
        scale.append({'c':ctext,'kappa':ktext,'K':K,'T':T,'ell':ell,
                      'tau_over_T':str(tau/T),'log10_beta':str(mp.log10(beta)),
                      'log10_zeta':str(logzeta/mp.log(10))})

result={'status':'PASS','scope':'Finite-state diagnostics of coupling/Duhamel and numerical checks of symbolic sufficient inequalities; not a proof of Theorem A or an independent review.',
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__},
        'dynamical_cases':len(records),'observable_inequalities':sum(len(x['functions']) for x in records),
        'softening_cases':len(scale),'dynamics':records,'softening':scale}
(OUT/'new_results_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k not in ['dynamics','softening']},indent=2))
