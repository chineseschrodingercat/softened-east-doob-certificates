"""Rational reconstruction and PSD correction, independent of solver feasibility.

Every saved dual matrix is a rational matrix. For symmetric A, a rational
factor L certifies lambda_min(A) >= mu - ||A-mu I-LL^T||_infinity.
Integer matrix multiplication evaluates that error exactly. No floating
eigenvalue is used as a certified endpoint.
"""
from fractions import Fraction as F
from math import lcm
import numpy as np
from mpmath import iv

def endpoint(t):
    sign,m,e,b=t
    if b<0: raise ArithmeticError('nonfinite')
    return F((-1 if sign else 1)*m)*(F(2)**e)

def scalar_interval(x): return endpoint(x._mpi_[0]),endpoint(x._mpi_[1])
def I(x):
    x=F(x); return iv.mpf(x.numerator)/x.denominator
def decimal(x,up=False,p=35):
    a=F(x)*10**p; n=-((-a.numerator)//a.denominator) if up else a.numerator//a.denominator
    return ('-' if n<0 else '')+str(abs(n)//10**p)+'.'+str(abs(n)%10**p).zfill(p)
def rat(a,digits=14):
    a=np.asarray(a,dtype=float); a=(a+a.T)/2
    return np.array([[F(int(round(float(v)*10**digits)),10**digits) for v in row] for row in a],dtype=object)
def dumpmat(a): return [[str(v) for v in row] for row in a]
def loadmat(a): return np.array([[F(v) for v in row] for row in a],dtype=object)
def eye(d): return np.eye(d,dtype=int).astype(object)

def local_h(c,k,s):
    c,k,s=map(F,(c,k,s)); iv.dps=100
    alo,ahi=scalar_interval(iv.exp(-I(s))*iv.sqrt(I(c*(1-c))))
    a=F(((alo+ahi)*10**18/2).__floor__(),10**18)
    error=(k+1)*max(a-alo,ahi-a)
    h=np.kron(np.diag([k,k+1]),np.array([[c,-a],[-a,1-c]],dtype=object))
    return h,error,{'a':str(a),'a_lower':str(alo),'a_upper':str(ahi),'local_operator_error':str(error)}

def embed(a,sites,n):
    sites=list(sites); rest=[j for j in range(n) if j not in sites]
    d=1<<n; out=np.zeros((d,d),dtype=object)
    for i in range(d):
        bi=[(i>>(n-1-j))&1 for j in range(n)]
        ai=sum(bi[j]<<(len(sites)-1-t) for t,j in enumerate(sites))
        for aj in range(1<<len(sites)):
            bj=bi.copy()
            for t,j in enumerate(sites): bj[j]=(aj>>(len(sites)-1-t))&1
            q=sum(bj[j]<<(n-1-j) for j in range(n))
            out[i,q]=a[ai,aj]
    return out

def open_h(h,n): return sum((embed(h,[j,j+1],n) for j in range(n-1)),np.zeros((1<<n,1<<n),dtype=object))

def trace_indices(n,sites):
    sites=list(sites); rest=[j for j in range(n) if j not in sites]; rows=[]
    for r in range(1<<len(rest)):
        ids=[]
        for a in range(1<<len(sites)):
            bits=[0]*n
            for t,j in enumerate(sites):bits[j]=(a>>(len(sites)-1-t))&1
            for t,j in enumerate(rest):bits[j]=(r>>(len(rest)-1-t))&1
            ids.append(sum(bits[j]<<(n-1-j) for j in range(n)))
        rows.append(ids)
    return rows

def ptr(a,n,sites):
    return sum((a[np.ix_(ids,ids)] for ids in trace_indices(n,sites)))

def make_psd_certificate(a):
    af=np.asarray(a,float); af=(af+af.T)/2
    mineig=float(np.linalg.eigvalsh(af)[0]); attempts=[]
    for gap in [1e-12,1e-10,1e-8,1e-6]:
        mu=F(str(mineig-gap))
        try: l=np.linalg.cholesky(af-float(mu)*np.eye(len(a))); break
        except np.linalg.LinAlgError: attempts.append(gap)
    else:raise ArithmeticError('factor generation failed')
    scale=10**16
    li=np.array([[int(round(float(v)*scale)) for v in row] for row in l],dtype=object)
    cert={'mu':str(mu),'factor_scale':str(scale),'factor':[[str(v) for v in row] for row in li],
          'failed_factor_gaps':attempts,'eigenvalue_diagnostic':mineig}
    lower=check_psd_certificate(a,cert);cert['lower_exact']=str(lower)
    return cert

def check_psd_certificate(a,cert):
    a=np.asarray(a,dtype=object); d=len(a)
    if a.shape!=(d,d) or not np.array_equal(a,a.T):raise ValueError('not exactly symmetric')
    mu=F(cert['mu']); scale=int(cert['factor_scale']); li=np.array([[int(v) for v in row] for row in cert['factor']],dtype=object)
    if li.shape!=a.shape or scale<=0:raise ValueError('bad factor')
    b=a-mu*eye(d);den=1
    for x in b.flat:den=lcm(den,F(x).denominator)
    ai=np.array([[int(F(v)*den) for v in row] for row in b],dtype=object)
    residual=ai*scale**2-den*(li@li.T)
    radius=F(max(sum(abs(int(v)) for v in row) for row in residual),den*scale**2)
    lower=mu-radius
    if 'lower_exact' in cert and F(cert['lower_exact'])!=lower:raise ValueError('PSD residual mismatch')
    return lower

def mps_maps(A,k):
    D=len(A[0]); DD=D*D
    V=np.zeros((DD,1<<k),dtype=object)
    for w in range(1<<k):
        prod=eye(D)
        for j in range(k):prod=prod@A[(w>>(k-1-j))&1]
        V[:,w]=prod.reshape(-1)
    # L consumes (physical,left-bond,right-bond); R consumes (left,right,physical).
    L=np.zeros((DD,2*DD),dtype=object);R=np.zeros((DD,DD*2),dtype=object)
    for a in range(D):
      for b in range(D):
       for t in range(D):
        for x in range(2):
         L[a*D+b,x*DD+t*D+b]=A[x][a,t]
         R[a*D+b,(a*D+t)*2+x]=A[x][t,b]
    return V,L,R

def inverse(a):
    n=len(a);aug=np.concatenate([a.copy(),eye(n)],axis=1)
    for j in range(n):
        pivot=next(i for i in range(j,n) if aug[i,j]!=0)
        if pivot!=j:aug[[j,pivot]]=aug[[pivot,j]]
        aug[j,:]=[F(x)/F(aug[j,j]) for x in aug[j,:]]
        for i in range(n):
            if i!=j:aug[i,:]-=aug[i,j]*aug[j,:]
    return aug[:,n:]

def gauged_maps(A,k,gauge=None):
    V,L,R=mps_maps(A,k)
    if gauge is None:return V,L,R
    W=loadmat(gauge);Winv=inverse(W)
    return W@V,W@L@np.kron(eye(2),Winv),W@R@np.kron(Winv,eye(2))

def slacks(spec):
    h,error,meta=local_h(spec['c'],spec['kappa'],spec['s'])
    kind=spec['method']; dual=[loadmat(v) for v in spec['duals']]
    if kind=='lti':
        n=spec['m'];Y=dual[0]
        return [embed(h,[0,1],n)+np.kron(eye(2),Y)-np.kron(Y,eye(2))],[F(1)],F(1),error
    if kind=='eisert':
        m,t=spec['m'],spec['overlap'];yl,yr=dual
        sw=2*open_h(h,m)-embed(yl,range(m-t,m),m)-embed(yr,range(t),m)
        ss=embed(h,[t-1,t],2*t)+np.kron(yl,eye(1<<t))+np.kron(eye(1<<t),yr)
        return [sw,ss],[F(1),F(1)],F(2*m-1),error
    if kind=='kull_rg':
        k0,n=spec['k0'],spec['last_k'];A=[loadmat(a) for a in spec['mps']];D=len(A[0]);DD=D*D
        V,L,R=gauged_maps(A,k0,spec.get('gauge'));VxI=np.kron(V,eye(2));IxV=np.kron(eye(2),V)
        LxI=np.kron(L,eye(2));IxR=np.kron(eye(2),R)
        alpha=dual[0];pairs=[dual[j:j+2] for j in range(1,len(dual),2)]
        out=[embed(h,[0,1],k0+1)+np.kron(eye(2),alpha)-np.kron(alpha,eye(2))+VxI.T@pairs[0][0]@VxI+IxV.T@pairs[0][1]@IxV]
        traces=[F(1)]
        for k in range(k0,n+1):
            bl,br=pairs[k-k0];a=-np.kron(eye(2),bl)-np.kron(br,eye(2))
            if k<n:
                gl,gr=pairs[k-k0+1];a+=LxI.T@gl@LxI+IxR.T@gr@IxR
            out.append(a)
            Vk,_,_=gauged_maps(A,k,spec.get('gauge'));traces.append(sum(F(x)**2 for x in Vk.flat))
        return out,traces,F(1),error
    if kind=='ffs_t0':
        n=spec['m'];Y=dual[0];zfactor=loadmat(spec['zfactor']);Z=zfactor@zfactor.T
        a=embed(h,[0,1],n)+np.kron(eye(2),Y)-np.kron(Y,eye(2))
        ops=ffs_operators(h,n)
        for i in range(4):
            for j in range(4):a-=Z[i,j]*ops[i][j]
        for y,b in zip(spec['stationary_duals'],ffs_stationary(h,n)):a+=F(y)*b
        a=(a+a.T)/2
        ys=list(map(F,spec['stationary_duals']))
        correction=error*(4*sum(abs(v) for v in Z.flat)+4*sum(abs(v) for v in ys[:4])+8*sum(abs(v) for v in ys[4:]))
        if F(spec['coefficient_correction'])!=correction:raise ValueError('invalid FFS coefficient correction')
        return [a],[F(1)],F(1),error+correction
    raise ValueError(kind)

def ffs_operators(h,n):
    H=open_h(h,n); basis=[]
    for i in range(2):
        for j in range(2):
            a=np.zeros((2,2),dtype=object);a[i,j]=1;basis.append(embed(a,[1],n))
    return [[(a.T@(H@b-b@H)+(a.T@(H@b-b@H)).T)/2 for b in basis] for a in basis]
def ffs_stationary(h,n):
    H=open_h(h,n);out=[]
    for i in range(2):
        for j in range(2):
            b=np.zeros((2,2),dtype=object);b[i,j]=1;a=embed(b,[1],n);B=H@a-a@H;out.append((B+B.T)/2)
    ops=ffs_operators(h,n)
    for i in range(4):
        for j in range(i):out.append(ops[i][j]-ops[j][i])
    return out

def verify_sdp(spec):
    kind=spec['method']
    if kind in ['lti','eisert'] and spec['m']<2:raise ValueError('at least two physical sites required')
    if kind=='eisert' and not 1<=spec['overlap']<=spec['m']:raise ValueError('invalid overlap support')
    if kind=='ffs_t0' and spec['m']<3:raise ValueError('FFS requires both neighbors of the perturbed site')
    if kind=='kull_rg' and not 1<=spec['k0']<=spec['last_k']:raise ValueError('invalid compressed hierarchy range')
    required={'lti':lambda:spec['m'],'eisert':lambda:max(spec['m'],2*spec['overlap']),
              'kull_rg':lambda:spec['last_k']+2,'ffs_t0':lambda:max(3,spec['m'])}[spec['method']]()
    if not isinstance(spec['n_min'],int) or spec['n_min']<required:
        raise ValueError('claimed n_min is below the support required by the SDP')
    mats,tracebounds,div,error=slacks(spec)
    if len(mats)!=len(spec['psd_certificates']):raise ValueError('wrong certificate count')
    lows=[check_psd_certificate(a,c) for a,c in zip(mats,spec['psd_certificates'])]
    value=lows[0]
    for j in range(1,len(lows)):
        # Eisert states are trace one; Kull compressed states have bounded trace.
        value+=lows[j] if spec['method']=='eisert' else min(lows[j],F(0))*tracebounds[j]
    elower=value/div-error;upper=-elower
    return {'status':'EXACT_DUAL_CERTIFICATE','upper_density':decimal(upper,True),
            'energy_lower_exact':str(elower),'matrix_lower_exact':[str(x) for x in lows],
            'trace_bounds':[str(x) for x in tracebounds],'coefficient_error_exact':str(error),
            'n_min':spec['n_min'],'scope':'All periodic rings n>=n_min and their thermodynamic limit'}
