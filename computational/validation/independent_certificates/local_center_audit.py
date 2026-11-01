"""Independent larger-state interval transfer and exhaustive residual audit."""
from independent_audit import *

def encode(bits): return sum(x<<j for j,x in enumerate(bits))

def local_functions(tr):
    R=tr['R']; m=2*R-1;center=R-1;c,k,s=map(F,[tr['c'],tr['kappa'],tr['s']]);v=list(map(F,tr['density']))
    phi=[];rates=[];cond=[]
    for word in range(2**m):
        bits=[(word>>j)&1 for j in range(m)]; changed=bits.copy();changed[center]^=1
        delta=sum(v[encode(changed[j:j+R])]-v[encode(bits[j:j+R])] for j in range(R))
        rate=(k+bits[center-1])*(c if bits[center]==0 else 1-c)
        b=fraciv(rate)*iv.exp(fraciv(delta-s))
        phi.append(bounds(b-fraciv(rate)));rates.append(bounds(b))
        odds=fraciv(c/(1-c))*iv.exp(fraciv(2*delta*(1-2*bits[center])))
        cond.append(bounds(odds/(1+odds)))
    return phi,rates,cond

def finite_center(tr,phi,n=18):
    # This transfer uses 2^(2R-2) states; production uses 2^(R-1).
    # It carries the entire residual window in one transition, so one
    # trace insertion suffices and no R-edge bridge formula is used.
    R=tr['R'];m=2*R-1;D=2**(m-1);c=F(tr['c']);v=list(map(F,tr['density']))
    te=[fraciv(c if w>>(m-1) else 1-c)*iv.exp(fraciv(2*v[w>>(m-R)])) for w in range(2**m)]
    P=[[iv.mpf(int(i==j)) for j in range(D)] for i in range(D)]
    for step in range(n-1):
        Q=[[iv.mpf(0) for j in range(D)] for i in range(D)]
        for i in range(D):
            for tail in range(D):
                z=P[i][tail]
                if z==0:continue
                for x in range(2):
                    w=tail+x*D;Q[i][w>>1]+=z*te[w]
        P=Q
    num=iv.mpf(0);den=iv.mpf(0)
    for w in range(2**m):
        tail=w&(D-1);head=w>>1;weight=te[w]*P[head][tail]
        f=iv.mpf([fraciv(phi[w][0]).a,fraciv(phi[w][1]).b])
        num+=weight*f;den+=weight
    assert bounds(den)[0]>0
    return bounds(num/den),D

def local_audit(path):
    tr=load(path/'trial_graph.json');r=load(path/'result.json');R=tr['R'];m=2*R-1;D=2**(m-1)
    phi,b,p=local_functions(tr);zl=list(map(F,tr['z_lower']));zu=list(map(F,tr['z_upper']))
    graph_lo=min(a+zl[w&(D-1)]-zl[w>>1] for w,(a,bv) in enumerate(phi))
    graph_hi=max(bv+zu[w&(D-1)]-zu[w>>1] for w,(a,bv) in enumerate(phi))
    assert F(r['graph']['lower_density'])<=graph_lo and graph_hi<=F(r['graph']['upper_density'])
    alo=F(0);ahi=F(0)
    for j in range(m):
        if j==R-1:continue
        lows=[];highs=[]
        for w in range(2**m):
            if (w>>(R-1))&1 or (w>>j)&1:continue
            a,aa=p[w];z,zz=p[w^(1<<j)]
            lows.extend([a-zz,z-aa]);highs.extend([aa-z,zz-a])
        alo+=max(F(0),max(lows));ahi+=max(F(0),max(highs))
    blo=min(b[w][0]+b[w^(1<<(R-1))][0] for w in range(2**m) if not((w>>(R-1))&1))
    bhi=min(b[w][1]+b[w^(1<<(R-1))][1] for w in range(2**m) if not((w>>(R-1))&1))
    # Exhaustive long-window maximum over every 4R-3 bit configuration,
    # independent of the production max-plus dynamic program.
    dlo=F(0);dhi=F(0)
    for x in range(2**(2*m-1)):
        dl=F(0);dh=F(0)
        for offset in range(m):
            w=(x>>offset)&(2**m-1);ww=w^(1<<(m-1-offset))
            dl+=phi[ww][0]-phi[w][1];dh+=phi[ww][1]-phi[w][0]
        dlo=max(dlo,dl);dhi=max(dhi,dh)
    co=r['constants']
    for name,lo,hi in [('alpha',alo,ahi),('beta',blo,bhi),('d',dlo,dhi)]:
        assert F(co[name]['lower'])<=lo<=hi<=F(co[name]['upper']),name
    possible=(ahi<1 and blo>0);assert possible==(co['correction_upper'] is not None)
    if possible:assert dhi*dhi/(4*blo*(1-ahi)**2)<=F(co['correction_upper'])
    (ql,qh),states=finite_center(tr,phi)
    assert F(r['finite_center']['lower'])<=ql<=qh<=F(r['finite_center']['upper'])
    return {'run':path.name,'range':R,'exact_graph_edges':2**m,'exhaustive_residual_words':2**(2*m-1),
            'quadratic_hypotheses_satisfied':possible,'alpha_upper':roundout(ahi,True),'beta_lower':roundout(blo),'total_residual_flip_upper':roundout(dhi,True),
            'strict_n18_center':{'lower':roundout(ql,p=70),'upper':roundout(qh,True,p=70),'width':float(qh-ql),'transfer_states':states},
            'status':'INDEPENDENT_INTERVAL_ENCLOSURE_WITHIN_ALL_FROZEN_ENDPOINTS'}

def monthus():
    base=PACKAGE/'monthus/results';pub=load(base/'published_results.json'); n=3;D=2**n
    B=np.zeros((D,D),dtype=object)
    for x in range(D):
        B[x,x]=12-2*sum(1+((x>>((i-1)%n))&1) for i in range(n))
        for i in range(n):B[x,x^(1<<i)]=1+((x>>((i-1)%n))&1)
    v=np.ones(D,dtype=object)
    for _ in range(100):v=B@v
    bounds3=[(F(x,y)-12)/12 for x,y in zip(B@v,v)]
    assert min(bounds3)==F(pub['finite_ring']['lower_exact']) and max(bounds3)==F(pub['finite_ring']['upper_exact'])
    a,b,d=F('.11'),F('.187'),F('.516');assert a+2*b+d==1
    assert F('.1434')**2<a*b and F('.3106')**2<b*d
    Fsimple=iv.sqrt(fraciv(a*b))+iv.sqrt(fraciv(b*d))/2-fraciv(a+3*b/2+d/2)
    assert bounds(Fsimple)[0]>F('-.3498')
    mu=[F('.384')]+[F('.132') if x.bit_count()==1 else F('.055') for x in range(1,8)]
    assert sum(mu)==1
    for i in range(3):
        cells={(l,r):sum(mu[x] for x in range(8) if ((x>>i)&1)==l and ((x>>((i+1)%3))&1)==r) for l,r in itertools.product(range(2),repeat=2)}
        assert cells=={(0,0):d,(1,0):b,(0,1):b,(1,1):a}
    output={'finite_n3':{'lower_exact':str(min(bounds3)),'upper_exact':str(max(bounds3))},'feasible_F_lower':roundout(bounds(Fsimple)[0]),'global_pair_realizable':True,'certificates':[]}
    for name in ['uniform_R3_simple_candidate.json','uniform_R4_candidate.json']:
        tr=load(base/name);R=tr['R'];m=2*R-1;fs=list(map(F,tr['factor']));zl=list(map(F,tr['z_lower']));zu=list(map(F,tr['z_upper']));loc=[];low=[];high=[]
        assert min(fs)>0 and tr['n_min']==m
        for w in range(2**m):
            bits=[(w>>j)&1 for j in range(m)];b2=bits.copy();b2[R-1]^=1
            ratio=math.prod(fs[encode(b2[j:j+R])]/fs[encode(bits[j:j+R])] for j in range(R))
            phi=F(1+bits[R-2],2)*(ratio/2-1);loc.append(phi)
            low.append(phi+zl[w&(2**(m-1)-1)]-zl[w>>1]);high.append(phi+zu[w&(2**(m-1)-1)]-zu[w>>1])
        expected=pub[f'R{R}_all_n_at_least_{m}'];assert min(low)==F(expected['lower_exact']) and max(high)==F(expected['upper_exact'])
        # Full exact product trial and original-generator action on the
        # minimum covered ring, including every physical flip.
        h=[]
        for x in range(2**m):h.append(math.prod(fs[sum(((x>>((i+j)%m))&1)<<j for j in range(R))] for i in range(m)))
        for x in range(2**m):
            residual=sum(F(1+((x>>((i-1)%m))&1),2)*(h[x^(1<<i)]/(2*h[x])-1) for i in range(m))
            local=sum(loc[sum(((x>>((i+j)%m))&1)<<j for j in range(m))] for i in range(m))
            assert residual==local and m*min(low)<=residual<=m*max(high)
        output['certificates'].append({'R':R,'n_min':m,'edges':2**m,'lower_exact':str(min(low)),'upper_exact':str(max(high)),'minimal_ring_all_flips_exact':True})
    assert F(output['certificates'][0]['upper_exact'])<F('-.3499')<F('-.3498')
    return output

def legacy_scope():
    rows=load(PACKAGE/'certificates/published/certificates.json');cases=load(PACKAGE/'certificates/cases.json');assert len(rows)==22 and sum(x['quadratic'] for x in rows)==20
    cert=[]
    for row in rows:
        case=next(x for x in cases if x['candidate']==row['candidate']);p=(PACKAGE/'certificates'/case['input_path']).resolve()
        assert hashlib.sha256(p.read_bytes()).hexdigest()==case['sha256']==row['sha256']
        tr=load(p);g=row['graph'];assert g['all_integer_n_at_least']==2*tr['R']-1
        if row['quadratic']:
            co=row['constants'];a,b,d=F(co['alpha']['upper']),F(co['beta']['lower']),F(co['d']['upper']);assert a<1 and b>0
            correction=d*d/(4*b*(1-a)**2)
            al,bu,dl=F(co['alpha']['lower']),F(co['beta']['upper']),F(co['d']['lower'])
            correction_low=dl*dl/(4*bu*(1-al)**2)
            assert correction_low<=F(co['correction_upper'])<=correction+F(1,10**65)
            inf=row['infinite'];ql=F(inf['Q_infinity']['lower']);qu=F(inf['Q_infinity']['upper']);upper=F(inf['spectral_density_limit']['upper'])
            assert F(inf['spectral_density_limit']['lower'])==ql
            assert upper>=qu+F(co['correction_upper'])
            assert F(inf['spectral_width_upper'])==upper-ql
        cert.append({'candidate':row['candidate'],'n_min':2*tr['R']-1,'quadratic':row['quadratic'],'graph_volume_scope':'all n >= n_min','quadratic_scope':'thermodynamic limit with Q_infinity, not a finite-n center'})
    return cert

def main():
    configure_output()
    out=[];start=time.perf_counter()
    for p in sorted((EXT/'results').glob('*/trial_graph.json')):
        result=local_audit(p.parent);out.append(result);save('independent_local_and_finite_centers.json',out);print(p.parent.name,result['status'],flush=True)
    save('independent_monthus.json',monthus());print('monthus exact PASS',flush=True)
    save('legacy_scope_audit.json',legacy_scope());print('legacy scope PASS',flush=True)
    save('local_center_summary.json',{'new_trials':len(out),'finite_center_interval_digits':130,'elapsed_seconds':time.perf_counter()-start,'status':'PASS'})

if __name__=='__main__':main()
