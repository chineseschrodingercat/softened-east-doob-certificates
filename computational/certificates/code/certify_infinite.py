"""Rigorous thermodynamic candidate center from exact decimal Perron seeds.

Power iteration generates seeds only. Validation uses positive-power Hilbert
contraction and outward interval operations. No unverified eigenvector accuracy
or finite-n substitution enters the certified center.
"""
import argparse,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path
from mpmath import mp,iv
from candidate_inputs import candidate_path
from phase0 import ROOT,I,bounds,rec,dec,local_tables


def edges(data,ctx):
    c=F(data['c']);R=data['R'];v=list(map(F,data['density']))
    def conv(x):
        x=F(x);return ctx.mpf(x.numerator)/x.denominator
    return [conv(c if w>>(R-1) else 1-c)*ctx.exp(conv(2*v[w])) for w in range(1<<R)]


def apply(vec,te,R,transpose=False):
    D=1<<(R-1);z=[vec[0]*0 for _ in range(D)]
    for w,t in enumerate(te):
        a=w&(D-1);b=w>>1
        if transpose:z[b]+=t*vec[a]
        else:z[a]+=t*vec[b]
    return z


def seed(data,digits=105,steps=700):
    mp.dps=digits+30;R=data['R'];te=edges(data,mp);D=1<<(R-1)
    vectors=[]
    for transpose in [False,True]:
        u=[mp.mpf(1)]*D
        for _ in range(steps):
            u=apply(u,te,R,transpose);scale=u[0];u=[x/scale for x in u]
        vectors.append([mp.nstr(x,digits,strip_zeros=False) for x in u])
    return {'right':vectors[0],'left':vectors[1],'seed_digits':digits,'power_steps':steps,
            'purpose':'Exact decimal trial vectors; their accuracy is established only by the interval validator.'}


def verify(data,seeds,dps=100):
    t0=time.perf_counter();iv.dps=dps
    R=data['R'];D=1<<(R-1);p=R-1
    assert R>=2 and len(data['density'])==1<<R
    assert 0<F(data['c'])<1 and F(data['kappa'])>0
    assert data.get('model','periodic_soft_east_activity')=='periodic_soft_east_activity'
    assert data.get('density_convention','U=sum(density); no implicit s factor')=='U=sum(density); no implicit s factor'
    assert data.get('bit_order','leftmost_is_least_significant')=='leftmost_is_least_significant'
    te=edges(data,iv)
    emin=min(bounds(x)[0] for x in te);emax=max(bounds(x)[1] for x in te)
    assert emin>0
    # Each ordered pair of states has exactly one path of p=R-1 edges.
    smin=emin**p;smax=emax**p
    tau=(smax-smin)/(smax+smin)
    assert 0<=tau<1
    errors=[];rho_bounds=[];trials=[];residuals=[]
    for key,transpose in [('right',False),('left',True)]:
        v=list(map(F,seeds[key]));assert len(v)==D and min(v)>0
        assert v[0]==1
        vi=list(map(I,v));Tv=apply(vi,te,R,transpose)
        ratios=[bounds(Tv[j]/I(v[j])) for j in range(D)]
        rho_bounds.append((min(x[0] for x in ratios),max(x[1] for x in ratios)))
        Sv=vi
        for _ in range(p):Sv=apply(Sv,te,R,transpose)
        ratios=[bounds(Sv[j]/I(v[j])) for j in range(D)]
        lower=min(x[0] for x in ratios);upper=max(x[1] for x in ratios)
        assert lower>0 and upper>=lower
        hilbert=bounds(iv.log(I(upper)/I(lower)))[1]
        err=hilbert/(1-tau)
        errors.append(err);trials.append(vi);residuals.append(hilbert)
    rho=(max(x[0] for x in rho_bounds),min(x[1] for x in rho_bounds))
    assert 0<rho[0]<=rho[1]
    error=sum(errors,F(0))
    factor=iv.mpf([iv.exp(-I(error)).a,iv.exp(I(error)).b])
    right,left=trials
    phi=local_tables(data,True)[0];m=2*R-1
    weights=[];num=I(0);den=I(0);insert0=I(0)
    for word in range(1<<m):
        path=I(1)
        for j in range(R):path*=te[(word>>j)&((1<<R)-1)]
        base=path*left[word&(D-1)]*right[word>>R]
        weight=base*factor
        ph=iv.mpf([I(phi[word][0]).a,I(phi[word][1]).b])
        den+=weight;num+=weight*ph
        insert0+=base
    assert bounds(den)[0]>0
    Q=bounds(num/den)
    return {'status':'VERIFIED_THERMODYNAMIC_CANDIDATE_CENTER','interval_dps':dps,'R':R,'states':D,
            'positive_power':p,'positive_power_entry_lower_exact':str(smin),'positive_power_entry_upper_exact':str(smax),
            'tau_upper':dec(tau,True),'tau_upper_exact':str(tau),
            'right_Hilbert_error_upper_exact':str(errors[0]),'left_Hilbert_error_upper_exact':str(errors[1]),
            'right_Hilbert_error_upper_float':float(errors[0]),'left_Hilbert_error_upper_float':float(errors[1]),
            'perron_root':rec(*rho),'Q_infinity':rec(*Q),'Q_raw_width_exact':str(Q[1]-Q[0]),
            'Q_raw_width_float':float(Q[1]-Q[0]),'seconds':time.perf_counter()-t0,
            'scope':'True infinite-chain expectation for the fixed range-R candidate. No finite-size extrapolation. Does not by itself enclose the tilted generator spectral density.'}


def exact_product_check(dps):
    data={'R':2,'c':'0.3','kappa':'0.1','s':'-0.03','density':['0','0','0.04','0.04']}
    se=seed(data);r=verify(data,se,dps)
    c=I(F(3,10));k=I(F(1,10));eta=I(F(1,25));s=I(F(-3,100))
    z=1-c+c*iv.exp(2*eta);pr=c*iv.exp(2*eta)/z
    exact=(k+pr)*(c*(1-pr)*(iv.exp(-s+eta)-1)+(1-c)*pr*(iv.exp(-s-eta)-1))
    exactlo,exacthi=bounds(exact);q=r['Q_infinity']
    assert F(q['lower'])<=exacthi and F(q['upper'])>=exactlo
    return {'status':'PASS_EXACT_PRODUCT_FORMULA','exact':rec(exactlo,exacthi),'verified':q,
            'perron_root_exact':rec(*bounds(z)),'perron_root_verified':r['perron_root']}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--dps',type=int,default=100)
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--seeds',type=Path,required=True)
    ap.add_argument('--constants',type=Path,required=True);ap.add_argument('--generate-seeds',action='store_true')
    args=ap.parse_args();out=args.output.resolve();sp=args.seeds.resolve()
    assert out.is_relative_to(ROOT) and sp.is_relative_to(ROOT)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use a new result directory')
    out.mkdir(parents=True,exist_ok=True);sp.mkdir(parents=True,exist_ok=True)
    (out/'product_crosscheck.json').write_text(json.dumps(exact_product_check(args.dps),indent=2)+'\n')
    rows=[]
    for constpath in sorted(args.constants.glob('*_constants.json')):
        constants=json.loads(constpath.read_text());path=candidate_path(constants['candidate'])
        assert constants['status']=='INTERVAL_CONSTANTS_VERIFIED'
        assert 0<=F(constants['alpha']['upper'])<1 and F(constants['beta']['lower'])>0
        assert F(constants['d']['upper'])>=0 and F(constants['correction_upper'])>=0
        data=json.loads(path.read_text());sha=hashlib.sha256(path.read_bytes()).hexdigest()
        assert sha==constants['candidate_sha256']
        sf=sp/(path.stem+'_seed.json')
        if args.generate_seeds and not sf.exists():
            s=seed(data);s['candidate_sha256']=sha
            sf.write_text(json.dumps(s,indent=2)+'\n')
        seeds=json.loads(sf.read_text());assert seeds['candidate_sha256']==sha
        result=verify(data,seeds,args.dps)
        result.update(candidate=path.name,candidate_sha256=sha,seed_sha256=hashlib.sha256(sf.read_bytes()).hexdigest(),
                      kappa=data['kappa'],s=data['s'],reference_N=data['reference_N'])
        qlo=F(result['Q_infinity']['lower']);qhi=F(result['Q_infinity']['upper'])
        err=F(constants['correction_upper']);interval=rec(qlo,qhi+err)
        # Use the published decimal interval width, including all endpoint rounding.
        width=F(interval['upper'])-F(interval['lower'])
        old=json.loads(path.with_name(path.stem+'_verified.json').read_text())
        oldlo=F(old['lower_density']);oldhi=F(old['upper_density'])
        assert max(qlo,oldlo)<=min(qhi+err,oldhi)
        result.update(correction_upper=constants['correction_upper'],spectral_density_limit=interval,
                      spectral_width_upper=dec(width,True),old_common_graph_width=dec(oldhi-oldlo,True),
                      old_width_over_new_width=float((oldhi-oldlo)/width),
                      spectral_scope='Thermodynamic spectral density Lambda=lim lambda_n/n; uses the quadratic theorem, uniform constants, and the almost-additive finite-range limit proof.')
        (out/(path.stem+'_Q_infinity.json')).write_text(json.dumps(result,indent=2)+'\n')
        rows.append(result)
        print(json.dumps({k:result[k] for k in ['candidate','Q_raw_width_float','spectral_width_upper','seconds']}),flush=True)
    (out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    (out/'runtime.json').write_text(json.dumps({'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'interval_dps':args.dps,'method':'Positive-power Hilbert residual enclosure with exact decimal seeds'},indent=2)+'\n')


if __name__=='__main__':main()
