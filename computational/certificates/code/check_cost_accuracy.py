"""Check main-table accuracy/cost claims against frozen endpoint and timing records."""
from pathlib import Path
from fractions import Fraction as F
import csv,hashlib,json,statistics
P=Path(__file__).resolve().parents[2]
def read(s):return json.loads((P/s).read_text(encoding="utf-8-sig"))
def main():
    if not __debug__:raise RuntimeError("Do not disable assertions")
    manifest=read("MANIFEST.json")
    for name,item in manifest.items():
        path=(P/name).resolve();assert path.is_relative_to(P.resolve())
        assert hashlib.sha256(path.read_bytes()).hexdigest()==item["sha256"],name
    canonical={r["candidate"]:r for r in read("certificates/published/certificates.json")}
    timings={r["candidate"]:r for r in read("certificates/published/range_study.json")}
    rows=read("certificates/published/table1_cost_accuracy.json")
    endpoints=read("certificates/published/table1_endpoints.json")
    with (P/"certificates/published/table1_endpoints.csv").open(newline="",encoding="utf-8-sig") as f:csv_endpoints=list(csv.DictReader(f))
    with (P/"certificates/published/table1_cost_accuracy.csv").open(newline="",encoding="utf-8-sig") as f:csv_costs=list(csv.DictReader(f))
    names={r["candidate"] for r in rows}
    assert len(rows)==len(names)==7
    assert names=={r["candidate"] for r in endpoints}=={r["candidate"] for r in csv_endpoints}=={r["candidate"] for r in csv_costs}
    for r in rows:
        x=canonical[r["candidate"]];assert r["sha256"]==x["sha256"]
        for k in ["c","kappa","s","N","R"]:assert r[k]==x[k]
        width=F(x["graph"]["upper_density"])-F(x["graph"]["lower_density"])
        assert F(r["graph_width_exact"])==width
        q=F(x["infinite"]["spectral_width_upper"]) if x["quadratic"] else None
        assert (F(r["quadratic_width_exact"]) if r["quadratic_width_exact"] else None)==q
        source=timings.get(r["candidate"])
        for method in ["graph","quadratic"]:
            expected=statistics.median(source["timing"][method+"_seconds"]) if source else None
            if source:assert len(source["timing"][method+"_seconds"])==3
            assert r[method+"_seconds_median"]==expected
        e=next(v for v in csv_endpoints if v["candidate"]==r["candidate"])
        assert e["candidate_sha256"]==x["sha256"]
        for k in ["lower","upper"]:assert e["graph_"+k]==x["graph"][k+"_density"]
        if x["quadratic"]:
            for k in ["lower","upper"]:
                assert e["Q_infinity_"+k]==x["infinite"]["Q_infinity"][k]
                assert e["Lambda_"+k]==x["infinite"]["spectral_density_limit"][k]
    moderate=sorted([x for x in timings.values() if x["kappa"]=="0.1" and x["N"]==18],key=lambda x:x["R"])
    target=F(1,10**6)
    first_graph=next(x for x in moderate if F(x["graph"]["upper_density"])-F(x["graph"]["lower_density"])<=target)
    first_quadratic=next(x for x in moderate if F(x["infinite"]["spectral_width_upper"])<=target)
    assert (first_graph["R"],first_quadratic["R"])==(9,7)
    r9=next(x for x in rows if x["kappa"]=="0.1" and x["R"]==9)
    ratio=float(F(r9["graph_width_exact"])/F(r9["quadratic_width_exact"]))
    cost=r9["quadratic_seconds_median"]/r9["graph_seconds_median"]
    assert round(ratio)==1994 and round(cost,2)==9.48
    report=dict(status="PASS",table_rows=7,quadratic_rows=5,first_reported_graph_range_at_1e_6=9,first_reported_quadratic_range_at_1e_6=7,width_ratio_R9=ratio,verification_time_ratio_R9=cost,timing_scope="original three-run medians; no candidate/seed generation; no new benchmark")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
