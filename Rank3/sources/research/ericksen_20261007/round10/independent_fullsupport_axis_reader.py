"""Original highest sources through all eight integer guard certificates."""
from pathlib import Path
from itertools import product
import time,json
from sage.all import load,save,QQ,ZZ,PolynomialRing,prod
out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
d=load(str(out/'general_curve_two_pair_offdiagonal.sobj'))
lead=load(str(out/'general_curve_dimensionless_leading.sobj'))
probe=load(str(out/'general_curve_fullsupport_axis_probe.sobj'))
elim=load(str(out/'general_curve_fullsupport_axis_certificates.sobj'))
portable=json.loads((out/'general_curve_fullsupport_axis_integer.json').read_text())
R=d['ring'];K=R.base_ring();H,Z,S,N=K.gens();im=d['images']
C=PolynomialRing(QQ,names=('H','Z','S'));h,z,s=C.gens();CK=C.fraction_field()
guards=[h,z,h-z,z-1,h+1,h-z+1,s-1,h*s-z+1,(h+1)*s-z];guard=prod(guards)
def unit(c):
    if not c:return False
    for p in [c.numerator(),c.denominator()]:
        for f,e in C(p).factor():
            f=C(f)
            if f.is_constant():continue
            if not any(all((CK(f)/g).derivative(v)==0 for v in C.gens()) for g in guards):return False
    return True
def cc(co):return K(co(0,0,H,Z,S,N))
mixed=[]
for (i,j),block in lead['blocks'].items():mixed.append(sum(cc(block[u,v])*im[2*i+u]*im[2*j+v] for u in range(2) for v in range(2)))
trilinear=[]
for coeffs in lead['cubic_rows']:
    trilinear.append(sum(cc(co)*im[u]*im[2+v]*im[4+w] for co,(u,v,w) in zip(coeffs,product(range(2),repeat=3))))
def from_terms(ts):
    assert all(type(c) is int and len(m)==3 and all(type(e) is int and e>=0 for e in m) for m,c in ts)
    return C({tuple(m):ZZ(c) for m,c in ts})
seen=[]
for i,bits in enumerate(product((0,1),repeat=3)):
    tag=''.join(map(str,bits));seen.append(bits)
    point=[K.zero()]*6
    for j,b in enumerate(bits):point[2*j+b]=K.one()
    rows=[K(f(*point)) for f in mixed+trilinear]
    pc=next(v for v in probe['cases'] if tuple(v['axes'])==bits)
    ec=next(v for v in elim['cases'] if tuple(v['axes'])==bits)
    js=next(v for v in portable['cases'] if tuple(v['axes'])==bits)
    checks[tag+'_all_original_leaders']=rows==pc['sources'] and len(rows)==7
    aa=rows[0].derivative(N);bb=rows[0].subs({N:0})
    checks[tag+'_linear_N']=rows[0]==aa*N+bb and not aa.derivative(N)
    ap=CK(aa.numerator()(h,z,s,0))/CK(aa.denominator()(h,z,s,0))
    bp=CK(bb.numerator()(h,z,s,0))/CK(bb.denominator()(h,z,s,0))
    nval=-bp/ap
    checks[tag+'_protected_N_pivot']=unit(ap)
    checks[tag+'_original_N_graph']=nval==ec['N_graph'] and unit(nval.denominator())
    source=[]
    for j,row in enumerate(rows[1:]):
        val=CK(row.numerator()(h,z,s,nval))/CK(row.denominator()(h,z,s,nval))
        if not val:continue
        p=ec['sources'][len(source)];source.append(p)
        checks[tag+'_primitive_original_'+str(j)]=unit(val/CK(p))
    checks[tag+'_complete_six_sources']=len(source)==len(ec['sources'])==len(js['sources'])==6
    ints=[from_terms(p) for p in js['sources']];multip=[from_terms(p) for p in js['multipliers']]
    for j,(integer,p) in enumerate(zip(ints,source)):
        ratio=CK(integer/p)
        checks[tag+'_integer_source_link_'+str(j)]=bool(ratio) and all(ratio.derivative(v)==0 for v in C.gens())
    target=from_terms(js['target']);D=ZZ(js['integer_multiplier'])
    checks[tag+'_original_guard']=ec['guard']==guard
    checks[tag+'_integer_target']=D>0 and target==D*guard
    checks[tag+'_integer_identity']=sum((aa*bb for aa,bb in zip(ints,multip)),C.zero())==target
    E=ec['local_ring']
    checks[tag+'_localized_identity']=sum((aa*bb for aa,bb in zip(ec['multipliers'],ec['localized_sources'])),E.zero())==1
checks['all_eight_axes']=len(seen)==8 and len(portable['cases'])==8 and len(elim['cases'])==8 and elim['complete']
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'scope':'Pointwise exclusion of every full-support triple of original opposite-pair coordinate-axis directions on all nine guards, from original mixed/trilinear highest sources and explicit integer guard identities. Not a full-support amplitude dimension result.'},str(out/'independent_fullsupport_axis_reader.sobj'))
report={'status':'EXACT_EIGHT_FULLSUPPORT_AXIS_SOURCE_CERTIFICATES_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_fullsupport_axis_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
