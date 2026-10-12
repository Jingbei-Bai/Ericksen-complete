"""Independent sparse Laurent arithmetic and genuine four-row source audit.

No LaurentSeriesRing, polynomial quotient, or polynomial remainder pipeline
is used. Every series carries an absolute Big-O bound; the quadratic
extension is represented by a pair and explicit multiplication.
"""
from sage.all import *
from pathlib import Path
from itertools import permutations
import time,json,os,gc
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(label,v):
    if not bool(v):raise AssertionError(label)
    checks.append(label)
R=PolynomialRing(QQ,names=('lr_a','lr_b'));a,b=R.gens();K=R.fraction_field();INF=10**8
inverted=[]
class Jet:
    def __init__(self,terms=None,prec=INF):
        self.prec=int(prec);self.d={int(e):K(c) for e,c in (terms or {}).items() if c and int(e)<self.prec}
    @staticmethod
    def constant(c):return Jet({0:K(c)}) if c else Jet()
    def val(self):return min(self.d) if self.d else self.prec
    def __neg__(self):return Jet({e:-c for e,c in self.d.items()},self.prec)
    def __add__(self,o):
        if not isinstance(o,Jet):o=Jet.constant(o)
        pr=min(self.prec,o.prec);dd=dict(self.d)
        for e,c in o.d.items():dd[e]=dd.get(e,K.zero())+c
        return Jet(dd,pr)
    __radd__=__add__
    def __sub__(self,o):return self+(-o if isinstance(o,Jet) else -K(o))
    def __rsub__(self,o):return (-self)+o
    def __mul__(self,o):
        if not isinstance(o,Jet):o=Jet.constant(o)
        if (not self.d and self.prec==INF) or (not o.d and o.prec==INF):return Jet()
        pr=min(self.prec+o.val(),o.prec+self.val(),INF);dd={}
        for e,c in self.d.items():
            for f,d in o.d.items():
                if e+f<pr:dd[e+f]=dd.get(e+f,K.zero())+c*d
        return Jet(dd,pr)
    __rmul__=__mul__
    def inv(self):
        if not self.d:raise ValueError('cannot invert an unknown/zero series')
        v=self.val();c=self.d[v];inverted.append(c)
        if self.prec==INF and len(self.d)==1:return Jet({-v:1/c})
        length=self.prec-v
        if length>=INF//2:raise ValueError('unsupported infinite rational inverse')
        gg=[1/c]
        for n in range(1,length):gg.append(-sum((self.d.get(v+j,K.zero())*gg[n-j] for j in range(1,n+1)),K.zero())/c)
        return Jet({-v+n:cc for n,cc in enumerate(gg) if cc},self.prec-2*v)
    def __truediv__(self,o):return self*(o.inv() if isinstance(o,Jet) else Jet.constant(o).inv())
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        ans=Jet.constant(1);v=self
        while n:
            if n&1:ans=ans*v
            n//=2
            if n:v=v*v
        return ans
    def record(self):return {'terms':self.d,'precision_absolute':self.prec,'valuation_lower_bound':self.val()}
def eqjet(A,B,prec):
    ck('precision sufficient '+str(len(checks)),A.prec>=prec and B.prec>=prec)
    return all(A.d.get(i,0)==B.d.get(i,0) for i in set(A.d)|set(B.d) if i<prec)
poly_groups={};evcache={}
def ray_polynomial(f,sign,precision):
    key=(f,sign)
    if key not in poly_groups:
        grouped={}
        for ex,c in f.dict().items():
            i,j,h,z=map(int,ex);power=h+z;dd=grouped.setdefault(power,{})
            dd[(i,j)]=dd.get((i,j),QQ.zero())+QQ(c)*sign**z
        poly_groups[key]={n:R(dd) for n,dd in grouped.items() if R(dd)}
    dd=poly_groups[key]
    if not dd:return Jet()
    degree=max(dd);return Jet({-n:K(c) for n,c in dd.items() if n>degree-precision},-degree+precision)
def ev(c,sign,precision):
    key=(c,sign,precision)
    if key not in evcache:
        evcache[key]=ray_polynomial(c.numerator(),sign,precision)/ray_polynomial(c.denominator(),sign,precision)
    return evcache[key]
def add2(u,v):return (u[0]+v[0],u[1]+v[1])
def neg2(u):return (-u[0],-u[1])
def mul2(u,v,p0,p1):return (u[0]*v[0]-p0*u[1]*v[1],u[0]*v[1]+u[1]*v[0]-p1*u[1]*v[1])
def poly_pair(coeffs,p0,p1):
    ans=(Jet(),Jet())
    for c in reversed(coeffs):ans=add2(mul2(ans,(Jet(),Jet.constant(1)),p0,p1),(c,Jet()))
    return ans
def exact_pair(coeffs,p0,p1):
    ans=(p0.parent().zero(),p0.parent().zero())
    for c in reversed(coeffs):ans=(-p0*ans[1]+c,ans[0]-p1*ans[1])
    return ans
def determinant3(rows,p0,p1):
    ans=(Jet(),Jet())
    for perm in permutations(range(3)):
        parity=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3));term=(Jet.constant(parity),Jet())
        for i in range(3):term=mul2(term,rows[i][perm[i]],p0,p1)
        ans=add2(ans,term)
    return ans

co=load(str(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.sobj'));OK=co['ring'].base_ring();dd=-co['relation'][0]
cm={int(it['source_row']):it for it in co['items']};basis={rr:load(str(p/('independent_ratio_incoming12_all_third_raw.row%d.sobj'%rr))) for rr in (1,17,5)}
streams={};rawrecords={};source_input_fractions={dd};zero8=(OK.zero(),)*8
def sqrtmul(cf,vv):
    out=[]
    for i in range(4):out.extend((cf[0]*vv[2*i]+dd*cf[1]*vv[2*i+1],cf[0]*vv[2*i+1]+cf[1]*vv[2*i]))
    return tuple(out)
def rawcombine(terms,ratios,key,slot=None):
    ans=[OK.zero() for unused in range(8)]
    for term,cf in zip(terms,ratios):
        vec=term[key] if slot is None else term[key][slot]
        val=sqrtmul(cf,vec)
        ans=[x+y for x,y in zip(ans,val)]
    return tuple(ans)
for row in (0,7,11,13):
    d=load(str(p/('independent_ratio_incoming12_third_stream_cokernel.row%d.sobj'%row)));rr=load(str(p/('independent_ratio_incoming12_all_third_raw.row%d.sobj'%row)));terms=[rr,basis[1],basis[17],basis[5]];ratios=d['source_range_ratios'];cf=cm[row]['cofactors'];parity=int(d['parity']);PP=d['quadratic'];p0=PP[0]/PP[2];p1=PP[1]/PP[2]
    ck('original range rows '+str(row),len(ratios)==4)
    for j,u in enumerate(ratios):
        prod=(u[0]*cf[0][0]+dd*u[1]*cf[0][1],u[0]*cf[0][1]+u[1]*cf[0][0]);ck('original range cofactor ratio '+str((row,j)),prod==(cf[j][0],cf[j][1]));source_input_fractions.update(u)
    for col in range(6):ck('all third range corrections vanish '+str((row,col)),rawcombine(terms,ratios,'third_range_coefficients',col)==zero8)
    for slot in (0,3):
        vec=rawcombine(terms,ratios,'unknown_coefficients',slot)
        ck('eliminated source foreign l parity '+str((row,slot)),all(vec[2*i+1-parity]==0 for i in range(4)))
        reduced=exact_pair([vec[2*i+parity] for i in range(4)],p0,p1)
        ck('genuine complete k2 du2 columns zero modulo P '+str((row,slot)),reduced==(OK.zero(),OK.zero()))
        ck('stored zero columns '+str((row,slot)),d['coefficients'][slot]==0)
    # A parity check uses exact original arithmetic for the remaining columns;
    # coefficients themselves are retained as four-term source combinations
    # for the separate Laurent reconstruction, avoiding huge generic sums.
    for slot in (1,2,4,5):
        key='constant' if slot==5 else 'unknown_coefficients';ix=None if slot==5 else slot
        vec=rawcombine(terms,ratios,key,ix)
        if slot in (1,2):vec=sqrtmul((OK.zero(),OK.one()),vec)
        ck('complete source foreign l parity zero '+str((row,slot)),all(vec[2*i+1-parity]==0 for i in range(4)))
        # Restoring row_pivot is checked BEFORE quadratic reduction on the
        # already stripped common l^parity. This is an exact K identity.
        pol=d['unreduced_coefficients'][slot]
        ck('raw source and restored row pivot '+str((row,slot)),all(vec[2*i+parity]==d['row_pivot']*pol[i] for i in range(4)))
        qs=d['quotients'][slot];red=d['coefficients'][slot]
        ck('exact stored polynomial division '+str((row,slot)),all(pol[i]==sum((qs[j]*PP[i-j] for j in range(max(0,i-2),min(i,qs.degree())+1)),OK.zero())+red[i] for i in range(max(pol.degree(),qs.degree()+2,red.degree())+1)))
    streams[row]={'coefficients':[list(f) for f in d['coefficients']],'row_pivot':d['row_pivot'],'quadratic':list(PP),'parity':parity}
    print('FULL RAW SOURCE',row,'PASS seconds',round(time.monotonic()-st,2),flush=True)
    del rr,terms,d;gc.collect()

results=[];all_inverse_factors={}
for sign,cycle in ((-1,0),(1,1),(-1,2)):
    evcache.clear()  # Every cycle audits its own complete inverse list.
    begin=len(inverted);precision=1;raw={};Pcoeff=streams[0]['quadratic'];Pseries=[ev(c,sign,precision) for c in Pcoeff];p0=Pseries[0]/Pseries[2];p1=Pseries[1]/Pseries[2]
    # This explicit inverse also records that the stripped sqrt(l)^parity
    # is a unit on both embeddings, and audits actual source-range inverses.
    ev(dd,sign,precision).inv()
    for c in source_input_fractions:ev(c,sign,precision)
    for row,d in streams.items():
        ck('all same actual quadratic '+str((sign,cycle,row)),d['quadratic']==Pcoeff)
        pivot=ev(d['row_pivot'],sign,precision)
        raw[row]=[tuple(v*pivot for v in poly_pair([ev(c,sign,precision) for c in f],p0,p1)) for f in d['coefficients']]
    ref=raw[7];odd=[]
    for row in (0,11,13):odd.append([add2(mul2(ref[2],raw[row][j],p0,p1),neg2(mul2(raw[row][2],ref[j],p0,p1))) for j in (1,4,5)])
    f=determinant3(odd,p0,p1);norm=f[0]**2-p1*f[0]*f[1]+p0*f[1]**2
    ck('strict known leading order '+str((sign,cycle)),norm.val()==-30 and norm.prec>=-29)
    lead=norm.d[-30]
    expected=load(str(p/('independent_incoming_third_base_infinity_probe.sobj' if sign<0 else 'independent_incoming_third_positive_ray_probe.sobj')))
    ER=expected['material_ring'];mp=ER.hom([a,b],R)
    def oldmap(c):return K(mp(c.numerator()))/K(mp(c.denominator()))
    ck('independent exact leading agrees '+str((sign,cycle)),lead==oldmap(expected['leading_coefficient']))
    ck('independent exact absolute BigO '+str((sign,cycle)),norm.prec==int(expected['source_norm'].precision_absolute()))
    material_values=[lead]+inverted[begin:];factors=set()
    for value in material_values:
        for pol in (K(value).numerator(),K(value).denominator()):
            factors.update(factor for factor,mul in pol.factor())
    S=PolynomialRing(QQ,names=('positive_s','positive_t'));ss,tt=S.gens();shift=R.hom([1+ss,1+ss+tt],S)
    trans=R.hom(([K(a),K(b)],[K(b/a),K(1/a)],[K(1/b),K(a/b)])[cycle],K);factor_records=[]
    C=2352*a**3*b**2+961*a**2*b**3-672*a**3*b-3970*a**2*b**2-1984*a*b**3+48*a**3+577*a**2*b+1664*a*b**2+1024*b**3
    decomposition=3313*ss**2*(ss+1)**3+ss*tt*(ss+1)**2*(7587*ss+610)+tt**2*(ss+1)*(5235*ss**2+548*ss+49)+(31*ss-1)**2*tt**3
    ck('exceptional C exact positive decomposition '+str((sign,cycle)),shift(C)==decomposition)
    for fac in sorted(factors,key=lambda f:(f.total_degree(),str(f))):
        ff=trans(fac);num=shift(ff.numerator());den=shift(ff.denominator())
        simple=bool(num) and (all(c>0 for c in num.coefficients()) or all(c<0 for c in num.coefficients()))
        special=(cycle==0 and K(fac)/K(C) in QQ and K(fac)/K(C)!=0)
        ck('every actual inverse and leading factor nonzero '+str((sign,cycle,fac)),(simple or special) and bool(den) and all(c>0 for c in den.coefficients()))
        factor_records.append({'factor':fac,'cyclic_numerator':num,'cyclic_denominator':den,'C_decomposition_used':special})
    rec={'ray_sign':sign,'cycle':cycle,'p0':p0.record(),'p1':p1.record(),'norm':norm.record(),'leading_coefficient':lead,'factors':factor_records,'inverse_count':len(inverted)-begin,'raw_rows':{row:[(u.record(),v.record()) for u,v in rr] for row,rr in raw.items()},'l_parities':{row:d['parity'] for row,d in streams.items()}}
    results.append(rec);print('RAY/CYCLE',sign,cycle,'valuation',norm.val(),'BigO',norm.prec,'FACTORS',len(factors),'PASS seconds',round(time.monotonic()-st,2),flush=True)
    save({'results':results,'checks':checks,'scope':__doc__},str(p/'independent_incoming_laurent_source_reader.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'rays_cycles':[(it['ray_sign'],it['cycle']) for it in results],'factor_counts':[len(it['factors']) for it in results],'scope':__doc__,'source_scope':'Four original raw third rows with all stored second-range corrections; exact range cancellation, parity, restored row pivots, exact P divisions, true zero columns. Determinant computed without a 2x2 pivot. All material inverses audited after each genuine cyclic relabeling.'}
(p/'independent_incoming_laurent_source_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
