"""Original full three squares, all tau-linear terms, integral range bound."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
orig=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
own=load(str(p/'independent_even_cubic_RR_linear_module.sobj'))
new=load(str(p/'independent_delta_even_L_module.sobj'))
K=orig['ring'].base_ring();a,b,h,z=K.gens();r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q
NK=new['field'];back=NK.ring().hom((a,r,H,Z),K)
def backcoef(c):
    c=NK(c);return back(c.numerator())/back(c.denominator())
R=PolynomialRing(K,names=tuple('ir_p%d'%i for i in range(6))+tuple('ir_v%d'%i for i in range(3))+('ir_tau',));x=R.gens();pp=x[:6];v=x[6:9];tau=x[9]
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));ell=vector(R,(pp[0]*v[0]-alpha*pp[5]*v[2],pp[1]*v[0]-alpha/beta*pp[3]*v[1],pp[2]*v[1]-beta*pp[4]*v[2]))
unit={1:2*Z**2/((Z-H)**2*(Z-1)**2),5:2*H**2*Z**2/((Z-H)**2*(Z-1)**2),6:2*H**2*Z**2/((Z-H)**2*(H+1)**2)}
idx={1:2,5:0,6:1};oldguard=load(str(p/'independent_ratio_weighted_degreefall.sobj'))['oldguard_numerator'];oldguard*=oldguard.parent()(a*b*(a-1)*(b-1)*(a-b));items=[]
def copyold(f,n=6):return R({tuple(map(int,e))+(0,)*(10-n):K(c) for e,c in f.dict().items() if c})
def copynew(f,n=9):return R({tuple(map(int,e))+(0,)*(10-n):backcoef(c) for e,c in f.dict().items() if c})
for qo,no in zip(own['items'],new['items']):
    qi=int(qo['source']);ck('source_order_%s'%qi,qi==no['Q_index'])
    A=[copyold(f) for f in qo['polynomial_cofactors']]
    AA=[copynew(f,6) for f in no['A_p_quadratic']]
    B=[copynew(f) for f in no['B_RR_linear']]
    for j in range(3):ck('independent_two_p3RR_cofactor_agreement_%d_%d'%(qi,j),A[j]==AA[j])
    cf=[A[j]+B[j] for j in range(3)]
    raw=orig['quartic_source_rows'][qi];hd={};bd={};td={};ntop={}
    for ee,c in raw.dict().items():
        e=tuple(map(int,ee));pd=sum(e[:6]);rd=sum(e[6:]);power=6-pd-2*rd
        ck('all_original_weights_%d_%s'%(qi,e),power>=0 and pd+rd<=4)
        hd[e+(power,)]=K(c)
        if power>=2:bd[e+(power-2,)]=K(c)
        elif power==1:td[e+(0,)]=K(c)
        else:ntop[e+(0,)]=K(c)
    hom=R(hd);tail=R(bd);tlinear=R(td)
    ck('original_top_exact_square_%d'%qi,R(ntop)==unit[qi]*ell[idx[qi]]**2)
    ck('ALL_original_tau_linear_exact_module_%d'%qi,tlinear==sum((cf[j]*ell[j] for j in range(3)),R.zero()))
    ck('FULL_original_identity_%d'%qi,hom==unit[qi]*ell[idx[qi]]**2+tau*sum((cf[j]*ell[j] for j in range(3)),R.zero())+tau**2*tail)
    for j,F in enumerate(cf):
        for c in F.coefficients():
            for fac,power in K(c).denominator().factor():
                ck('old_only_coefficient_den_%d_%d_%s'%(qi,j,fac),oldguard.quo_rem(fac)[1]==0)
    for fac,power in unit[qi].numerator().factor():ck('square_unit_num_%d_%s'%(qi,fac),oldguard.quo_rem(fac)[1]==0)
    for fac,power in unit[qi].denominator().factor():ck('square_unit_den_%d_%s'%(qi,fac),oldguard.quo_rem(fac)[1]==0)
    items.append({'source':qi,'square_index':idx[qi],'square_unit':unit[qi],'homogeneous_source':hom,'tau_linear_coefficients':cf,'tau_squared_tail':tail,'monic_linear_coefficients':[f/unit[qi] for f in cf],'monic_constant':tail/unit[qi]})
    print('FULL ORIGINAL Q',qi,'PASS',round(time.monotonic()-st,2),flush=True)
save({'ring':R,'coefficient_field':K,'ell':ell,'items':items,'checks':checks,'source_files':['round8/rank3_joint_fiber_leading.sobj','independent_even_cubic_RR_linear_module.sobj','independent_delta_even_L_module.sobj'],'scope':'Exact full weighted homogenizations Qh=U ell_j^2+tau sum Aji ell_i+tau^2 B. No p/RR/Delta/C division; e=ell/tau is integral and bounded on all valued arcs over old protected base.'},str(p/'independent_even_integral_range_reader.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'full_original_sources':[it['source'] for it in items],'tail_term_counts':[len(it['tau_squared_tail'].dict()) for it in items],'scope':'Integral range bound only, not RR/tau boundedness or no-lift on Delta=0.'}
(p/'independent_even_integral_range_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
