"""Spectrum-uniform generic-base kernel, with ALL u/v tangent values kept."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p12_axis_generic_normal.sobj'));el=load(str(p/'independent_ratio_p12_axis_tau_elimination.sobj'));re=load(str(p/'independent_ratio_axis_tau_univariate.sobj'));M=d['matrix'];U=d['ring'];K=d['coefficient_field'];u,v=U.gens()
ceta=K(M[0,5].monomial_coefficient(u));assert M[0,5]==ceta*u
eta_p32=U({(int(e[0])-1,int(e[1])):c for e,c in M[0,4].dict().items()});assert M[0,4]==eta_p32*u
lam=-eta_p32/ceta;beta=-M[13,3]/K(M[13,0]);assert lam in U and beta in U
vec1=vector(U,[beta,0,0,1,0,0,0]);vec2=vector(U,[0,0,0,0,1,lam,0]);vec3=vector(U,[0,0,0,0,0,1,0])
assert M*vec1==0 and M*vec2==0
M0=M.apply_map(lambda f:U(f.subs({u:0})));assert M0*vec3==0
Ks=re['records'][0]['cyclic_determinant'].parent();a,r=Ks.gens();ph=K.ring().hom([a,a*(1+r)/(1+a*r),1/r,-(1+r)/r],Ks)
def cyc(c):return ph(c.numerator())/ph(c.denominator())
base_factors={'cu':el['cu'],'A_p23':el['A_coefficient'],'Q_u2':re['u2_coefficient'],'Q1_p21':K(M[1,1]),'H0_p13':K(M[13,0]),'Q0_RR1_over_u':ceta}
certs=[]
for name,c in base_factors.items():
 val=cyc(c);assert val
 fac=val.factor();certs.append({'name':name,'factor':c,'cyclic_value':val,'cyclic_factorization':fac});print(name,str(fac),flush=True)
for it in re['records']:
 assert it['cyclic_determinant'];print('RESULTANT',it['chart'],it['cyclic_determinant'].factor(),flush=True)
save({'ring':U,'coefficient_field':K,'normal_positions':d['normal_positions'],'kernel_u_nonzero':(vec1,vec2),'kernel_u_zero':(vec1.subs({u:0}),vec2.subs({u:0}),vec3),'beta':beta,'lambda':lam,'base_nonzero_certificates':certs,'resultant_sources':re['records'],'source_files':['independent_ratio_p12_axis_generic_normal.sobj','independent_ratio_p12_axis_tau_elimination.sobj','independent_ratio_axis_tau_univariate.sobj'],'scope':__doc__,'qualification':'Rank5 for u!=0 and rank4 for u=0 on a material-uniform generic base open. Every kernel has tau=0, but higher-order ramified lifts have NOT been excluded.'},str(p/'independent_ratio_axis_normal_kernel_certificate.sobj'),compress=False)
(p/'independent_ratio_axis_normal_kernel_certificate.json').write_text(json.dumps({'checks':{'two_symbolic_kernel_vectors':True,'additional_u_zero_vector':True,'base_factors_nonzero_for_each_material':True,'two_cyclic_resultant_units':True},'rank_u_nonzero':5,'rank_u_zero':4,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
