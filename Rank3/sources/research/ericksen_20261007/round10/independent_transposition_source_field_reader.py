"""Independent exact reader for full original transposition covariance.

Substitution is reconstructed using a polynomial homomorphism. All source,
base-field, and full state-field identities, involution, and unit guards
are checked without trusting the creator's booleans.
"""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'joint_fields_direct.sobj'));R=d['ring'];K=R.base_ring();a,b,h,z=K.gens();x=R.gens();A=K.ring()
cre=load(str(p/'independent_transposition_even_field_covariance.sobj'));fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
def ck(n,v):
 assert bool(v),n;checks.append(n)
pp=(2,3,0,1,5,4,7,6,8);dirs=(1,0,2);par=(1/a,b/a,1/h,z/h);ph=A.hom(par,K);cache={}
def cf(c):
 c=K(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def cp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
statehom=R.hom([x[j]/h for j in pp],R)
def pull(f):return statehom(R({tuple(map(int,e)):cf(c) for e,c in f.dict().items()}))
guard=A(a*b*(a-1)*(b-1)*(a-b)*h*z*(h-1)*(z-1)*(h-z));denominators=[]
def unit(n,f):
 f=K(f);ck(n+' nonzero',f!=0)
 for mode,g in (('numerator',f.numerator()),('denominator',f.denominator())):
  for v,mul in g.factor():
   if not v.is_constant():ck(n+' guarded '+mode+' '+str(v),guard%v==0)
 denominators.append((n,f))
for j,v in enumerate(par):ck('parameter involution '+str(j),cf(v)==K.gens()[j])
for j,v in enumerate(x):ck('state involution '+str(j),pull(pull(v))==v)
for j,g in enumerate((a,b,a-1,b-1,a-b,h,z,h-1,z-1,h-z)):unit('guard image '+str(j),cf(g))
Q0=cp(fd['Q0']);qfactor=-a*(b-1)/(a-b)/h**4;unit('Q0 unit',qfactor);ck('original Q0 transport',pull(Q0)==qfactor*Q0)
Rs=[cp(it['full']) for it in fd['items']];rpp=(1,0,3,2,5,4);rf=(h**-4,h**-4,qfactor,qfactor,-qfactor,-qfactor)
for j,(jj,v) in enumerate(zip(rpp,rf)):unit('R unit '+str(j),v);ck('original complete R transport '+str(j),pull(Rs[j])==v*Rs[jj])
Qs=[cp(oq['quartic_source_rows'][qi]) for qi in (1,5,6)];eqco=(K.zero(),K.zero(),b*(a-1)/(a-b)/h**4)
for j,jj in enumerate((1,0,2)):ck('original complete even transport '+str(j),pull(Qs[j])==h**-4*Qs[jj]+eqco[j]*Q0)
F=d['fields'];B=d['base_fields']
for i,ii in enumerate(dirs):
 ck('complete base h transport '+str(i),pull(B[i,0])==B[ii,0]/h**3)
 ck('complete base z transport '+str(i),pull(B[i,1])==-B[ii,1]/h**2+z*B[ii,0]/h**3)
 for j,jj in enumerate(pp):
  co=next(it['Q0_cofactor'] for it in cre['field_Q0_cofactors'] if it['direction']==i and it['variable']==j)
  ck('complete state transport '+str((i,j)),pull(F[i,j])==-F[ii,jj]/h**2+x[jj]*B[ii,0]/h**3+co*Q0)
print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
save({'ring':R,'parameter_images':par,'state_images':[x[j]/h for j in pp],'direction_permutation':dirs,'source_R_permutation':rpp,'source_R_units':rf,'Q0_unit':qfactor,'even_Q0_cofactors':eqco,'unit_checks':denominators,'checks':checks,'scope':__doc__},str(p/'independent_transposition_source_field_reader.sobj'),compress=False)
(p/'independent_transposition_source_field_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n')
