"""First genuine H-derived subcritical initials on the full-amplitude line.

P_H=m det(C) H-A_H adj(C) W R cancels the (9,2) layer. For
1<gamma<3/2 the (10,1) initial is m[det(C) H41-A_H adj(C)d].
The stored line row is this initial divided by the protected m*A*B*C.
"""
from sage.all import *
from pathlib import Path
import time,json,os
p=Path(__file__).resolve().parent;st=time.monotonic();mod=load(str(p/'independent_delta_H_global_L_module.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));q=load(str(p/'independent_delta_subcritical_Q0_direction.sobj'))
old=mod['original_field'];pool=full['polynomial_pool'];B=PolynomialRing(QQ,names=('ha','hr','hH','hZ','ht'));K=B.fraction_field();a,r,H,Z,t=K.gens();x=-(t+1)/t;z=t
ph=old.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],K);cache={}
def cv(c):
 c=old(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
dh=q['field'].ring().hom([a,r,H,Z,x,z],K)
def dc(c):c=q['field'](c);return dh(c.numerator())/dh(c.denominator())
P=PolynomialRing(K,names=('hA','hB','hC'));A,Bb,C=P.gens();dr=tuple(P(w)*s for w,s in zip((x,x*z,1,z,1,1),(A,Bb,A,C,Bb,C)))
v=tuple(P(dc(w))*s for w,s in zip(q['RR_kernel'],(C,Bb,A)));det=dc(q['C_determinant']);en=vector(P,[P(dc(w))*s for w,s in zip(q['e_numerator'],(A*C,Bb*C,A*Bb))]);vals=dr+v
def face(dag,bi):
 out={}
 for sign,li,ri in dag:
  for ee,c in pool[li].items():
   e=tuple(map(int,ee))
   for ff,b in pool[ri].items():
    f=tuple(map(int,ff));g=tuple(e[j]+f[j] for j in range(9))
    if (sum(g[:6]),sum(g[6:]))==bi:out[g]=out.get(g,old.zero())+sign*c*b
 return {e:c for e,c in out.items() if c}
def pull(dd):return sum((P(cv(c))*prod(vals[j]**int(e[j]) for j in range(9)) for e,c in dd.items()),P.zero())
checkpoint=p/'independent_delta_subcritical_H_initials.sobj'
previous=load(str(checkpoint)) if checkpoint.exists() else None
checks=list(previous['checks']) if previous else [];items=list(previous['items']) if previous else []
def ck(label,yes):assert yes,label;checks.append(label)
for j,it in enumerate(mod['items']):
 if j<len(items):continue
 ck('global L module proved '+str(j),it['remainder']==0)
 rec=full['records'][j];ck('source indexing '+str(j),rec['direction']==it['direction'] and rec['source_position']==it['source_position'])
 ff=face(rec['full_signed_product_DAG'],(4,1));cf=[]
 for f in it['cofactors']:
  cf.append({tuple(map(int,e)):old(c)/prod(mod['original_to_binomial_variable_scales'][k]**int(e[k]) for k in range(9)) for e,c in f.dict().items()})
 Ap=vector(P,[pull(dd) for dd in cf]);hp=pull(ff);row=det*hp+Ap*en
 ck('full row degree five '+str(j),not row or all(sum(e)==5 for e in row.dict()))
 parity=set(tuple(int(v)%2 for v in e) for e in row.dict());ck('one common amplitude parity '+str(j),len(parity)<=1)
 pa=next(iter(parity)) if parity else (0,0,0);primitive=P({tuple(int(e[k])-pa[k] for k in range(3)):c for e,c in row.dict().items()})
 ck('amplitude factor exact '+str(j),row==primitive*prod(P.gens()[k]**pa[k] for k in range(3)))
 print('H',j,'face41terms',len(ff),'rowterms',len(row.dict()),'parity',pa,'coeffterms',[len(c.numerator().dict()) for c in row.coefficients()],'seconds',time.monotonic()-st,flush=True)
 items.append({'row':j,'direction':it['direction'],'source_position':it['source_position'],'full_source_DAG':rec['full_signed_product_DAG'],'H41_original':ff,'L_cofactors_original':cf,'full_line_initial_divided_mABC':row,'amplitude_parity':pa,'amplitude_reduced_row':primitive})
 tmp=p/'independent_delta_subcritical_H_initials.tmp.sobj';save({'source_files':['independent_delta_H_global_L_module.sobj','independent_ratio_global_R_derivative_degreefall.sobj','independent_delta_subcritical_Q0_direction.sobj'],'field':K,'ring':P,'p_images':dr,'RR_kernel':v,'e_numerator':en,'C_determinant':det,'items':items,'checks':checks,'scope':__doc__},str(tmp),compress=False);os.replace(tmp,p/'independent_delta_subcritical_H_initials.sobj')
(p/'independent_delta_subcritical_H_initials.json').write_text(json.dumps({'checks':len(checks),'rows':len(items),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
