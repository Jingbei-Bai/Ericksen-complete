from sage.all import load,save,PolynomialRing,QQ,matrix
from pathlib import Path
import json,time,os
out=Path(__file__).resolve().parent;start=time.monotonic();old=load(str(out/'independent_ratio_RR_normal_old_sources.sobj'));new=load(str(out/'independent_ratio_RR_normal_all_derivatives.sobj'));oldM=old['normal_matrix'];oldK=oldM.base_ring().base_ring();oldB=oldK.ring();B=PolynomialRing(QQ,names=('aa','rho'));a,r=B.gens();K=B.fraction_field();q=1+r;mp=oldB.hom([K(a),K(a*q/(1+a*r)),K(1/r),K(-q/r)],K);T=PolynomialRing(K,names=('R1','R2','R3'));rr=T.gens();guard=B(a*r*q*(a-1)*(1+a*r));cache={}
def cv(c):
 c=oldK(c)
 if c not in cache:cache[c]=mp(c.numerator())/mp(c.denominator())
 return cache[c]
def protected(c):return all(not guard.quo_rem(f)[1] for v in (c.numerator(),c.denominator()) for f,e in v.factor())
rawrows=[];names=[]
for i in range(7):rawrows.append([{tuple(e):c for e,c in oldM[i,j].dict().items()} for j in range(7)]);names.append(old['source_names'][i])
for it in new['items']:rawrows.append(it['columns']);names.append('V'+str(it['direction'])+'_'+it['source_name'])
rows=[];scales=[];records=[]
for i,raw in enumerate(rawrows):
 row=[T({tuple(e):cv(c) for e,c in col.items()}) for col in raw];poss=[c for f in row for c in f.coefficients() if c and protected(c)];scale=poss[0] if poss else K.one();row=[T({tuple(e):c/scale for e,c in f.dict().items()}) for f in row];rows.append(row);scales.append(scale);rec={'row':i,'source':names[i],'terms':sum(len(f.dict()) for f in row),'parameter_scale_protected':True,'coefficient_a_derivatives_zero':all(c.derivative(a)==0 for f in row for c in f.coefficients()),'seconds':time.monotonic()-start};records.append(rec);print(rec,flush=True)
 payload={'parameter_ring':B,'coefficient_field':K,'RR_ring':T,'parameters':(a,r),'RR_variables':rr,'old_parameter_map':list(mp.im_gens()),'matrix':matrix(T,rows),'source_names':names[:len(rows)],'row_scales':scales,'records':records,'guard':guard,'source_files':['independent_ratio_RR_normal_old_sources.sobj','independent_ratio_RR_normal_all_derivatives.sobj'],'scope':'Full a,r spectral family at H=1,Z=-1, exact source row specialization and protected parameter scaling. No rank assertion.'};tmp=out/'independent_ratio_RR_cyclic_full_matrix.tmp.sobj';save(payload,str(tmp),compress=False);os.replace(tmp,out/'independent_ratio_RR_cyclic_full_matrix.sobj');jt=out/'independent_ratio_RR_cyclic_full_matrix.json.tmp';jt.write_text(json.dumps({'status':'BUILDING' if len(rows)<28 else 'EXACT_ALL_SPECTRA_BASE_SECTION_MATRIX','records':records},indent=2)+'\n');os.replace(jt,out/'independent_ratio_RR_cyclic_full_matrix.json')
print('DONE',time.monotonic()-start,flush=True)

