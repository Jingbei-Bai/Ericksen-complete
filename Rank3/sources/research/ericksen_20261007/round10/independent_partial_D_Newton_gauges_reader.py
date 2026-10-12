"""Each field needs only its own active p row: six exact partial-D gauges.

This proves regular Newton fields and horizontal rank on inactive rank2
RR-torus points. It is not a generic no-lift certificate for that layer.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));full=load(str(p/'independent_ratio_p_boundary_support_gauges.sobj'));ref=load(str(p/'independent_support_Newton_field_reader.sobj'));R=full['ring'];K=R.base_ring();v=R.gens();checks=[]
def ck(label,b):assert b,label;checks.append(label)
def cp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
hat=matrix(R,3,9,[cp(f) for f in d['first_gauged_fields'].list()]);base=matrix(R,3,2,[cp(f) for f in d['base_fields'].list()]);src=list(map(cp,d['six_full_sources']));units=d['source_unit_coefficients'];perm=d['source_permutation'];ck('source_permutation',tuple(perm)==(2,1,0))
choices=[((3,0,1),(4,1,1)),((2,2,1),(5,3,-1)),((0,4,1),(1,5,-1))];items=[]
for i,ls in enumerate(choices):
 for row,var,sign in ls:
  den=v[var];f=vector(R,[den*hat[i,j] for j in range(9)]);f[6+i]-=units[i]*sign*src[row];bb=vector(R,[den*base[i,j] for j in range(2)])
  ck('one_column_inverse',full['D'][row,perm[i]]==sign*den and all(full['D'][row,k]==0 for k in range(3) if k!=perm[i]))
  for j,ff in enumerate(f):
   for ex in ff.dict():
    pd,rd=sum(map(int,ex[:6])),sum(map(int,ex[6:]));off=1 if j<6 else 0;wg=0 if j<6 else 1
    ck('all_gamma_field_%d_%d_%d_%s'%(i,var,j,ex),pd+rd-1-off-wg<=1 and pd+2*rd-1-off-2*wg<=2)
  for ff in bb:
   for ex in ff.dict():
    pd,rd=sum(map(int,ex[:6])),sum(map(int,ex[6:]));ck('all_gamma_base',pd+rd-1<=1 and pd+2*rd-1<=2)
  for chart in full['items']:
   sel=chart['selection'][perm[i]]
   if tuple(sel)!=(row,var,sign):continue
   mult=chart['denominator']//den
   ck('exact_full_chart_overlap',vector(R,list(chart['full_field_numerators'].row(i)))==mult*f)
   ck('exact_base_chart_overlap',vector(R,list(chart['full_base_numerators'].row(i)))==mult*bb)
  items.append({'direction':i,'source_row':row,'p_index':var,'source_sign':sign,'denominator':den,'field_numerators':f,'base_numerators':bb})
# Exact finite support argument; only necessary binomial support equality
# is used, valid because every coefficient of L is a protected nonzero unit.
support=[]
for rm in range(1,8):
 rr=[bool(rm&(1<<i)) for i in range(3)]
 for pm in range(1,64):
  pp=[bool(pm&(1<<i)) for i in range(6)]
  if any([rr[0] and pp[0],rr[0] and pp[1],rr[1] and pp[2]][i] != [rr[2] and pp[5],rr[1] and pp[3],rr[2] and pp[4]][i] for i in range(3)):continue
  active=[pp[0] or pp[1],pp[2] or pp[3],pp[4] or pp[5]];horizontal=[i for i in range(3) if active[i] and rr[i]]
  if sum(rr)==3:ck('torus_two_available_fields_%s'%pm,len(horizontal)>=2)
  if sum(rr)==2 and len(horizontal)<2:
   missing=rr.index(False);ck('two_RR_exception_only_missing_row',sum(active)==1 and active[missing])
  support.append({'RR_mask':rm,'p_mask':pm,'D_rank':sum(active),'available_horizontal_directions':horizontal})
save({'source_files':['independent_ratio_p_boundary_second_gauge.sobj','independent_ratio_p_boundary_support_gauges.sobj'],'ring':R,'items':items,'support_records':support,'checks':checks,'scope':__doc__},str(p/'independent_partial_D_Newton_gauges_reader.sobj'),compress=False)
(p/'independent_partial_D_Newton_gauges_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'partial_gauges':len(items),'incidence_support_masks':len(support),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),'SUPPORTS',len(support),'seconds',time.monotonic()-st,flush=True)
