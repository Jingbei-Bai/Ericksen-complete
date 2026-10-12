"""Global full-source pair syzygies, ordinary initials and small-span screen."""
from pathlib import Path
from sage.all import load,save,PolynomialRing,matrix,vector
import time,json
out=Path(__file__).resolve().parent;start=time.monotonic()
fd=load(str(out/'independent_ratio_weighted_degreefall.sobj'))
oldpair=load(str(out.parent/'round9/joint_fiber_pair_combinations.sobj'))
orig=load(str(out.parent/'round8/rank3_joint_fiber_leading.sobj'))
K=fd['ring'].base_ring();R=PolynomialRing(K,names=tuple('gp_p%d'%i for i in range(6))+tuple('gp_rr%d'%i for i in range(3)));v=R.gens()
def cp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items() if c})
def top(f,n):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(map(int,e))==n})
src={int(it['source_index']):cp(it['full']) for it in fd['items']};q0=cp(fd['Q0']);spec=[((10,5,1),(11,4,1)),((13,3,1),(19,2,1)),((15,1,1),(18,0,-1))]
pairvars=((0,2),(1,4),(3,5));items=[];records=[];checks=0
colors=((1,0),(0,1),(1,0),(1,1),(0,1),(1,1),(1,1),(0,1),(1,0))
def char(e):return tuple(sum(int(e[j])*colors[j][k] for j in range(9))%2 for k in range(2))
for poly in [q0]+[cp(orig['quartic_source_rows'][j]) for j in (1,5,6)]:
 assert all(char(e)==(0,0) for e in poly.dict());checks+=1
for i,parts in enumerate(spec):
 f=sum((sgn*v[j]*src[idx] for idx,j,sgn in parts),R.zero());assert f==cp(oldpair['polynomials'][i]);checks+=1
 hist={}
 for e in f.dict():hist[sum(map(int,e))]=hist.get(sum(map(int,e)),0)+1
 assert max(hist)==4;checks+=1
 assert all(char(e)==colors[pairvars[i][0]] for e in f.dict());checks+=1
 comparisons=[]
 for label,target,basis in [('full',f,[v[j]*q0 for j in pairvars[i]]),('ordinary_top',top(f,4),[v[j]*top(q0,3) for j in pairvars[i]])]:
  dd=[{tuple(map(int,e)):c for e,c in g.dict().items()} for g in basis];td={tuple(map(int,e)):c for e,c in target.dict().items()};support=set().union(*[set(z) for z in dd]);outside=sorted(set(td)-support)
  if outside:
   witness=outside[0];record={'scope':label,'member':False,'reason':'target monomial absent from both parity-compatible p*Q0 sources','witness_exponent':witness,'witness_coefficient':td[witness]};checks+=1
  else:
   exps=sorted(support);A=matrix(K,[[z.get(e,K.zero()) for z in dd] for e in exps]);ids=list(A.transpose().pivots());assert len(ids)==2;checks+=1
   y=vector(K,[td.get(exps[j],K.zero()) for j in ids]);co=A.matrix_from_rows(ids).inverse()*y;rem=target-sum((co[j]*basis[j] for j in range(2)),R.zero());record={'scope':label,'member':not bool(rem),'coefficients':co,'remainder':rem}
  comparisons.append(record)
 rec={'pair':i,'terms':len(f.dict()),'ordinary_layers':hist,'RR_top_degree':max(sum(map(int,e[6:])) for e in top(f,4).dict()),'small_span_membership':[bool(z['member']) for z in comparisons],'seconds':time.monotonic()-start};records.append(rec)
 items.append({'pair':i,'full':f,'ordinary_initial':top(f,4),'full_source_combination':parts,'comparison_candidates':pairvars[i],'small_span_comparisons':comparisons,'record':rec});print(rec,flush=True)
save({'ring':R,'coefficient_field':K,'items':items,'Q0':q0,'checks':checks,'scope':'Exact global full identities, no p/state division, including inactive D and Delta=0. The p*Q0 screen is over the generic coefficient field and is not a full ideal test. Original Q1,Q5,Q6 have trivial Klein character and cannot supply a constant-coefficient contribution in these three nontrivial parity sectors.'},str(out/'independent_ratio_p_boundary_global_pair_syzygies.sobj'),compress=False)
(out/'independent_ratio_p_boundary_global_pair_syzygies.json').write_text(json.dumps({'checks':checks,'records':records,'seconds':time.monotonic()-start},indent=2)+'\n');print('DONE',checks,time.monotonic()-start,flush=True)
