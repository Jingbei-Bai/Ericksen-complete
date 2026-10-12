"""Separate base content and perform the exact finite Klein invariant map.

Unprotected base factors are recorded, never silently discarded in a
uniform-spectrum claim. Neither det nor a determinant norm is divided out.
"""
from pathlib import Path
import time,json,os
from sage.all import QQ,PolynomialRing,load,save
from cysignals.alarm import alarm,cancel_alarm
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'pencil_fast_numerators.sobj'))
CR=d['base_ring'];K=CR.fraction_field();a,q,T,Hb=CR.gens();r=q-1;m=1+a*r
oldguard=a*q*r*m*T*(T-1)*Hb*(Hb+T)*(Hb+1)*(1-Hb-T)
original_spectral_guard=a-1
audit={};unexpected=[]
for label,den in [('graph',d['common_graph_denominator'])]+[(f'source_{i}',p) for i,p in enumerate(d['source_denominators'])]:
    factors=[]
    for f,e in den.factor():
        protected=not oldguard.quo_rem(f)[1]
        also_original=not (oldguard*original_spectral_guard).quo_rem(f)[1]
        factors.append({'factor':str(f),'power':int(e),'listed_guard':protected,'original_guard':also_original})
        if not also_original:unexpected.append((label,str(f)))
    audit[label]=factors
assert not unexpected,unexpected
print('denominator audit PASS',round(time.monotonic()-start,2),flush=True)
grouped=[]
for p in d['rows']:
    temp={}
    for exp,c0 in p.dict().items():temp.setdefault(tuple(exp[:3]),{})[tuple(exp[3:])]=c0
    grouped.append({mon:CR(coeff) for mon,coeff in temp.items()})
    print('grouped',len(temp),round(time.monotonic()-start,2),flush=True)
R=PolynomialRing(K,names=('u','v','s'));RI=PolynomialRing(K,names=('U','V','S','W'));primitive=[];contents=[];invariants=[];records=[]
for i,coeff in enumerate(grouped):
    ordered=sorted(coeff,key=lambda mon:coeff[mon].number_of_terms());content=coeff[ordered[0]]
    alarm(90)
    try:
        for mon in ordered[1:]:
            content=content.gcd(coeff[mon])
            if content.is_constant():break
    finally:cancel_alarm()
    content=CR(content/content.lc())
    reduced={}
    for mon,c0 in coeff.items():
        quo,rem=c0.quo_rem(content);assert not rem;reduced[mon]=quo
    poly=R(reduced);primitive.append(poly);contents.append(content)
    ff=[(str(f),int(e),not oldguard.quo_rem(f)[1],not(oldguard*original_spectral_guard).quo_rem(f)[1]) for f,e in content.factor()]
    invdict={}
    for mon,c0 in reduced.items():
        ex=list(mon)
        if i==0:ex[2]+=1 # p23 != 0, so multiply the anti-invariant H by s.
        parity=ex[0]%2;assert all(e%2==parity for e in ex)
        key=tuple((e-parity)//2 for e in ex)+(parity,)
        invdict[key]=c0
    inv=RI(invdict);invariants.append(inv)
    # Direct embedding verifies the invariant representation without a norm.
    u,v,s=R.gens();mapped=inv(u*u,v*v,s*s,u*v*s)
    assert mapped==(s*poly if i==0 else poly)
    rec={'row':i,'fiber_terms':len(reduced),'joint_degree':int(poly.total_degree()),
         'base_terms_after_content':sum(c0.number_of_terms() for c0 in reduced.values()),
         'content_factors':[{'factor':f,'power':e,'listed_guard':ok,'original_guard':ok2} for f,e,ok,ok2 in ff],
         'invariant_terms':len(inv.dict()),'invariant_W_degree':int(inv.degree(RI.gen(3)))}
    records.append(rec);print(json.dumps(rec),'elapsed',round(time.monotonic()-start,2),flush=True)
    payload={'base_ring':CR,'fiber_ring':R,'invariant_ring':RI,'parameters':tuple(map(K,CR.gens())),
       'base_contents':contents,'primitive_rows':primitive,'invariant_rows':invariants,
       'denominator_audit':audit,'records':records,'old_guard':oldguard,
       'source_file':'pencil_fast_numerators.sobj','invariant_relation':'W^2=U*V*S',
       'scope':'Exact base-content split and Klein map. Any content factor outside original guards needs an all-spectrum/c generic-Hb audit. No dimension or coprimality assertion.'}
    tmp=out/'pencil_base_content_and_klein_writing.sobj';save(payload,str(tmp));os.replace(tmp,out/'pencil_base_content_and_klein.sobj')
    (out/'pencil_base_content_and_klein.json').write_text(json.dumps({'status':'PARTIAL' if i==0 else 'COMPLETE_EXACT_SPLIT','denominator_audit':audit,'records':records,'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']},indent=2)+'\n')
print('DONE',round(time.monotonic()-start,2),flush=True)
