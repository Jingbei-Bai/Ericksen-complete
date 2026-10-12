"""Complete support-gauged filtration for every rational growth 1<gamma<2."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'))
sg=load(str(p/'independent_ratio_p_boundary_support_gauges.sobj'))
K=g['coefficient_field'];a,b,h,z=K.gens();R=PolynomialRing(K,names=tuple('nf_p%d'%i for i in range(6))+tuple('nf_R%d'%i for i in range(3)));vs=R.gens()
def cp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items() if c})
def part(f,pd,rd):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(map(int,e[:6]))==pd and sum(map(int,e[6:]))==rd})
hat=matrix(R,3,9,[cp(f) for f in g['first_gauged_fields'].list()]);base=matrix(R,3,2,[cp(f) for f in g['base_fields'].list()]);Q=cp(g['source_Q0']);sources=vector(R,[cp(f) for f in g['six_full_sources']]);units=g['source_unit_coefficients'];perm=g['source_permutation']
Dp=matrix(R,3,6,[part(hat[i,j],1,1) for i in range(3) for j in range(6)]);Db=matrix(R,3,2,[part(base[i,j],0,1) for i in range(3) for j in range(2)])
records=[]
for i in range(3):
    radial=[]
    for j in range(6):
        for k in range(j+1,6):
            v=Dp[i,j]*vs[k]-Dp[i,k]*vs[j]
            if v:radial.append((j,k,v))
    ck('nonzero_nonradial_Newton_p_'+str(i),bool(radial))
    rec={'direction':i,'p_entries':[(j,str(f)) for j,f in enumerate(Dp.row(i)) if f],'first_nonradial_wedge':(radial[0][0],radial[0][1],str(radial[0][2]))};records.append(rec)
    print('TRUE PRINCIPAL p',rec,flush=True)
    for j in range(6):
        for e in hat[i,j].dict():
            pd,rd=sum(map(int,e[:6])),sum(map(int,e[6:]))
            # Degree(hat p)-1 <= gamma throughout [1,2] iff both endpoints.
            ck('full_p_filtration_%d_%d_%s'%(i,j,e),pd+rd-1<=1 and pd+2*rd-1<=2)
    for j in range(2):
        for e in base[i,j].dict():
            pd,rd=sum(map(int,e[:6])),sum(map(int,e[6:]))
            ck('full_base_filtration_%d_%d_%s'%(i,j,e),pd+rd<=1 and pd+2*rd<=2)
out=[]
for it in sg['items']:
    ci=int(it['chart']);m=cp(it['denominator']);W=matrix(R,3,6,[cp(v) for v in it['left_inverse_numerator'].list()]);D=matrix(R,6,3,[cp(v) for v in sg['D'].list()]);ck('support_leftinverse_'+str(ci),W*D==m*identity_matrix(R,3));ck('support_m_degree_'+str(ci),all(sum(map(int,e[:6]))==3 and not any(e[6:]) for e in m.dict()))
    fields=m*hat;wr=W*sources
    for i in range(3):fields[i,6+i]-=units[i]*wr[perm[i]]
    principal=matrix(R,3,9);basepr=matrix(R,3,2)
    for i in range(3):
        for j in range(9):
            ck('full_support_source_%d_%d_%d'%(ci,i,j),fields[i,j]==cp(it['full_field_numerators'][i,j]))
            for e in fields[i,j].dict():
                pd,rd=sum(map(int,e[:6])),sum(map(int,e[6:]))
                # numerator/m contributes degree pd+gamma*rd-3;
                # target component weighs 1 for p, gamma for RR.
                off=1 if j<6 else 0;rg=0 if j<6 else 1
                ck('allgamma_degree_%d_%d_%d_%s'%(ci,i,j,e),pd+rd-3-off-rg<=1 and pd+2*rd-3-off-2*rg<=2)
            if j<6:
                principal[i,j]=part(fields[i,j],4,1);ck('p_principal_unchanged_%d_%d_%d'%(ci,i,j),principal[i,j]==m*Dp[i,j])
            else:principal[i,j]=part(fields[i,j],3,2)
        for j in range(2):basepr[i,j]=m*Db[i,j];ck('base_principal_unchanged_%d_%d_%d'%(ci,i,j),basepr[i,j]==part(cp(it['full_base_numerators'][i,j]),3,1))
    out.append({'chart':ci,'denominator':m,'principal_field_numerators':principal,'principal_base_numerators':basepr})
    print('CHART',ci,'ALL COMPONENTS / WHOLE GAMMA INTERVAL PASS',round(time.monotonic()-st,2),flush=True)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'p_principal_records':records,'scope':'All eight full necessary-source gauges have increment<=gamma for EVERY rational1<gamma<2, with correct principal pRR,pureRR^2,baseRR components. The p principal is nonzero and nonradial in every direction; delta(p)=0 is false. No detL,detC,RR or unselected p division. No deltaL/exclusion result asserted.'}
save({'checks':checks,'ring':R,'coefficient_field':K,'p_principal':Dp,'base_principal':Db,'items':out,'records':records},str(p/'independent_support_Newton_field_reader.sobj'),compress=False)
(p/'independent_support_Newton_field_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
