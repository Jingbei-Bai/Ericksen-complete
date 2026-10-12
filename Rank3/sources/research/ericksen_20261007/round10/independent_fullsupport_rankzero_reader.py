"""Complete five-wall actual-open exclusion of rank-zero amplitude matrices."""
from pathlib import Path
from itertools import product
import time,json
from sage.all import load,save,matrix,vector,prod
out=Path(__file__).resolve().parent;start=time.monotonic();checks={};cert={}
raw=load(str(out/'general_curve_two_pair_offdiagonal.sobj'))
forms=load(str(out/'general_curve_full_amplitude_forms.sobj'))
norm=load(str(out/'general_curve_homogeneous_normalized.sobj'))
lead=load(str(out/'general_curve_dimensionless_leading.sobj'))
pc=load(str(out/'independent_fullsupport_rankzero_probe.sobj'))
ac=load(str(out/'independent_fullsupport_rankzero_axis_probe.sobj'))
wc=load(str(out/'independent_fullsupport_rankzero_five_walls.sobj'))
axis=load(str(out/'independent_fullsupport_axis_reader.sobj'))
R=raw['ring'];K=R.base_ring();C=K.ring();H,Z,S,N=K.gens();xs=R.gens();im=raw['images']
L=H*S-Z;P=L+S+1;M=L+S-1;D=L-S+1
old=[C(v) for v in [H,Z,H-Z,Z-1,H+1,H-Z+1,S-1,L+1,L+S]]
def unit(v,guardlist=old):
    if not v:return False
    for p in [v.numerator(),v.denominator()]:
        for f,e in C(p).factor():
            f=C(f)
            if f.is_constant():continue
            if not any(all((K(f)/g).derivative(c)==0 for c in C.gens()) for g in guardlist):return False
    return True
AA=matrix(R,6,3)
for row,f in enumerate(raw['cubics']):
    odd=norm['odd_pairs'][row]
    for mon,co in f.dict().items():
        degrees=[mon[2*i]+mon[2*i+1] for i in range(3)];degrees[odd]-=1
        assert sorted(degrees)==[0,0,2]
        AA[row,degrees.index(2)]+=R({tuple(mon):co})
    checks['original_amplitude_source_'+str(row)]=sum(AA.row(row))==f
checks['same_amplitude_matrix']=AA==forms['A']
delta=[L*P,S*M,D]
for i,case in enumerate(pc['cases']):
    phi=R.hom([x if j//2==i else R.zero() for j,x in enumerate(xs)],R)
    pure=[phi(f) for f in raw['cubics'] if phi(f)];B=case['linear_matrix'];x,y=xs[2*i:2*i+2]
    checks['original_pure_matrix_'+str(i)]=pure==[x*y*(B[j,0]*x+B[j,1]*y) for j in range(2)]
    checks['pure_delta_factor_'+str(i)]=B.det()==case['determinant'] and unit(B.det()/delta[i])
# Explicit ideal identities prove the three delta loci pairwise disjoint
# on the guards. No generic resultant or rank sample is used.
syzygies=[(L,S,-1,H,Z),(L,M,-1,1,S-1),(P,S,-1,H+1,Z-1),(P,M,1,-1,K(2)),
          (L,D,1,-1,S-1),(P,D,(H-1)/2,-(H+1)/2,Z-1),
          (S,D,H-1,-1,Z-1),(M,D,K(1)/2,-K(1)/2,S-1)]
for i,(f,g,a,b,target) in enumerate(syzygies):
    checks['delta_pair_identity_'+str(i)]=a*f+b*g==target and unit(target)
cert['delta_pair_syzygies']=syzygies
checks['all_axis_source_certificate']=len(axis['checks'])==169 and all(axis['checks'].values())
remaining=[];chosen=[]
for case in ac['cases']:
    free=case['free'];axes=tuple(case['axes']);tag=str(free)+''.join(map(str,axes));phi=R.hom(case['images'],R)
    if not case['hits']:
        remaining.append((free,axes));continue
    i,j,saved=case['hits'][0];f=phi(AA[i,j]);mon,co=next(iter(f.dict().items()))
    checks['axis_chart_source_'+tag]=f==saved and len(f.dict())==1 and unit(co)
    checks['axis_chart_only_nonzero_coordinates_'+tag]=all(not mon[k] for k in range(6) if k//2!=free)
    chosen.append({'free':free,'axes':axes,'row':i,'column':j,'source':f})
checks['all_nine_axis_charts']=len(chosen)==9 and remaining==[(0,(1,1)),(1,(1,0)),(2,(0,0))]
cert['nine_axis_chart_sources']=chosen

def cc(c):return K(c(0,0,H,Z,S,N))
mixed=[sum(cc(B[i,j])*im[2*u+i]*im[2*v+j] for i in range(2) for j in range(2)) for (u,v),B in lead['blocks'].items()]
trilinear=[sum(cc(co)*im[i]*im[2+j]*im[4+k] for co,(i,j,k) in zip(row,product(range(2),repeat=3))) for row in lead['cubic_rows']]
wall_polys=[L,P,S,M,D];derivative_units=[H,H+1,K.one(),H+1,H-1]
for i,(wall,derivunit) in enumerate(zip(wall_polys,derivative_units)):
    derivative=wall.derivative(H)+S*wall.derivative(Z)+N*wall.derivative(S)
    checks['actual_wall_contact_'+str(i)]=derivative==derivunit*N
checks['Hminusone_pivot_covered']=D==(H-1)*S-(Z-1)
newcases=[]
for i,case in enumerate(wc['cases']):
    free=case['free'];args=case['parameter_images'];label=case['wall']
    def spec(v):return K(v.numerator()(*args))/K(v.denominator()(*args))
    B=matrix(K,2,2,[spec(v) for v in pc['cases'][free]['linear_matrix'].list()])
    pulled=[]
    for g in old:
        value=spec(K(g));pulled.extend(C(f) for pp in [value.numerator(),value.denominator()] for f,e in C(pp).factor() if not C(f).is_constant())
    ratio=-B[0,0]/B[0,1]
    point=[K(v) for v in ([1,ratio,0,1,0,1] if free==0 else [0,1,1,ratio,1,0] if free==1 else [1,0,1,0,1,ratio])]
    checks[label+'_protected_linear_graph']=unit(B[0,1],pulled) and unit(ratio,pulled) and not B*vector(K,[1,ratio])
    checks[label+'_same_graph']=B==case['pure_matrix'] and ratio==case['ratio'] and point==case['direction']
    A=matrix(K,6,3,[spec(K(f(*point))) for f in AA.list()])
    vals=[spec(K(f(*point))) for f in mixed+trilinear]
    checks[label+'_original_all_source_values']=A==case['A'] and vals==case['direction_values']
    if label=='L':
        checks[label+'_contradiction']=A[2,0]==2*(H-Z) and unit(A[2,0],pulled)
    elif label=='S':
        checks[label+'_contradiction']=A[0,1]==2*(Z-1) and unit(A[0,1],pulled)
    elif label=='delta2':
        checks[label+'_contradiction']=vals[0]==-4*H**2*((H-1)*S+1)**2/(S*(S-1)) and unit(vals[0],pulled)
    elif label=='Eplus':
        checks[label+'_mixed_square']=vals[0]==4*H*(H+Z)**2/(Z-1)
        checks[label+'_unit_square_coefficient']=unit(vals[0]/(H+Z)**2,pulled)
        checks[label+'_forced_intersection']=K(args[2]).subs({Z:-H})==-1
        checks[label+'_final_contradiction']=A[2,0].subs({Z:-H})==4*H
    elif label=='Eminus':
        checks[label+'_mixed_square']=vals[0]==-4*H**4*(Z+1)**2/(H-Z)
        checks[label+'_unit_square_coefficient']=unit(vals[0]/(Z+1)**2,pulled)
        checks[label+'_forced_intersection']=K(args[2]).subs({Z:-1})==0
        checks[label+'_final_contradiction']=A[0,1].subs({Z:-1})==-4
    newcases.append({'wall':label,'pure_linear_matrix':B,'opposite_ratio':ratio,'A':A,'old_direction_sources':vals,'pulled_guards':pulled})
checks['all_five_walls']=[v['wall'] for v in newcases]==['L','Eplus','S','Eminus','delta2']
assert all(checks.values()),[k for k,v in checks.items() if not v]
cert.update({'checks':checks,'five_wall_cases':newcases,'scope':'No actual H interval can have any full-support necessary direction with amplitude matrix rank zero. Uses only old direction sources and six cubic amplitude entries, never a quartic amplitude condition. Includes every pure-column determinant wall and all two-axis charts.'})
save(cert,str(out/'independent_fullsupport_rankzero_reader.sobj'))
report={'status':'EXACT_COMPLETE_FULLSUPPORT_RANKZERO_OPEN_EXCLUSION_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_fullsupport_rankzero_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
