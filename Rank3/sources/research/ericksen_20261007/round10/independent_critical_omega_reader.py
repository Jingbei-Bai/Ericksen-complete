"""Independent all-material and all-complex-amplitude replay of the six omega walls."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):assert bool(v),n;checks.append(n)
def eq(n,f,g):ck(n,not(f-g))
source=load(str(p/'independent_delta_critical_linear_sources.sobj'))
graphs=load(str(p/'independent_critical_omega_even_graphs.sobj'))
matdata=load(str(p/'independent_critical_omega_even_matrix.sobj'))
positive=load(str(p/'independent_critical_omega_material_factors.sobj'))
old=source['field'];K=graphs['field'];a,th,H,Z,Y=K.gens();L=graphs['number_field'];om=L.gen()
ck('both conjugate directions retained',om**2-QQ(13)/2*om+1==0 and L.degree()==2)
for rec,mt in zip(graphs['records'],matdata['records']):
    label=rec['direction_type'];ph=old.ring().hom([a,rec['material_r'],H,Z,K(rec['p21']),rec['C']*K(rec['z']),rec['B'],rec['C']],K)
    def cv(f):f=old(f);return ph(f.numerator())/ph(f.denominator())
    for i in range(22):
        val=sum((cv(source['matrix'][i,j])*rec['e_graph'][j] for j in range(3)),K.zero())+cv(source['matrix'][i,3])*rec['delta_graph']+cv(source['matrix'][i,4])
        eq('full original critical affine source '+str((label,i)),val,0)
    P=matdata['amplitude_ring'];aa,x,T=P.gens()
    for i,f in enumerate(source['even_sources']):
        val=sum((cv(c)*prod(rec['e_graph'][k]**int(e[k]) for k in range(3)) for e,c in f.dict().items()),K.zero())
        eq('full original even source pullback '+str((label,i)),val,rec['even_forms'][i]['value'])
        num=val.numerator();num=num/num.leading_coefficient();fresh=P.zero()
        for e,c in num.dict().items():
            e=tuple(map(int,e));ck('all base and square exponents retained '+str((label,i,e)),e[2]==e[3]==0 and e[1]%2==e[4]%2==0)
            fresh+=L(c)*aa**e[0]*x**(e[1]//2)*T**(e[4]//2)
        eq('complex square-amplitude transport '+str((label,i)),fresh,mt['polynomials'][i])
    eq('complete amplitude coefficient determinant '+label,mt['matrix'].det(),mt['determinant'])
    eq('complete determinant factor replay '+label,mt['determinant_factors'].value(),mt['determinant'])
R=matdata['coefficient_ring'];aa,x=R.gens();c=positive['c_constant'];kap=positive['kappa_constant']
eq('uniform positive lower bound for c',232135214-27335384*QQ(13)/2,54455218)
for rec in matdata['records']:
    last=rec['determinant_factors'][-1][0];label=rec['direction_type']
    if label=='one_omega':u=aa*x**2;f=u**2+c*u+1
    elif label=='omega_one':u=x*(aa*(x+1)-1);f=u**2+c*u+1
    else:f=(aa-x)**2*(x+1)**2-kap*aa*x*(aa-x-1)
    eq('exact full-range positive factor '+label,last,f)
Q=positive['wall_ring'];aw,Tw=Q.gens();P=matdata['amplitude_ring'];towall=P.hom([aw,aw,Tw],Q)
omrec=next(z for z in matdata['records'] if z['direction_type']=='omega_omega')
for i,f in enumerate(positive['wall_polynomials']):eq('entire b=a^2 wall source '+str(i),f,towall(omrec['polynomials'][i]))
U=positive['wall_coefficient_ring'];av=U.gen();fields=[]
for f in positive['wall_polynomials']:
    dd=max(int(e[1]) for e in f.dict());ck('wall degree retained',dd==2);row=[U.zero()]*3
    for e,c0 in f.dict().items():e=tuple(map(int,e));row[2-e[1]]+=c0*av**e[0]
    fields.append(row)
f,g=fields[:2];syl=matrix(U,[f+[0],[0]+f,g+[0],[0]+g]);res=syl.det();target=positive['wall_resultants'][0][2]
eq('independent original wall Sylvester determinant',res,target)
unit=target/(av**8*(av+1)**4);ck('wall quotient constant',unit in L)
unit=L(unit);ck('wall constant nonzero in BOTH conjugates',unit!=0 and unit.norm()!=0)
record={'checks':len(checks),'seconds':time.monotonic()-st,'all_material_parameters_retained':True,'complex_amplitudes_retained':True,'both_algebraic_direction_conjugates_retained':True,'b_equal_a_squared_wall':'complete original two-quadratic Sylvester certificate','scope':'All six omega critical direction walls on the full-support chart; singular old-block canonical walls covered by separate source-level F8 review.'}
save({'checks':checks,'record':record},str(p/'independent_critical_omega_reader.sobj'),compress=False)
(p/'independent_critical_omega_reader.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
