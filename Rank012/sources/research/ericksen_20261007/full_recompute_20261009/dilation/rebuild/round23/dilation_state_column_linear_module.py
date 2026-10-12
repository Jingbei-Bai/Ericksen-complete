"""Full differential-ideal tangent module at a cylindrical state.

The first numeric run is diagnostic only.  No material specialization is a
positive theorem.  Closure under all three linearized E_i uses their exact
vanishing at the reference point, rather than a finite-order truncation.
"""
from sage.all import *
from pathlib import Path
import pickle,json,sys,time,signal

OUT=Path(__file__).absolute().parent
START=time.time()
signal.alarm(150)
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
MODE=sys.argv[1] if len(sys.argv)>1 else 'diagnostic'
if MODE=='diagnostic':
    K=QQ;a=K(1)/2;b=K(1)/3;p=K(3)/5;q=K(4)/5
elif MODE=='material':
    K=FractionField(PolynomialRing(QQ,names=('a','b')))
    a,b=K.gens();p=K(3)/5;q=K(4)/5
elif MODE in ('pzero','qzero'):
    K=FractionField(PolynomialRing(QQ,names=('a','b')))
    a,b=K.gens();p=K(MODE=='qzero');q=K(MODE=='pzero')
elif MODE=='resonance':
    K=FractionField(PolynomialRing(QQ,names=('a','z')))
    a,z=K.gens();p=2*z/(1+z*z);q=(1-z*z)/(1+z*z)
    b=a*p*p+q*q
else:
    K=FractionField(PolynomialRing(QQ,names=('a','b','z')))
    a,b,z=K.gens();p=2*z/(1+z*z);q=(1-z*z)/(1+z*z)

vs=list(map(K,[0,0,p,0,0,q,0,0,0,p*q,0,0,0,0,0,0,0]))
v=vector(K,[q,-p,0]);ax=vector(K,[0,0,1]);w=ax.cross_product(v)
weights=[1]*9+[2]*6+[3]*2
base=vector(K,vs+list(v)+list(ax))
def cf(c):
    # The stored coefficient field uses the same material names.
    return K(c(a=a,b=b)) if hasattr(c,'subs') else K(c)
def val(poly):
    return sum(cf(c)*prod(vs[i]**ex[i] for i in range(17)) for ex,c in poly.items())
def jac(poly):
    out=[K(0)]*23
    for ex,c in poly.items():
        cc=cf(c)
        for j in range(17):
            if ex[j]:
                out[j]+=cc*ex[j]*prod(vs[i]**(ex[i]-(i==j)) for i in range(17))
    return vector(K,out)

F0=[[val(z) for z in row] for row in src['Fraw']]
L=[]
checks={}
# Frame derivative is tiny and evaluated by a separate degree-three ring.
PR=PolynomialRing(K,names=tuple('P'+str(i) for i in range(9))+('v0','v1','v2','a0','a1','a2'))
gg=PR.gens();pp=matrix(PR,3,3,gg[:9]);vv=vector(PR,gg[9:12]);aa=vector(PR,gg[12:15]);ww=aa.cross_product(vv)
framepoint=list(vs[:9])+list(v)+list(ax)
for i in range(3):
    rows=[]
    for j in range(17):
        z=jac(src['Fraw'][i][j])
        z[j]+=weights[j]*v[i]
        z[17+i]+=weights[j]*vs[j]
        rows.append(z)
        checks['equilibrium_state_%s_%s'%(i,j)]=F0[i][j]+weights[j]*v[i]*vs[j]==0
    Wi=pp.row(i)-pp.trace()/2*vector(PR,[int(i==j) for j in range(3)])
    frame=list(ww[i]*ww-Wi.cross_product(vv))+list(-Wi.cross_product(aa))
    for j,z in enumerate(frame):
        checks['equilibrium_frame_%s_%s'%(i,j)]=z(*framepoint)==0
        grad=[K(z.derivative(g)(*framepoint)) for g in gg]
        rows.append(vector(K,grad[:9]+[K(0)]*8+grad[9:]))
    L.append(matrix(K,rows))
assert all(checks.values()),[n for n,x in checks.items() if not x]

rows=[]
for name,polys in [('C3',src['C3']),('Q4',src['Q4']),('G5',src['basis'][:27])]:
    for j,poly in enumerate(polys):
        assert val(poly)==0
        rows.append(jac(poly))
rows.extend([vector(K,[0]*17+list(2*v)+[0]*3),vector(K,[0]*20+list(2*ax)),vector(K,[0]*17+list(ax)+list(v))])
Jraw=matrix(K,rows)
J=Jraw.row_space().basis_matrix()
history=[J.nrows()]
print('initial rank',history[-1],flush=True)
for step in range(24):
    JJ=J.stack(J*L[0]).stack(J*L[1]).stack(J*L[2]).row_space().basis_matrix()
    history.append(JJ.nrows());print('closure',step,history[-1],flush=True)
    if JJ.nrows()==J.nrows(): break
    J=JJ
else: raise AssertionError('finite-dimensional closure did not stabilize')
N=J.right_kernel().basis_matrix().transpose()
E=[]
for i in range(3):
    Ei=N.solve_right(L[i]*N)
    assert L[i]*N==N*Ei
    E.append(Ei)
Dr=sum((v[i]*E[i] for i in range(3)),matrix.zero(K,N.ncols()))
Da=sum((w[i]*E[i] for i in range(3)),matrix.zero(K,N.ncols()))
Dz=sum((ax[i]*E[i] for i in range(3)),matrix.zero(K,N.ncols()))
checks['RA_commute']=Dr*Da-Da*Dr==0
checks['AZ_commute']=Da*Dz-Dz*Da==0
# [R,Z]=Z gives matrix convention DZ*DR-DR*DZ=DZ.
checks['RZ_bracket']=Dz*Dr-Dr*Dz==Dz
rec={'mode':MODE,'scope':'diagnostic' if MODE=='diagnostic' else 'generic parameter field only; pivot walls not covered',
     'closure_ranks':history,'tangent_dimension':N.ncols(),'R_charpoly':str(Dr.charpoly().factor()),
     'A_charpoly':str(Da.charpoly().factor()),'Z_charpoly':str(Dz.charpoly().factor()),
     'checks':checks,'seconds':time.time()-START}
print(json.dumps(rec,indent=2),flush=True)
save({'field':K,'base':base,'ambient_jacobians':L,'raw_source_jacobian':Jraw,'source_jacobian':J,'tangent_basis':N,'E':E,'R':Dr,'A':Da,'Z':Dz,'checks':checks},str(OUT/('dilation_state_column_linear_module_'+MODE+'.sobj')))
(OUT/('dilation_state_column_linear_module_'+MODE+'.json')).write_text(json.dumps(rec,indent=2))
