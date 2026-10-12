"""Independent all-parameter source replay and conic finiteness audit.

This reader constructs the graph directly, rather than importing its creator's
directions, conics, determinants, or checking routines.
"""
from pathlib import Path
import json, time
from sage.all import *

out = Path(__file__).resolve().parent
started = time.monotonic()
full = load(str(out / 'general_curve_dimensionless_full.sobj'))
norm = load(str(out / 'general_curve_homogeneous_normalized.sobj'))
K = full['ring'].base_ring()
a,r,H,Z,S,N = K.gens()
L = H*S-Z
point = [S/(2*Z), K.one(), -L/(2*H*Z), K.one(), -L/(2*H*Z), K.one()]
nv = S*L/(H*Z)
checks = {}
def check(name, value):
    checks[name] = bool(value)
    assert checks[name], name
def subst_n(c):
    return K(c.numerator()(a,r,H,Z,S,nv))/K(c.denominator()(a,r,H,Z,S,nv))
for i, f in enumerate(full['rows']):
    degree = 2 if i < 3 else 3
    v = sum(subst_n(c)*prod(point[j]**int(m[j]) for j in range(6))
            for m,c in f.dict().items() if sum(m)==degree)
    check('original_full_leading_graph_'+str(i), v==0)

R = PolynomialRing(K, names=('X','Y','T'))
X,Y,T = R.gens()
qactual=[]
for f in norm['quartics']:
    v=R.zero()
    for m,c in f.dict().items():
        powers=[sum(m[2*j:2*j+2]) for j in range(3)]
        check('even_pair_degrees_'+str(len(checks)), all(d%2==0 for d in powers))
        cc=K(c.numerator()(H,Z,S,nv))/K(c.denominator()(H,Z,S,nv))
        v+=cc*prod(point[j]**int(m[j]) for j in range(6))*prod(R.gen(j)**(int(powers[j])//2) for j in range(3))
    qactual.append(v)

C=PolynomialRing(QQ,'k'); k=C.gen()
P=PolynomialRing(C,names=('x','y','z'));x,y,z=P.gens()
qq=[(k*k+1)*x*x+(5*k-1)*x*y+k*(1-k)*y*y+(k*k-5*k)*x*z-(k-1)**2*y*z+(k-1)*z*z,
    k*(k-1)*x*x+(5*k-4)*x*y+(-k*k+2*k-2)*y*y+k*k*x*z+(-k*k-3*k+4)*y*z+k*z*z,
    (k-1)*x*x+x*y-k*y*y+(-4*k*k+5*k)*x*z+(-4*k*k+3*k+1)*y*z+(2*k*k-2*k+1)*z*z]
scales=[S/(8*H**3*Z),L/(8*H**4*Z),-S*L/(8*H**3*Z**2)]
for i,f in enumerate(qq):
    mapped=R({tuple(m):K(c(H*S/Z)) for m,c in f.dict().items()})
    check('all_parameter_quartic_identity_'+str(i),qactual[i]==scales[i]*mapped)
matrices=[matrix(C,3,3,lambda i,j:f.monomial_coefficient(P.gen(i)**2) if i==j else f.monomial_coefficient(P.gen(i)*P.gen(j))/2) for f in qq]
determinants=[-4*(k-2)*(k-QQ(1)/2)*k*(k-1)**2,
              -4*(k-1)*(k-QQ(1)/2)*(k+1)*k*k,
              2*(k-2)*(k-1)*k*(k+1)]
for i in range(3):check('conic_determinant_'+str(i),matrices[i].det()==determinants[i])
check('nonproportional_q0_q1',qq[0].monomial_coefficient(y*y)*qq[1].monomial_coefficient(z*z)-qq[0].monomial_coefficient(z*z)*qq[1].monomial_coefficient(y*y)==-2*(k-1)**2)
check('q1_nondegenerate_k2',matrices[1].det()(2)!=0)
check('q2_nondegenerate_khalf',matrices[2].det()(QQ(1)/2)!=0)
check('q0_nonzero_k2',qq[0].monomial_coefficient(x*x)(2)!=0)
check('q0_nonzero_khalf',qq[0].monomial_coefficient(x*x)(QQ(1)/2)!=0)

# Reconstruct six amplitude linear rows from independently normalized cubics.
odd=[2,1,2,0,1,0]
AA=[]
for i,f in enumerate(norm['cubics']):
    row=[K.zero()]*3
    for m,c in f.dict().items():
        d=[sum(m[2*j:2*j+2]) for j in range(3)];d[odd[i]]-=1
        check('cubic_amplitude_parity_'+str(len(checks)),sorted(d)==[0,0,2])
        j=d.index(2)
        cc=K(c.numerator()(H,Z,S,nv))/K(c.denominator()(H,Z,S,nv))
        row[j]+=cc*prod(point[jj]**int(m[jj]) for jj in range(6))
    AA.append(row)
AM=matrix(K,AA)
walls=[]
for label,sv in [('S_zero',K.zero()),('L_zero',Z/H)]:
    def ev(c):return K(c.numerator()(a,r,H,Z,sv,0))/K(c.denominator()(a,r,H,Z,sv,0))
    B=matrix(K,6,3,lambda i,j:ev(AM[i,j]))
    qs=[R({tuple(m):ev(c) for m,c in f.dict().items()}) for f in qactual]
    expected=matrix(K,6,3)
    if label=='S_zero':
        expected[0,1]=(Z-1)/(4*H**3)
        expected[5,1]=-(Z-1)/(2*H**3)
    else:
        expected[2,0]=(H-Z)/(4*H**3)
        expected[4,0]=(H-Z)/(2*H**3)
    for ii in range(6):
        for jj in range(3):check(label+'_amplitude_entry_'+str(ii)+'_'+str(jj),B[ii,jj]==expected[ii,jj])
    record={'wall':label,'A':B,'quartics':qs}
    walls.append(record)
    print(label,'A',[[str(c.factor()) if c else '0' for c in row] for row in B.rows()],flush=True)
    print(label,'Q',[str(f.factor()) if f else '0' for f in qs],flush=True)

save({'checks':checks,'amplitude_matrix':AM,'quartic_amplitudes':qactual,'exceptional_walls':walls},str(out/'independent_simple_direction_reader.sobj'))
report={'all_pass':all(checks.values()),'checks_count':len(checks),'checks':checks,'seconds':time.monotonic()-started,'scope':'Independent full source and all-parameter conic identities; k=0,1 excluded for full support by exact amplitude entries using only H,Z-1,H-Z old guards. Entire explicit direction graph amplitude is finite; direction completeness not claimed.'}
(out/'independent_simple_direction_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
