"""Independent full-spectrum audit from the original 17-state vector fields.

Derive the four jet reconstructions anew before reading root certificates.
No fixed-spectrum substitution and no numerical solution are used.
"""
from pathlib import Path
from functools import lru_cache
import pickle,json,time,hashlib,itertools
import sympy as sy
from sage.all import QQ,PolynomialRing,matrix,vector,load,save

OUT=Path(__file__).resolve().parent; BASE=OUT.parent
start=time.time(); path=BASE/'generic_frame_complete_api.pkl'
source=pickle.load(path.open('rb'))
variables=source['variables']; parameters=source['parameters']; F=source['vectorfields']
P=PolynomialRing(QQ,names=('a','b','h','z')); a,b,h,z=P.gens(); K=P.fraction_field()
R=PolynomialRing(K,names=('y','w','t0','t1','t2','t3','t4','t5'))
y,w,t0,t1,t2,t3,t4,t5=R.gens()
A=a*(b-1); C=b*(a-1); D=a-b
checks={}
def check(label,test):
    assert test,label
    checks[label]=True

# This checks the normalization law from the unnormalized source itself.
weights=[1]*9+[2]*6+[3]*2
for i in range(3):
    for j in range(9):
        polynomial=sy.Poly(F[i][j],*variables)
        check('weight2_F%d_pi%d'%(i,j),all(sum(e*k for e,k in zip(mon,weights))==2
              for mon,coefficient in polynomial.terms() if coefficient!=0))

pi=[R(1),y,R(0),R(D*h/A)*y,R(h),w,R(0),R(-A*z/(C*h))*w,R(z)]
context=dict(zip(variables[:15],pi+[t0,t1,t2,t3,t4,t5]))
context.update({parameters[0]:R(a),parameters[1]:R(b)})
@lru_cache(None)
def evaluate(expression):
    if expression in context:return context[expression]
    if expression.is_Rational:return R(QQ(int(expression.p))/int(expression.q))
    if expression.is_Add:return sum((evaluate(x) for x in expression.args),R.zero())
    if expression.is_Mul:
        result=R.one()
        for x in expression.args:result*=evaluate(x)
        return result
    if expression.is_Pow:
        value=evaluate(expression.base); exponent=int(expression.exp)
        if exponent<0:
            assert value.is_constant() and value!=0
            return R(K(value.constant_coefficient())**exponent)
        return value**exponent
    raise ValueError(expression)

# Verify the three Bcal equations, without imposing any jet reconstruction.
gamma=[C/A,D/C,-A/D]
Bcal=[pi[8]*pi[5]+R(gamma[0])*pi[4]*pi[7],
      pi[0]*pi[6]+R(gamma[1])*pi[8]*pi[2],
      pi[4]*pi[1]+R(gamma[2])*pi[0]*pi[3]]
check('Bcal_all_zero',all(v==0 for v in Bcal))
beta=[K(1),K(a),K(b)]; d_original=[K(D),K(A),K(-C)]
for permutation in itertools.permutations(range(3)):
    inversions=sum(permutation[i]>permutation[j] for i in range(3) for j in range(i+1,3))
    sign=(-1)**inversions
    rebased=[beta[j]/beta[permutation[0]] for j in permutation]
    new_d=[rebased[i]*(rebased[(i+1)%3]-rebased[(i+2)%3]) for i in range(3)]
    common=sign/beta[permutation[0]]**2
    check('spectral_relabeling_'+'_'.join(map(str,permutation)),
          new_d==[common*d_original[j] for j in permutation])

# All six derivatives of the zero pair are retained. First solve four.
tangent={(i,j):evaluate(F[i][j]) for i in range(3) for j in (2,6)}
four=[tangent[i,j] for i in (0,2) for j in (2,6)]
eliminated=[t0,t1,t2,t4]
M=matrix(K,4,4,[K(q.monomial_coefficient(t)) for q in four for t in eliminated])
remainder=vector(R,[q-sum(R(M[i,j])*eliminated[j] for j in range(4))
                    for i,q in enumerate(four)])
check('four_tangencies_linear_in_reconstructed_jets',
      all(not any(q.degree(t)>0 for t in eliminated) for q in remainder))
detM=M.det(); check('reconstruction_determinant_nonzero',detM!=0)
solution=-(M.inverse().change_ring(R)*remainder)
expected=[-pi[3]*z,R(0),2*pi[3]*h-R(2*A/C)*pi[3]*z,R(0)]
check('four_reconstructions_match',list(solution)==expected)
jet_map=dict(zip(eliminated,solution))
sub=R.hom([jet_map.get(v,v) for v in R.gens()],R)
check('four_coordinate_tangencies_vanish',all(sub(q)==0 for q in four))

qY=sub(tangent[1,2]); qW=sub(tangent[1,6])
cy=qY.monomial_coefficient(y*y); cw=qW.monomial_coefficient(w*w)
check('y_square_equation_affine',qY==cy*y*y+qY.constant_coefficient())
check('w_square_equation_affine',qW==cw*w*w+qW.constant_coefficient())
check('square_equation_coefficients_nonzero',cy!=0 and cw!=0)
Y=-qY.constant_coefficient()/cy; W=-qW.constant_coefficient()/cw

# For V_i=E_i/p11 and normalized pi_j=original_pi_j/p11,
# V_i(pi_j)=F_i(pi_j)-pi_j F_i(p11), using weight two above.
V=[]
for i in (0,2):
    sigma=sub(evaluate(F[i][0]))
    V.append([sub(evaluate(F[i][j]))-pi[j]*sigma for j in (4,8,1,5)])
v1h,v1z,v1y,v1w=V[0]; v3h,v3z,v3y,v3w=V[1]
mu2=1-D*z/C
check('normalized_V1_h',v1h==R(-2*A*mu2/D)*w)
check('normalized_V1_z',v1z==R(-2*C*mu2/D)*w)
check('normalized_V3_h',v3h==R(2*h*mu2*(A+D*h)/A)*y)
check('normalized_V3_z',v3z==R(2*D*h*z*mu2/A)*y)
eY=2*y*v3y-R(Y.derivative(h))*v3h-R(Y.derivative(z))*v3z
eW=2*w*v1w-R(W.derivative(h))*v1h-R(W.derivative(z))*v1z
check('first_tangency_divisible_by_y',all(m[0]>0 for m in eY.dict()))
check('second_tangency_divisible_by_w',all(m[1]>0 for m in eW.dict()))
eY=R(eY/y); eW=R(eW/w)
def square_reduce(q):
    result=K.zero()
    for exponent,coefficient in q.dict().items():
        assert exponent[0]%2==exponent[1]%2==0 and not any(exponent[2:])
        result+=coefficient*Y**(exponent[0]//2)*W**(exponent[1]//2)
    return result
f,g=square_reduce(eY),square_reduce(eW)
NY=a*b*(a*b-1)-2*a*b*(b-1)*h-b*D*h**2+2*a*b*(a-1)*z-3*a*D*z**2
N3=a*b*(a*b-1)+2*a*b*(b-1)*h+3*b*D*h**2-2*a*b*(a-1)*z+a*D*z**2
check('compact_Y_formula',Y==-NY/(4*h*b*D))
check('compact_V3_y_formula',v3y==R(1-2*D**2*h*z/(A*C))*y*y+R(N3/(4*b*A)))
compact_f=2*(1-2*D**2*h*z/(A*C))*Y+N3/(2*b*A)-2*h*mu2/A*((A+D*h)*Y.derivative(h)+D*z*Y.derivative(z))
check('compact_necessary_f_formula',f==compact_f)
num=f.numerator(); denominator=f.denominator()
coefficient=sum((c*a**m[0]*b**m[1] for m,c in num.dict().items()
                 if tuple(m)[2:]==(3,1)),P.zero())
check('necessary_polynomial_nonzero_each_spectrum',coefficient==4*b*(a-b)**3)
L=C*(A+D*h)-A*D*z
det_velocity=v1h*v3z-v1z*v3h
check('ratio_velocity_determinant',det_velocity==R(4*h*mu2**2*L/(A*D))*y*w)
zline=K(C*(A+D*h)/(A*D))
restrict=P.hom([K(a),K(b),K(h),zline],K)
restricted=restrict(num)
constant=restricted.numerator()(a=a,b=b,h=0,z=0)/restricted.denominator()(a=a,b=b,h=0,z=0)
constant=K(constant)
constant_expected=K(4*a**3*(b-1)**3*(a-1)*b**2/(b-a))
check('restricted_constant_nonzero_each_spectrum',constant==constant_expected)
check('restricted_polynomial_in_h',restricted.denominator().degree(h)==0)

# Compare only after independently deriving the quantities above.
root_open=load(str(OUT/'zero_pair_open_branch.sobj'))
root_rank=load(str(OUT/'zero_pair_ratio_rank_obstruction.sobj'))
check('root_Y_W_match',Y==K(root_open['Y']) and W==K(root_open['W']))
check('root_necessary_f_g_match',f==K(root_open['ratio_equations'][0]) and
      g==K(root_open['ratio_equations'][1]))
check('root_nonzero_polynomial_match',num==P(root_rank['necessary_numerator']))
check('root_restricted_polynomial_match',restricted==K(root_rank['restricted_numerator']))

# Independently derive the all-offdiagonal-zero Ricci obstruction. Keep all
# nine D_i h_j formal: no constant-h assumption is permitted.
S=PolynomialRing(QQ,names=['H1','H2','H3']+[f'd{i}{j}' for i in range(3) for j in range(3)])
H=S.gens()[:3]; derivatives=[S.gens()[3+3*i:6+3*i] for i in range(3)]
def eps(i,j,k):
    if len({i,j,k})<3:return 0
    return 1 if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)] else -1
def structure(k,i,j):return -eps(i,j,k)*H[k]
Gamma={(k,i,j):(structure(k,i,j)-structure(i,j,k)+structure(j,k,i))/2
       for k in range(3) for i in range(3) for j in range(3)}
def bare_derivative(i,q):return sum(q.derivative(H[j])*derivatives[i][j] for j in range(3))
def curvature(l,k,i,j):
    return bare_derivative(i,Gamma[l,j,k])-bare_derivative(j,Gamma[l,i,k])+sum(
        Gamma[m,j,k]*Gamma[l,i,m]-Gamma[m,i,k]*Gamma[l,j,m]-
        structure(m,i,j)*Gamma[l,m,k] for m in range(3))
Ric=[sum(curvature(k,i,k,i) for k in range(3)) for i in range(3)]
for i in range(3):
    j,k=(i+1)%3,(i+2)%3
    check('variable_h_diagonal_Ricci_%d'%i,Ric[i]==(H[i]**2-(H[j]-H[k])**2)/2)
# q2-q3=2(H2-H3)(H2+H3-H1), and q1=0 with H1!=0
# protects H2-H3. Then H1=H2+H3 and q1=4H2H3, a contradiction.
q=[2*v for v in Ric]
check('flat_diagonal_difference_factor',q[1]-q[2]==2*(H[1]-H[2])*(H[1]+H[2]-H[0]))
check('flat_diagonal_final_factor',q[0].subs({H[0]:H[1]+H[2]})==4*H[1]*H[2])
payload={'P':P,'K':K,'R':R,'four_tangency_matrix':M,'four_tangency_determinant':detM,
         'reconstructions':solution,'square_coefficients':[cy,cw],'Y':Y,'W':W,
         'V':V,'f':f,'g':g,'f_numerator':num,'f_denominator':denominator,
         'nonzero_h3z_coefficient':coefficient,'L':L,'restricted_numerator':restricted,
         'restricted_h0_constant':constant,'variable_h_Ricci':Ric,'checks':checks}
save(payload,str(OUT/'zero_pair_independent_reader.sobj'))
report={'status':'PASS_ALL_PARAMETER_ORIGINAL_SOURCE','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'checks':checks,'four_tangency_determinant':str(detM.factor()),
        'square_coefficients':[str(c.factor()) for c in (cy,cw)],
        'necessary_denominator':str(denominator.factor()),
        'f_h3z_coefficient':str(coefficient.factor()),
        'restricted_h0_constant':str(constant.factor()),
        'ratio_velocities':[[str(v.factor()) for v in row[:2]] for row in V],
        'V3_y':str(v3y),
        'scope':'Actual open Bcal=0,p13=p31=0,y*w!=0 branch on the old h*mu protected chart; not a pointwise algebraic exclusion.',
        'seconds':time.time()-start}
(OUT/'zero_pair_independent_reader.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2),flush=True)
