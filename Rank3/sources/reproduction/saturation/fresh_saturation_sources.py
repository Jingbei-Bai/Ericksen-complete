"""Full polynomial source presentations for the inactive incidence charts."""
from sage.all import *
from pathlib import Path
import argparse,hashlib,json,time

parser=argparse.ArgumentParser()
parser.add_argument('--work-root',type=Path,required=True)
parser.add_argument('--output-dir',type=Path)
args=parser.parse_args(); root=args.work_root.resolve()
source_dir=root/'research/ericksen_20261007/full_recompute_20261009/rank3'
out=args.output_dir.resolve() if args.output_dir else root/'reproduction/saturation'
assert out.is_relative_to(root), 'Saturation output must be inside the work directory'
out.mkdir(parents=True,exist_ok=True); checks=[]; started=time.monotonic()
def ck(name,value):
    assert value,name
    checks.append(name)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
inputs=[source_dir/'fresh_global_polynomials.sobj',source_dir/'fresh_joint_fields.sobj']
print('Load complete global and joint sources',flush=True)
g,j=[load(str(p)) for p in inputs]
K=PolynomialRing(QQ,names=('ia','ir','iH','iZ')).fraction_field()
a,r,H,Z=K.gens();q=1+r;m=1+a*r
old=g['ring'].base_ring(); parameter_map=(a,a*q/m,H/r,q*Z/r)
ph=old.ring().hom(parameter_map,K)
R=PolynomialRing(K,names=('kp12','kp13','kp21','kp23','kp31','kp32','nv0','nv1','nv2'))
x=R.gens();vv=vector(R,x[6:])
def coefficient(c):
    c=old(c)
    return ph(c.numerator())/ph(c.denominator())
def polynomial(f):return R({tuple(map(int,e)):coefficient(c) for e,c in f.dict().items()})
def face(f,i,j):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(e[:6])==i and sum(e[6:])==j})
print('Adapt six R sources, Q0 and Q1,Q5,Q6',flush=True)
gad={key:matrix(R,value.nrows(),value.ncols(),[polynomial(v) for v in value.list()]) for key,value in g.items() if key in ('L','C','D')}
gad.update(Q0=polynomial(g['Q0']),sources=[polynomial(f) for f in g['sources']],offset=vector(R,[polynomial(f) for f in g['offset']]))
jad={'quartic_source_rows':{qi:polynomial(j['quartic_source_rows'][qi]) for qi in (1,5,6)}}
L,C,D= gad['L'],gad['C'],gad['D'];offset=gad['offset'];Q0=gad['Q0'];sources=gad['sources'];ell=L*vv
ck('all complete R sources',len(sources)==6)
alpha=(Z-H)/(Z*(H+1)); beta=H*(Z-1)/(Z*(H+1))
expected_L=matrix(R,[[x[0],0,-alpha*x[5]],[x[1],-alpha/beta*x[3],0],[0,x[2],-beta*x[4]]])
ck('original incidence matrix',L==expected_L)
support=x[4]*x[2]*x[0]
W=matrix(R,3,6,0);W[0,0]=x[2]*x[0];W[1,2]=x[4]*x[0];W[2,3]=x[4]*x[2]
ck('support identity',W*D==support*identity_matrix(R,3))
B=[]
for pindex,rindex in ((0,6),(1,6),(2,7)):
    terms={}
    for e,c in face(Q0,1,2).dict().items():
        if e[pindex]:
            ck('Q0 cofactor exact coordinate division',e[pindex]==1 and e[rindex]>=1)
            en=list(map(int,e));en[pindex]-=1;en[rindex]-=1;terms[tuple(en)]=c
    B.append(R(terms))
B=vector(R,B)
ck('complete Q0 top module',face(Q0,1,2)==B*ell)
ck('all complete R31 modules',vector(R,[face(f,3,1) for f in sources])==D*C*ell)
ck('all complete R40 offsets',vector(R,[face(f,4,0) for f in sources])==D*offset)
det=C.det();adj=C.adjugate()
initial=support*(det*face(Q0,2,1)-B*adj*offset)
ck('whole top cancellation',support*det*face(Q0,1,2)-B*adj*W*vector(R,[face(f,3,1) for f in sources])==0)
ck('whole initial identity',initial==support*det*face(Q0,2,1)-B*adj*W*vector(R,[face(f,4,0) for f in sources]))
bidegrees={(sum(e[:6])+6,sum(e[6:])) for e in Q0.dict()}
for f in sources:bidegrees.update((sum(e[:6])+4,sum(e[6:])+1) for e in f.dict())
for i,jj in bidegrees-{(7,2),(8,1)}:
    endpoint1=9-i-jj;endpoint2=QQ(19)/2-i-QQ(3)*jj/2
    ck('complete subcritical degree inequality '+str((i,jj)),endpoint1>=0 and endpoint2>=0 and (endpoint1>0 or endpoint2>0))
q_payload={'ring':R,'field':K,'full_Q0':Q0,'full_R_sources':sources,'D':D,'W':W,'support_denominator':support,
           'C':C,'C_determinant':det,'L':L,'Q0_L_cofactors':B,'subcritical_initial':initial,
           'source_DAG':'m*det(C)*Q0 - B^T*adj(C)*W*R_full','all_full_DAG_bidegrees':sorted(bidegrees)}
o_payload={'ring':R,'field':K,'L':L,'C':C,'offset':offset}
print('Q0 and six R polynomial identities verified',flush=True)
# The even source presentation retains the full polynomial identity.

def build_even_sources(fg_adapted, joint_adapted, K, R, checks):
    """Return complete adapted Q1,Q5,Q6 and their polynomial L cofactors."""
    from sage.all import matrix, vector

    def ck(name, value):
        if not bool(value):
            raise AssertionError(name)
        checks.append(name)

    ck('even coefficient ring and nine ordered state variables',
       R.base_ring() == K and R.ngens() == 9)
    x = R.gens()
    L = matrix(R, fg_adapted['L'])
    C = matrix(R, fg_adapted['C'])
    D = matrix(R, fg_adapted['D'])
    offset = vector(R, fg_adapted['offset'])
    ck('even fresh range/source matrix dimensions',
       L.nrows() == L.ncols() == C.nrows() == C.ncols() == 3
       and D.nrows() == 6 and D.ncols() == 3 and len(offset) == 3)
    ell = L * vector(R, x[6:])
    detL = L.det()
    ck('even range determinant is a nonzero polynomial', detL != 0)
    adj_transpose = L.adjugate().transpose()

    def face(f, pd, rd):
        return R({tuple(map(int, e)): c for e, c in R(f).dict().items()
                  if sum(e[:6]) == pd and sum(e[6:]) == rd})

    def exact_polynomial_division(num, den):
        # Lexicographic monomial division over K; no state denominator remains.
        work = {tuple(map(int, e)): K(c) for e, c in R(num).dict().items() if c}
        dd = {tuple(map(int, e)): K(c) for e, c in R(den).dict().items() if c}
        if not dd:
            raise ZeroDivisionError('zero polynomial divisor')
        lead = max(dd)
        lc = dd[lead]
        quotient = {}
        remainder = {}
        while work:
            e = max(work)
            c = work[e]
            if all(e[i] >= lead[i] for i in range(9)):
                qe = tuple(e[i] - lead[i] for i in range(9))
                qc = c / lc
                quotient[qe] = quotient.get(qe, K.zero()) + qc
                for ee, cc in dd.items():
                    key = tuple(ee[i] + qe[i] for i in range(9))
                    value = work.get(key, K.zero()) - qc * cc
                    if value:
                        work[key] = value
                    else:
                        work.pop(key, None)
            else:
                remainder[e] = work.pop(e)
        return R(quotient), R(remainder)

    quartics = joint_adapted['quartic_source_rows']
    items = []
    for qi, component in ((1, 2), (5, 0), (6, 1)):
        full = R(quartics[qi])
        q22 = face(full, 2, 2)
        square = ell[component] ** 2
        square_dict = {tuple(map(int, e)): c for e, c in square.dict().items()}
        ck('Q%d range square is nonzero' % qi, bool(square_dict))
        exponent = max(square_dict)
        q22_dict = {tuple(map(int, e)): c for e, c in q22.dict().items()}
        unit = K(q22_dict.get(exponent, K.zero()) / square_dict[exponent])
        ck('Q%d exact full square coefficient' % qi,
           unit != 0 and q22 == unit * square)

        q31 = face(full, 3, 1)
        qv = vector(R, [q31.derivative(v) for v in x[6:]])
        ck('Q%d complete R-linear cubic face' % qi,
           q31 == sum((qv[j] * x[6+j] for j in range(3)), R.zero())
           and all(all(sum(e[6:]) == 0 for e in f.dict()) for f in qv))
        numerators = adj_transpose * qv
        cofactors = []
        for j, numerator in enumerate(numerators):
            quotient, remainder = exact_polynomial_division(numerator, detL)
            ck('Q%d exact polynomial adjugate division %d' % (qi, j),
               remainder == 0 and numerator == quotient * detL)
            ck('Q%d pure p quadratic cofactor %d' % (qi, j),
               all(sum(e[:6]) == 2 and sum(e[6:]) == 0
                   for e in quotient.dict()))
            cofactors.append(quotient)
        AA = vector(R, cofactors)
        ck('Q%d complete cubic-linear L identity' % qi,
           q31 == AA * ell)
        pure_p4 = face(full, 4, 0)
        residual = full - unit * square - AA * ell - pure_p4
        for e in residual.dict():
            pd, rd = sum(e[:6]), sum(e[6:])
            # Twice the second endpoint bound avoids a Python float.
            ck('Q%d complete lower subcritical layer %s' % (qi, tuple(e)),
               4 - pd - rd >= 0 and 8 - 2 * pd - 3 * rd >= 0
               and (4 - pd - rd > 0 or 8 - 2 * pd - 3 * rd > 0))
        items.append({
            'Q_index': qi,
            'square_component': component,
            'full_source': full,
            'full_Q': full,
            'square_unit': unit,
            'A_p_quadratic': AA,
            'pure_p4': pure_p4,
            'initial_DAG': 'm^2*(U*(adj(C)*d)[j]^2-det(C)*A*(adj(C)*d)+det(C)^2*Q_p4)',
            'source_DAG': 'Let T=adj(C)*W*R_full and M=m*det(C). Source is M^2*Q-2*U*M*ell[j]*T[j]+U*T[j]^2-M*A*T.',
        })
    ck('three full even sources in original Q order',
       [(it['Q_index'], it['square_component']) for it in items]
       == [(1, 2), (5, 0), (6, 1)])
    support_denominator = x[4] * x[2] * x[0]
    W = matrix(R, 3, 6, 0)
    W[0, 0] = x[2] * x[0]
    W[1, 2] = x[4] * x[0]
    W[2, 3] = x[4] * x[2]
    ck('even polynomial support inverse',
       W * D == matrix(R, 3, 3,
                       lambda i, j: support_denominator if i == j else 0))
    return {
        'ring': R, 'field': K, 'L': L, 'ell': ell, 'C': C,
        'C_determinant': C.det(), 'D': D, 'W': W,
        'support_denominator': support_denominator,
        'full_R_sources': [R(f) for f in fg_adapted['sources']],
        'offset': offset, 'items': items,
        'checks': list(checks),
        'source_files': ['fresh_joint_fields.sobj', 'fresh_global_polynomials.sobj'],
        'scope': 'Complete adapted original even polynomials and exact polynomial L cofactors; no historical source array or state localization.',
    }


print('Derive complete even polynomial cofactors',flush=True)
e_payload=build_even_sources(gad,jad,K,R,checks)
# The complete fraction-free source identity is a universal polynomial
# identity, so it applies to the full Q, not only to its leading layers.
UR=PolynomialRing(QQ,names=('M','u','Q','ell0','ell1','ell2','T0','T1','T2','A0','A1','A2'))
MM,uu,QQfull=UR.gens()[:3]
LL=vector(UR,UR.gens()[3:6]);TT=vector(UR,UR.gens()[6:9]);AAformal=vector(UR,UR.gens()[9:12])
ZZformal=MM*LL-TT
for i,it in enumerate(e_payload['items']):
    c=int(it['square_component']);AA=it['A_p_quadratic'];Q=it['full_source'];u=it['square_unit']
    Fformal=MM**2*QQfull-2*uu*MM*LL[c]*TT[c]+uu*TT[c]**2-MM*AAformal*TT
    rewrittenformal=uu*ZZformal[c]**2+MM*AAformal*ZZformal+MM**2*(QQfull-uu*LL[c]**2-AAformal*LL)
    ck('universal complete even source multiplyback '+str(i),Fformal==rewrittenformal)
    ck('whole even p4 initial '+str(i),face(Q-u*ell[c]**2-AA*ell,4,0)==it['pure_p4'])
payloads={
 'independent_delta_subcritical_Q0_full_source.sobj':q_payload,
 'independent_delta_Q0_kernel_compatibility.sobj':o_payload,
 'independent_delta_subcritical_even_full_source.sobj':e_payload}
print('Save three source presentations',flush=True)
outputs=[]
for name,data in payloads.items():
    target=out/name;temp=out/(name+'.tmp.sobj')
    save(data,str(temp),compress=False);temp.replace(target)
    outputs.append({'path':str(target.relative_to(root)),'sha256':digest(target),'keys':sorted(data)})
report={'status':'PASS','checks':checks,'count':len(checks),'seconds':time.monotonic()-started,
        'inputs':[{'path':str(p.relative_to(root)),'sha256':digest(p)} for p in inputs],
        'outputs':outputs,'even_count':3}
(out/'source_identities.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS full saturation source identities',len(checks),flush=True)
