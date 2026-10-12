"""Complete tangent carrier for B+epsilon(R+omega A)=0 at a cylinder.

This uses full-angular fields, not an axisymmetric perturbation ansatz.
The target angular rate and pressures are deliberately not fixed: imposing
their additional equations can only shrink the actual state subset.
"""
from sage.all import *
from pathlib import Path
import json,time,signal
OUT=Path(__file__).absolute().parent;START=time.time();signal.alarm(180)
def norm(x):
    if not x:return x
    lc=x.denominator().leading_coefficient()
    return x.parent()(x.numerator()/lc)/x.parent()(x.denominator()/lc)
def nm(M):return M.apply_map(norm)
def det_safe(M):
    W=[list(r) for r in M.rows()];n=M.nrows();d=M.base_ring()(1)
    for i in range(n):
        k=next((k for k in range(i,n) if W[k][i]),None)
        if k is None:return M.base_ring()(0)
        if k!=i:W[i],W[k]=W[k],W[i];d=-d
        z=W[i][i];d=norm(d*z)
        for j in range(i+1,n):
            if not W[j][i]:continue
            c=norm(W[j][i]/z)
            for k in range(i+1,n):W[j][k]=norm(W[j][k]-c*W[i][k])
            W[j][i]=M.base_ring()(0)
    return d
def protected(poly,guard):
    return bool(poly) and all(guard%f==0 for f,e in poly.factor())
checks={};records={};cert={}
for mode in ['general','pzero','qzero','resonance']:
    src=load(str(OUT/('dilation_state_column_linear_module_'+mode+'.sobj')))
    K=src['field'];P=K.ring();gg=K.gens();a=gg[0]
    if mode=='general':
        b,z=gg[1:];guard=P(a*b*(a-1)*(b-1)*(a-b)*z*(z-1)*(z+1)*(z*z+1))
    elif mode=='resonance':
        z=gg[1];b=a*(2*z/(1+z*z))**2+((1-z*z)/(1+z*z))**2
        guard=P(a*(a-1)*z*(z-1)*(z+1)*(z*z+1)*b.numerator())
    else:
        b=gg[1];guard=P(a*b*(a-1)*(b-1)*(a-b))
    L=list(map(nm,src['ambient_jacobians']));base=src['base'];v=vector(K,base[17:20]);ax=vector(K,base[20:23]);w=ax.cross_product(v)
    Ra=nm(sum((v[i]*L[i] for i in range(3)),matrix.zero(K,23)))
    Aa=nm(sum((w[i]*L[i] for i in range(3)),matrix.zero(K,23)))
    Ba=nm(sum((ax[i]*L[i] for i in range(3)),matrix.zero(K,23)))
    raw=nm(Ba.stack(src['raw_source_jacobian']))
    print(mode,'raw ready',flush=True)
    J=nm(raw.row_space().basis_matrix())
    cc=J.pivots();free=[j for j in range(23) if j not in cc]
    N=matrix.zero(K,23,len(free))
    for k,j in enumerate(free):
        N[j,k]=1
        for i,col in enumerate(cc):N[col,k]=-J[i,j]
    N=nm(N)
    print(mode,'RREF ready',flush=True)
    rr=raw.transpose().pivots();mi=det_safe(raw.matrix_from_rows_and_columns(rr,cc))
    print(mode,'minor ready',flush=True)
    checks[mode+'_dimension_9']=N.ncols()==9
    checks[mode+'_necessary_sources_zero']=raw*N==0
    checks[mode+'_minor_protected']=protected(mi.numerator(),guard) and protected(mi.denominator(),guard)
    checks[mode+'_basis_denominators_protected']=all(protected(x.denominator(),guard) for M in [N,J] for x in M.list())
    checks[mode+'_basis_identity_minor']=N.matrix_from_rows(free)==identity_matrix(K,9)
    checks[mode+'_complete_linear_differential_closure']=all(J*Li*N==0 for Li in L)
    R,A,B=[nm(N.solve_right(V*N)) for V in [Ra,Aa,Ba]]
    I=identity_matrix(K,9)
    checks[mode+'_B_zero']=B==0
    checks[mode+'_R_characteristic']=R.charpoly()==R.charpoly().parent()('x*(x+2)*(x-1)^2*(x+1)^5')
    checks[mode+'_R_annihilator']=R*(R-I)*(R+2*I)*(R+I)**2==0
    checks[mode+'_RA_commute']=R*A==A*R
    checks[mode+'_A_annihilator']=A*A*(A*A+I)==0
    # The R primary summand for eigenvalues 0 and 1 is obtained without any
    # parameter-dependent eigenbasis: a universal Hermite spectral projector.
    T=PolynomialRing(QQ,'x');x=T.gen();pplus=x*(x-1);pneg=(x+2)*(x+1)**2
    _,inv,_=xgcd(pneg,pplus);projector=(inv*pneg)%(pplus*pneg)
    Pcu=nm(sum((K(c)*R**i for i,c in enumerate(projector.list())),matrix.zero(K,9)))
    checks[mode+'_Pcu_idempotent']=Pcu*Pcu==Pcu
    checks[mode+'_Pcu_rank3']=Pcu.trace()==3
    checks[mode+'_A_center_unstable_zero']=A*Pcu==0
    checks[mode+'_R_positive_semisimple']=R*(R-I)*Pcu==0
    # A stable row-space basis can have apparent pivot factors.  Record them
    # honestly instead of promoting a generic-field rank calculation.
    records[mode]={'minor':str(mi.factor()),'minor_rows':list(rr),'minor_columns':list(cc),'denominator_factors':sorted(set(str(f) for M in [N,J] for y in M.list() if y for f,e in y.denominator().factor())),
       'R_charpoly':str(R.charpoly().factor()),'A_charpoly':str(A.charpoly().factor()),'Pcu':str(projector)}
    cert[mode]={'raw':raw,'J':J,'N':N,'ambient_R':Ra,'ambient_A':Aa,'ambient_B':Ba,'R':R,'A':A,'B':B,'Pcu':Pcu,'minor':mi,'minor_rows':rr,'minor_columns':cc,'guard':guard}
    print(mode,records[mode],flush=True)
rec={'status':'BOUNDARY_CARRIER_LINEAR_SOURCE_CERTIFICATE','checks':checks,'all_pass':all(checks.values()),'records':records,'seconds':time.time()-START}
save(cert,str(OUT/'dilation_state_homothety_boundary_carrier.sobj'))
(OUT/'dilation_state_homothety_boundary_carrier.json').write_text(json.dumps(rec,indent=2))
print('all',rec['all_pass'],len(checks),'seconds',rec['seconds'],flush=True)
assert all(checks.values()),[n for n,z in checks.items() if not z]
