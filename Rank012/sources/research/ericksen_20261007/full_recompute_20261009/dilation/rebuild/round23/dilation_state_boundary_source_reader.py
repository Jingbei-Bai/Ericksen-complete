"""Independent componentwise replay of the full angular boundary source."""
from sage.all import *
from pathlib import Path
import pickle,json,time
OUT=Path(__file__).absolute().parent;START=time.time()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
cert=load(str(OUT/'dilation_state_homothety_boundary_carrier.sobj'))
checks={};weights=[1]*9+[2]*6+[3]*2
for mode in ['general','pzero','qzero','resonance']:
    old=load(str(OUT/('dilation_state_column_linear_module_'+mode+'.sobj')))
    K=old['field'];g=K.gens();a=g[0]
    if mode=='general':b,z=g[1:]
    elif mode=='resonance':
        z=g[1];b=a*(2*z/(1+z*z))**2+((1-z*z)/(1+z*z))**2
    else:b=g[1]
    base=old['base'];S=base[:17];v=vector(K,base[17:20]);ax=vector(K,base[20:23]);w=ax.cross_product(v)
    pi=matrix(K,3,3,S[:9]);Id=identity_matrix(K,3)
    def coef(c):return K(c(a=a,b=b)) if hasattr(c,'subs') else K(c)
    def poly_jac(poly):
        return vector(K,[sum(coef(c)*ex[j]*prod(S[k]**(ex[k]-(j==k)) for k in range(17)) for ex,c in poly.items() if ex[j]) for j in range(17)])
    FJ=[[poly_jac(z) for z in row] for row in src['Fraw']]
    L=[]
    for i in range(3):
        Wi=pi.row(i)-pi.trace()/2*Id.row(i);cols=[]
        for j in range(23):
            delta=vector(K,[int(k==j) for k in range(23)])
            dS=delta[:17];dv=vector(K,delta[17:20]);da=vector(K,delta[20:23]);dpi=matrix(K,3,3,dS[:9])
            dWi=dpi.row(i)-dpi.trace()/2*Id.row(i);dw=da.cross_product(v)+ax.cross_product(dv)
            out=[FJ[i][k].dot_product(vector(K,dS))+weights[k]*(v[i]*dS[k]+dv[i]*S[k]) for k in range(17)]
            out+=list(dw[i]*w+w[i]*dw-dWi.cross_product(v)-Wi.cross_product(dv))
            out+=list(-dWi.cross_product(ax)-Wi.cross_product(da))
            cols.append(vector(K,out))
        Li=matrix(K,cols).transpose();L.append(Li)
        checks[mode+'_raw_field_'+str(i)]=Li==old['ambient_jacobians'][i]
    rows=[vector(K,list(poly_jac(z))+[0]*6) for seq in [src['C3'],src['Q4'],src['basis'][:27]] for z in seq]
    rows.extend([vector(K,[0]*17+list(2*v)+[0]*3),vector(K,[0]*20+list(2*ax)),vector(K,[0]*17+list(ax)+list(v))])
    J0=matrix(K,rows);checks[mode+'_raw_compatibility_source']=J0==old['raw_source_jacobian']
    B=sum((ax[i]*L[i] for i in range(3)),matrix.zero(K,23))
    checks[mode+'_full_raw_augmented_source']=cert[mode]['raw']==B.stack(J0)
    # At epsilon=0 every E_i vanishes; epsilon block is v_i, with no
    # mixed entries.  H=B+epsilon(R+omega A) contributes only DB.
    N=cert[mode]['N'];Ne=block_matrix(K,[[N,matrix.zero(K,23,1)],[matrix.zero(K,1,9),matrix(K,[[1]])]])
    Je=block_matrix(K,[[cert[mode]['J'],matrix.zero(K,14,1)]])
    Le=[block_diagonal_matrix(L[i],matrix(K,[[v[i]]])) for i in range(3)]
    checks[mode+'_complete_augmented_closure']=all(Je*M*Ne==0 for M in Le)
    R=block_diagonal_matrix(cert[mode]['R'],matrix(K,[[1]]));A=block_diagonal_matrix(cert[mode]['A'],matrix(K,[[0]]))
    Pcu=block_diagonal_matrix(cert[mode]['Pcu'],matrix(K,[[1]]))
    checks[mode+'_extended_nonnegative_rank4']=Pcu.trace()==4
    checks[mode+'_extended_angular_zero']=A*Pcu==0
    checks[mode+'_extended_characteristic']=R.charpoly()==R.charpoly().parent()('x*(x+2)*(x-1)^3*(x+1)^5')
rec={'status':'PASS_FULL_PARAMETER_ORIGINAL_SOURCE_REPLAY','checks':checks,'all_pass':all(checks.values()),'count':len(checks),'seconds':time.time()-START}
(OUT/'dilation_state_boundary_source_reader.json').write_text(json.dumps(rec,indent=2))
print(json.dumps({k:rec[k] for k in ['status','all_pass','count','seconds']},indent=2),flush=True)
assert all(checks.values()),[n for n,z in checks.items() if not z]
