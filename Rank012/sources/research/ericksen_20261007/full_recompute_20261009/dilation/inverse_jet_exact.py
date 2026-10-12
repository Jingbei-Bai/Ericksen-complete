"""General second-jet proof of the variable-Jacobian divergence identities.

Work on the dense LDU chart of GL(3); clear denominators to obtain identities
on every invertible jet. This is an exact chart argument, not sampling.
"""
from sage.all import *
from pathlib import Path
import json,time
HERE=Path(__file__).absolute().parent; start=time.time()
K=FractionField(PolynomialRing(QQ,names=('d0','d1','d2','l10','l20','l21','u01','u02','u12')))
d0,d1,d2,l10,l20,l21,u01,u02,u12=K.gens()
Lo=matrix(K,[[1,0,0],[l10,1,0],[l20,l21,1]])
Up=matrix(K,[[1,u01,u02],[0,1,u12],[0,0,1]])
F=Lo*diagonal_matrix(K,[d0,d1,d2])*Up; G=F.inverse();H=G*G.transpose()
checks={}; coefficients=[]
for comp in range(3):
 for j in range(3):
  for k in range(j,3):
   DF=[matrix.zero(K,3) for _ in range(3)]
   DF[k][comp,j]=1;DF[j][comp,k]=1
   DG=[-G*D*G for D in DF]
   DC=[D.transpose()*F+F.transpose()*D for D in DF]
   DH=[D*G.transpose()+G*D.transpose() for D in DG]
   divC=vector(K,[sum(DC[jj][i,jj] for jj in range(3)) for i in range(3)])
   divH=vector(K,[sum(DH[jj][i,jj] for jj in range(3)) for i in range(3)])
   deltaf=vector(K,[sum(DF[jj][ii,jj] for jj in range(3)) for ii in range(3)])
   trgrad=vector(K,[D.trace()/2 for D in DC])
   logJ=vector(K,[(G*D).trace() for D in DF])
   deltag=vector(K,[sum(G[kk,aa]*DG[kk][ii,aa] for kk in range(3) for aa in range(3)) for ii in range(3)])
   checks['divC_Hessian_%d%d%d'%(comp,j,k)]=divC==F.transpose()*deltaf+trgrad
   checks['divH_Hessian_%d%d%d'%(comp,j,k)]=divH==deltag-H*logJ
assert all(checks.values()),[n for n,v in checks.items() if not v]
rec={'status':'PASS_GENERAL_INVERTIBLE_SECOND_JET','checks':checks,'count':len(checks),
     'scope':'All Hessian coefficients, arbitrary three by three invertible derivative by dense-chart polynomial continuation; inverse identities follow by exchanging maps.',
     'division':'d0*d1*d2 on the computation chart only; denominator-cleared identity extends to det(F)!=0.',
     'seconds':time.time()-start}
(HERE/'inverse_jet_exact.json').write_text(json.dumps(rec,indent=2))
print(json.dumps({'status':rec['status'],'count':len(checks),'seconds':rec['seconds']}),flush=True)
