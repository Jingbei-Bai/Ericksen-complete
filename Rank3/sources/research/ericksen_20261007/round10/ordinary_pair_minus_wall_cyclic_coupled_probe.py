from sage.all import *
from pathlib import Path
import resource,os
p=Path(__file__).resolve().parent;resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
d=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'));P=d['ring'];ko=P.base_ring();K=PolynomialRing(QQ,names=('m_a','m_r')).fraction_field();a,r=K.gens();m=1+a*r;k=a*r;ph=ko.ring().hom([a,r,1,-1],K);T=PolynomialRing(K,'ll');l=T.gen();A=T.quotient(l*l+k/m,'lb');lb=A.gen()
def cv(f):
 out=A.zero()
 for ee,co in f.dict().items():out+=ph(co.numerator())/ph(co.denominator())*lb**int(ee)
 return out
M=d['matrix'].apply_map(cv);print('cols',d['columns'],flush=True);print('rank',M.rank(),flush=True);print('pivotrows',M.transpose().pivots(),flush=True);print('rightkernel',M.right_kernel().basis(),flush=True);print('range4rank',M[:,:4].rank(),flush=True);print('range4pivots',M[:,:4].transpose().pivots(),flush=True)
save({'ring':A,'matrix':M,'columns':d['columns'],'rank':M.rank(),'kernel':M.right_kernel().basis()},str(p/'ordinary_pair_minus_wall_cyclic_coupled_probe.sobj'),compress=False)
