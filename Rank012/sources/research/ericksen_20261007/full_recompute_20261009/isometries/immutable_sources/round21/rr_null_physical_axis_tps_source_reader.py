"""Replay the axis blowup by EXACT rational tail valuations, not series recurrences.

All differences are checked to have valuation >=6 in e, independently of the
truncated-series arithmetic used to discover the coefficients.
"""
from sage.all import *
from pathlib import Path
import time,signal,json,resource
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(120)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
d=load(str(O/'rr_null_physical_axis_binary_tps_plus.sobj'))
src=load(str(O/'rr_null_physical_axis_reduction_checkpoint_plus.sobj'))
binary=load(str(O/'rr_null_physical_axis_binary_resultant_plus.sobj'))
F=d['field'];P=F.ring();s,t=F.gens();E=src['e_ring'];e=E.gen();FE=E.fraction_field();R=d['ring'];a,b=R.gens()
S=PolynomialRing(FE,names=('a','b','l'));aa,bb,ll=S.gens();ev=R.hom([aa,bb],S)
def lift(series):return sum(ev(f)*e**(2*j) for j,f in enumerate(series))
def val(z):
 z=FE(z)
 return z.numerator().valuation()-z.denominator().valuation() if z else Infinity
def tail(f):return all(val(c)>=6 for c in S(f).coefficients())
checks={}
N=src['reduced_matrix'];vals=[min(val(c) for c in row) for row in N.rows()]
N=diagonal_matrix(FE,[e**(-v) for v in vals])*N
G=lift(d['G']);u=lift(d['u']);v=lift(d['v']);A=lift(d['A']);B=lift(d['B']);C=lift(d['C']);M=lift(d['M']);NN=lift(d['N'])
checks['G_exact_tail']=tail(N[1,5]-G)
checks['u_exact_tail']=tail(2*N[1,5]*u-N[1,2])
checks['v_exact_tail']=tail(2*N[1,5]*v-N[1,4])
xx=ll-u*aa-v*bb;mons=[aa*aa,aa*bb,aa*xx,bb*bb,bb*xx,xx*xx]
forms=[sum(co*m for co,m in zip(row,mons)) for row in N.rows()]
checks['paired_no_x_square']=N[0,5]==0
checks['first_source_exact_tail']=tail(forms[0]-(M*ll+A))
checks['square_source_exact_tail']=tail(forms[1]-(G*ll*ll+B))
checks['third_source_exact_tail']=tail(forms[2]-(S(N[2,5])/G*(G*ll*ll+B)+NN*ll+C))
checks['D_ideal_identity_tail']=tail(lift(d['D'])-(M*C-NN*A))
checks['H_ideal_identity_tail']=tail(lift(d['H'])-(G*A*A+B*M*M))
checks['A_C_zero_initial']=d['A'][0]==0 and d['C'][0]==0
checks['B_two_zero_initial_orders']=d['B'][0]==0 and d['B'][1]==0
checks['D_exact_first_order']=d['D'][0]==0 and d['D'][1]!=0
checks['H_exact_second_order']=d['H'][0]==0 and d['H'][1]==0 and d['H'][2]!=0
expectedG=-(s*s-t*t)*(s*s-1)*(s*s+t*t-2)**5/(8*(s*s+t*t)**2*((s*s-1)**2+(t*t-1)**2))
checks['protected_G0_formula']=d['G'][0]==expectedG
checks['D_primitive_source']=d['D'][1]==binary['D_scale']*binary['D']
checks['H_primitive_source']=d['H'][2]==binary['H_scale']*binary['H']
U=PolynomialRing(P,'a')
ff=U([P(binary['H'].monomial_coefficient(a**j*b**(4-j))) for j in range(5)])
gg=U([P(binary['D'].monomial_coefficient(a**j*b**(3-j))) for j in range(4)])
checks['fixed_degree_homogeneous_resultant']=ff.resultant(gg)==binary['resultant']
checks['resultant_factorization']=prod(f**m for f,m in binary['factorization'])*binary['factorization'].unit()==binary['resultant']
assert all(checks.values()),checks
rep={'status':'PASS_COMPLETE_AXIS_TRUNCATION_SOURCE_REPLAY','checks':checks,'count':len(checks),'seconds':time.time()-st,
 'note':'Exact rational tails, protected G0, source ideal consequences, homogeneous binary resultant; combine with material-wall and target-reflection readers for theorem A.'}
(O/'rr_null_physical_axis_tps_source_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
