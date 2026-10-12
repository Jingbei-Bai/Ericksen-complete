"""All 37 full-source third jets with the genuine second range correction.

The preceding first tangent corrections are zero on the separately recorded
generic chart. Second kernel corrections k2,eta2, second tangents dl2,du2,dv2
and all third range corrections are retained, rather than fixed.
"""
from pathlib import Path
_here=Path(__file__).resolve().parent
_engine=(_here/'independent_ratio_incoming12_raw_third_range.py').read_text()
exec(compile(_engine.split("F=raw['ring']",1)[0],str(_here/'independent_ratio_incoming12_raw_third_range.py'),'exec'))
# Delay reduction by P(y) until a complete source has been added.  In a
# third jet y-degree never exceeds epsilon-degree, so degrees >3 are unused.
_n2=n2
zero=(K(0),)*8;one=(K(1),)+(K(0),)*7;ell=(K(0),K(1))+(K(0),)*6;yy=(K(0),K(0),K(1))+(K(0),)*5
def add(x,y):return tuple(x[i]+y[i] for i in range(8))
def mul(x,y):
 out=[K(0)]*8
 for i,ci in enumerate(x):
  if not ci:continue
  for j,cj in enumerate(y):
   if not cj:continue
   lp=i%2+j%2;yp=i//2+j//2
   if yp>3:continue
   out[lp%2+2*yp]+=ci*cj*(dd if lp>=2 else 1)
 return tuple(out)
def uv(f):
 f=U(f)%c['relation'];return (K(f[0]),K(f[1]))+(K(0),)*6
def final_reduce(v):
 return (v[0]+v[4]*rel0+v[6]*rel0*rel1,v[1]+v[5]*rel0+v[7]*rel0*rel1,v[2]+v[4]*rel1+v[6]*(rel0+rel1**2),v[3]+v[5]*rel1+v[7]*(rel0+rel1**2))
n2=[tuple(v)+(K(0),)*4 for v in _n2];k1=mul(ell,yy);p13=add(scale(kap,ell),scale(-1,k1))
states=[{0:one},{1:p13,2:n2[0]},{2:n2[1]},{2:n2[2]},{1:k1},{0:ell},{0:scale(-2,ell)},{},{0:scale(-2,one)}];cache={};moncache={}
rawmode=os.environ.get('INCOMING_KEEP_CUBIC','0')=='1'
scope='All 37 full-source third jets with genuine second range correction; second kernel and tangent corrections and all third range corrections retained.';prefix=p/('independent_ratio_incoming12_all_third_raw' if rawmode else 'independent_ratio_incoming12_all_third_range')
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hd=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));ed=load(str(p/'independent_ratio_even_derivative_cyclic_direct.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));old=load(str(p/'independent_ratio_incoming12_generic_nullwall_all.sobj'));si={int(it['source_row']):it for it in old['items']}
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+800*1024**2,)*2)
pool=raw['polynomial_pool'];O=PolynomialRing(K,names=tuple('nt_p%d'%j for j in range(6))+tuple('nt_R%d'%j for j in range(3)));ov=O.gens()
def cp(f):return O({tuple(map(int,e)):K(cc) for e,cc in f.dict().items()})
def put(name,f):pool[name]=f.dict() if hasattr(f,'dict') else f;return name
Q={qi:cp(oq['quartic_source_rows'][qi]) for qi in (0,1,5,6)};R=[cp(f) for f in g['six_full_sources']]
for qi,f in Q.items():put('allQ%d'%qi,f)
for j,f in enumerate(R):put('allR%d'%j,f)
for name,ddict in enumerate(hd['polynomial_pool']):put('allH:'+str(name),ddict)
def one3(ref,degree):
 out=zero
 for (d,ep),vv in factor(ref).items():
  if degree-d>=0 and ep+degree-d==3:out=add(out,vv)
 return out
def source3(row):
 if row<4:
  qi=(0,1,5,6)[row];return one3('allQ%d'%qi,int(Q[qi].total_degree()))
 if row<10:return one3('allR%d'%(row-4),4)
 if row<28:
  out=zero
  for sign,aa,bb in hd['records'][row-10]['full_signed_product_DAG']:
   out=add(out,scale(sign,product(factor('allH:'+str(aa)),factor('allH:'+str(bb)))))
  return out
 it=ed['items'][row-28];i=int(it['direction']);f=Q[it['Q_source_index']];out=zero
 for j in range(9):
  aa=put('allEF%d_%d'%(i,j),g['first_gauged_fields'][i,j]);bb=put('allEQ%d_%d'%(row,j),f.derivative(ov[j]));out=add(out,product(factor(aa),factor(bb)))
 for j,bp in enumerate((h,z)):
  aa=put('allEB%d_%d'%(i,j),g['base_fields'][i,j]);bb=put('allEQb%d_%d'%(row,j),{tuple(map(int,e)):K(cc).derivative(bp) for e,cc in f.dict().items()});out=add(out,product(factor(aa),factor(bb)))
 for pos,(co,(_,rp,sign)) in enumerate(zip(it['global_support_coefficient_dicts'],ed['supports'][i])):
  aa=put('allEco%d_%d'%(row,pos),{tuple(map(int,e)):K(cc)*g['source_unit_coefficients'][i] for e,cc in co.items()});out=add(out,scale(-sign,product(factor(aa),factor('allR%d'%rp))))
 return out
def atomic(obj,target):
 temp=Path(str(target)+'.tmp.sobj');save(obj,str(temp),compress=False);os.replace(temp,target)
out=[]
for row in range(37):
 target=Path(str(prefix)+'.row%d.sobj'%row)
 if target.exists():out.append(load(str(target)));print('resume',row,flush=True);continue
 third=source3(row);rhs=si[row]['rhs']
 coeff=[add(scale(2,mul(uv(rhs[0]),k1)),uv(rhs[2])),add(mul(uv(rhs[1]),k1),uv(rhs[4]))]
 for j in range(3):coeff.append(add(mul(uv(rhs[6+3*j]),k1),uv(rhs[8+3*j])))
 if not rawmode:third=final_reduce(third);coeff=[final_reduce(v) for v in coeff]
 item={'source_row':row,'constant':third,'unknown_coefficients':coeff,'third_range_coefficients':[(uv(f) if rawmode else final_reduce(uv(f))) for f in si[row]['linear_row']]};out.append(item);atomic(item,target)
 print('row',row,'constant-support',[j for j,v in enumerate(third) if v],'linear-count',sum(any(v) for v in coeff),'sec',round(time.monotonic()-st,2),flush=True)
atomic({'coefficient_field':K,'l_square':dd,'y_quadratic':P,'basis':(['1','l','y','l*y','y^2','l*y^2','y^3','l*y^3'] if rawmode else ['1','l','y','l*y']),'items':out,'unknowns':['k2','eta2','dl2','du2','dv2'],'second_range':_n2,'scope':scope},Path(str(prefix)+'.sobj'))
Path(str(prefix)+'.json').write_text(json.dumps({'rows':37,'seconds':time.monotonic()-st,'constant_support':{str(it['source_row']):[j for j,v in enumerate(it['constant']) if v] for it in out},'scope':scope},indent=2)+'\n');print('DONE',time.monotonic()-st,flush=True)
