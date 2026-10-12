"""Execute the preserved reader 324 bytes, using full exact rational arithmetic.

The DAG backend retains complete finite numerator and denominator polynomials
over F=Q(i)(s,t), without expansion until a coefficient/zero query is needed.
The optional flat backend is the canonical isomorphism F(e)=Q(i)(s,t,e).
No e-adic truncation, AST edit, skipped assertion, or numerical specialization
is used. Only this source SHA is accepted. Runtime resource alarms are disabled
and replaced by the stated RSS monitor. The original source report remains the
primary assertion result; this runner records bytes, inputs and all 18 results.
"""
from sage.all import *
import sage.all as sg
from sage.rings.fraction_field import FractionField_generic
from sage.rings.fraction_field_element import FractionFieldElement
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from pathlib import Path
from functools import lru_cache
import argparse,ast,hashlib,importlib,json,os,resource,shutil,signal,sys,threading,time,traceback
ORIGINAL_SHA='ae8706b3814c6fc775571ff730e9c6837f8e9ba904f4843d26a365fcb4b16a5a'
NAMES=['rr_null_physical_axis_binary_tps_plus.sobj','rr_null_physical_axis_reduction_checkpoint_plus.sobj','rr_null_physical_axis_binary_resultant_plus.sobj']
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bundle-root',type=Path,required=True)
p.add_argument('--input-dir',type=Path,default=None)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--backend',choices=['dag','flat','native'],default='dag')
p.add_argument('--rss-limit-mib',type=int,default=4096)
args=p.parse_args();assert sys.flags.optimize==0,'Assertions must be enabled'
iso=args.bundle_root/'research'/'ericksen_20261007'/'full_recompute_20261009'/'isometries'
source=iso/'immutable_sources'/'round21'/'rr_null_physical_axis_tps_source_reader.py'
input_dir=args.input_dir or iso/'replay'/'round21'
O=args.output_dir;O.mkdir(parents=True,exist_ok=True)
source_bytes=source.read_bytes();assert hashlib.sha256(source_bytes).hexdigest()==ORIGINAL_SHA
executed=O/source.name;executed.write_bytes(source_bytes)
input_info={}
for name in NAMES:
 original=input_dir/name;target=O/name
 digest=hashlib.sha256(original.read_bytes()).hexdigest()
 shutil.copyfile(original,target)
 assert hashlib.sha256(target.read_bytes()).hexdigest()==digest
 input_info[name]={'fresh_input_path':str(original),'executed_input_path':str(target),'sha256':digest}
expected=[]
for node in ast.walk(ast.parse(source_bytes)):
 if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Subscript):
  target=node.targets[0]
  if isinstance(target.value,ast.Name) and target.value.id=='checks' and isinstance(target.slice,ast.Constant):expected.append(target.slice.value)
assert len(expected)==18
assert len([node for node in ast.walk(ast.parse(source_bytes)) if isinstance(node,ast.Assert)])==1
st=time.time();checks={};last_line=None;peak_rss=0;finished=False;load_records=[];backend_contract={};backend_info={}
def rss():return int(Path('/proc/self/statm').read_text().split()[1])*os.sysconf('SC_PAGE_SIZE')
def emit(event,status='RUNNING',extra=None):
 rep={'status':status,'event':event,'seconds':time.time()-st,'checks':dict(checks),'checks_completed':len(checks),'last_original_line':last_line,'peak_rss_bytes':peak_rss,'backend':args.backend,'backend_contract':backend_contract,'backend_info':backend_info,'loads':load_records,'original_source_sha256':ORIGINAL_SHA,'executed_source_sha256':hashlib.sha256(executed.read_bytes()).hexdigest(),'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'inputs':input_info}
 if extra:rep.update(extra)
 (O/'progress.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
 print(json.dumps({'event':event,'status':status,'seconds':rep['seconds'],'checks_completed':len(checks),'last_original_line':last_line,'peak_rss_mib':peak_rss//1024**2}),flush=True)
 return rep
def monitor():
 global peak_rss
 while not finished:
  peak_rss=max(peak_rss,rss())
  if args.rss_limit_mib and peak_rss>args.rss_limit_mib*1024**2:
   emit('outer_rss_limit_exceeded','RESOURCE_ABORT_NOT_MATHEMATICAL_FAILURE')
   os.kill(os.getpid(),signal.SIGTERM)
  time.sleep(.5)
threading.Thread(target=monitor,daemon=True).start()
# These are process-only runtime controls. The original source bytes retain
# their own resource calls, and the report discloses how they were mediated.
old_setrlimit=resource.setrlimit
resource.setrlimit=lambda kind,limits: None if kind==resource.RLIMIT_AS else old_setrlimit(kind,limits)
old_alarm=signal.alarm;signal.alarm=lambda seconds:old_alarm(0)
na=importlib.import_module('cysignals.alarm')
def disabled_native_alarm(seconds):na.cancel_alarm()
na.alarm=disabled_native_alarm
if hasattr(sg,'alarm'):sg.alarm=disabled_native_alarm
na.cancel_alarm();old_alarm(0)
loaded={}
class EValPolynomial:
 """Restricted numerator proxy: only reader324's no-argument valuation use.

 Native polynomial values and arithmetic remain inside the Cython fraction
 element. This is not advertised as a universal polynomial replacement.
 """
 def __init__(self,poly):self.poly=poly
 def valuation(self,*args,**kwds):
  if args or kwds:return self.poly.valuation(*args,**kwds)
  return min(m[2] for m in self.poly.dict()) if self.poly else Infinity
 def __getattr__(self,name):return getattr(self.poly,name)
 def _sage_(self):return self.poly
 def _polynomial_(self,R):return R(self.poly)
 def __repr__(self):return repr(self.poly)
class FlatFraction(FractionFieldElement):
 def numerator(self):return EValPolynomial(FractionFieldElement.numerator(self))
 def denominator(self):return EValPolynomial(FractionFieldElement.denominator(self))
def raw_num(z):return FractionFieldElement.numerator(z)
def raw_den(z):return FractionFieldElement.denominator(z)
def install_flat_backend():
 d=loaded[NAMES[0]];src=loaded[NAMES[1]];E=src['e_ring'];F=E.base_ring();P=F.ring();K=P.base_ring();oldFE=E.fraction_field()
 assert F==d['field']
 assert P.variable_names()==('s','t') and E.variable_name()=='e'
 T=PolynomialRing(K,names=('s','t','e'));ss,tt,ee=T.gens();FF=FractionField_generic(T,element_class=FlatFraction)
 def mapP(z):return T({(m[0],m[1],0):c for m,c in P(z).dict().items()})
 @lru_cache(maxsize=512)
 def mapF(z,*den):
  z=F(z);out=FF(mapP(z.numerator()))/FF(mapP(z.denominator()))
  if den:
   assert len(den)==1;return out/mapF(den[0])
  return out
 def mapE(z,*den):
  z=E(z);out=sum((mapF(c)*FF(ee)**j for j,c in z.dict().items()),FF.zero())
  if den:
   assert len(den)==1;return out/mapE(den[0])
  return out
 @lru_cache(maxsize=128)
 def mapFE(z,*den):
  z=oldFE(z);out=mapE(z.numerator())/mapE(z.denominator())
  if den:
   assert len(den)==1;return out/mapFE(den[0])
  return out
 for domain,fun in [(F,mapF),(E,mapE),(oldFE,mapFE)]:FF.register_coercion(SetMorphism(Hom(domain,FF,Rings()),fun))
 s,t=F.gens();e=E.gen()
 zz=oldFE((s*s+t*t)*e**6+(s/t)*e**7)/oldFE((s-t)*e**2+e**3)
 fzz=FF(zz);literal=((FF(ss)**2+FF(tt)**2)*FF(ee)**6+(FF(ss)/FF(tt))*FF(ee)**7)/(FF(ss-tt)*FF(ee)**2+FF(ee)**3)
 backend_contract['explicit_full_formula_transport']=raw_num(fzz)*raw_den(literal)==raw_num(literal)*raw_den(fzz)
 backend_contract['no_argument_e_valuation_preserved']=zz.numerator().valuation()-zz.denominator().valuation()==fzz.numerator().valuation()-fzz.denominator().valuation()
 backend_contract['F_generators']=FF(s)==FF(ss) and FF(t)==FF(tt)
 backend_contract['E_generator']=FF(e)==FF(ee)
 backend_contract['zero_valuation']=FF.zero().numerator().valuation()==Infinity
 backend_contract['mapped_E_zero_parent']=FF(E.zero()).parent() is FF and type(FF(E.zero())) is FlatFraction
 backend_contract['arithmetic_type_preserved']=type(fzz+fzz*fzz) is FlatFraction
 E.fraction_field=lambda:FF
 backend_contract['negative_e_power']=e**(-2)==FF(1)/FF(ee**2)
 R=d['ring'];a,b=R.gens();S=PolynomialRing(FF,names=('a','b','l'));aa,bb,ll=S.gens();ev=R.hom([aa,bb],S)
 backend_contract['actual_saved_R_hom']=ev(a*b*s)==aa*bb*FF(ss)
 for label,val in [('parameter_common_factor',oldFE((s-t)*(e**6+e**7))/oldFE((s-t)*e**2)),('e_common_factor',oldFE(e**5*(s+e))/oldFE(e**3*(t+e)))]:
  fv=FF(val);backend_contract[label]=val.numerator().valuation()-val.denominator().valuation()==fv.numerator().valuation()-fv.denominator().valuation()
 small=matrix(oldFE,[[zz,zz+1],[zz-1,zz*zz]])
 transported=diagonal_matrix(FF,[e**4,e**(-2)])*small
 backend_contract['matrix_transport_parent']=transported.base_ring() is FF
 backend_contract['matrix_transport_full_entry']=transported[1,1]==FF(1)/FF(ee**2)*FF(zz*zz)
 assert all(backend_contract.values()),backend_contract
 backend_info.update({'old_field':str(oldFE),'new_field':str(FF),'canonical_isomorphism':'Q(i)(s,t)(e) -> Q(i)(s,t,e)','full_rational_function_arithmetic':True,'series_truncation_used':False,'restricted_proxy_scope':'Only this original reader source SHA, which calls numerator/denominator only through their no-argument e valuation. Native full polynomial arithmetic is retained internally.'})
 emit('full_exact_backend_installed')
old_load=sg.load

def install_dag_backend():
 import reader324_full_dag_backend_v2 as dag_backend
 CompleteFractionField=dag_backend.CompleteFractionField
 actual_module=Path(dag_backend.__file__).resolve()
 expected_module=Path(__file__).resolve().with_name('reader324_full_dag_backend_v2.py')
 assert actual_module==expected_module,(actual_module,expected_module)
 d=loaded[NAMES[0]];src=loaded[NAMES[1]];E=src['e_ring'];F=E.base_ring();oldFE=E.fraction_field()
 assert F==d['field'] and E.variable_name()=='e'
 LF=CompleteFractionField(E,oldFE);s,t=F.gens();e=E.gen()
 literal=oldFE((s*s+t*t)*e**7+(s/t)*e**37)/oldFE((s-t)*e**2+e**3)
 z=LF(literal)
 backend_contract['explicit_full_formula_transport']=z.native()==literal
 backend_contract['exact_full_e_valuation']=z.numerator().valuation()-z.denominator().valuation()==literal.numerator().valuation()-literal.denominator().valuation()
 backend_contract['high_degree_nonzero_is_retained']=LF(e**123).numerator().valuation()==123
 backend_contract['full_polynomial_zero_scan']=not(z-z)
 backend_contract['coefficient_field_constants']=LF(s).native()==oldFE(s) and LF(t).native()==oldFE(t)
 backend_contract['complex_constant_retained']=LF(F.ring().base_ring().gen()).native()==oldFE(F.ring().base_ring().gen())
 backend_contract['mapped_E_zero_parent']=LF(E.zero()).parent() is LF and not LF(E.zero())
 backend_contract['exact_zero_valuation']=LF.zero().numerator().valuation()==Infinity
 backend_contract['old_field_all_saved_entries']=all(LF(c).native()==oldFE(c) for c in src['reduced_matrix'].list())
 E.fraction_field=lambda:LF
 backend_contract['negative_e_power']=e**(-2)==LF.one()/LF(e**2)
 R=d['ring'];a,b=R.gens();S=PolynomialRing(LF,names=('a','b','l'));aa,bb,ll=S.gens();ev=R.hom([aa,bb],S)
 backend_contract['actual_saved_R_hom']=ev(a*b*s)==aa*bb*LF(s)
 small=matrix(oldFE,[[literal,literal+1],[literal-1,literal*literal]])
 transported=diagonal_matrix(LF,[e**4,e**(-2)])*small
 backend_contract['mixed_matrix_transport_parent']=transported.base_ring() is LF
 backend_contract['mixed_matrix_transport_entry']=transported[1,1].native()==oldFE.one()/oldFE(e**2)*literal*literal
 assert all(backend_contract.values()),backend_contract
 backend_info.update({'old_field':str(oldFE),'new_field':str(LF),'canonical_isomorphism':'Frac(F[e]) represented by complete finite polynomial DAG fractions; F=Q(i)(s,t)','full_rational_function_arithmetic':True,'series_truncation_used':False,'fixed_truncation_order':None,'complete_polynomial_degree_bound':True,'zero_and_valuation_query':'Exact F coefficient scan to first nonzero or through the full finite degree upper bound; no precision argument','restricted_backend_scope':'Only original reader SHA '+ORIGINAL_SHA,'backend_module_actual_import_path':str(actual_module),'backend_module_sha256':hashlib.sha256(actual_module.read_bytes()).hexdigest(),'denominator_management':'Common structurally identical complete polynomial factors and multiplicities; exact common denominator and identical-factor cancellation only'})
 emit('full_exact_dag_backend_installed')
def tracked_load(filename,*a,**kw):
 path=Path(filename)
 assert path.name in NAMES and path.parent.resolve()==O.resolve(),path
 digest=hashlib.sha256(path.read_bytes()).hexdigest();assert digest==input_info[path.name]['sha256']
 obj=old_load(str(path),*a,**kw);loaded[path.name]=obj
 load_records.append({'path':str(path),'sha256':digest,'from_fresh_copy':True})
 if len(loaded)==3 and args.backend=='flat':install_flat_backend()
 if len(loaded)==3 and args.backend=='dag':install_dag_backend()
 return obj
sg.load=tracked_load
executed_name=str(executed)
def trace(frame,event,arg):
 global checks,last_line
 if frame.f_code.co_filename==executed_name and frame.f_code.co_name=='<module>' and event=='line':
  last_line=frame.f_lineno
  observed=frame.f_locals.get('checks')
  if isinstance(observed,dict):checks=dict(observed)
  emit('original_line_'+str(last_line))
 return trace
namespace={'__name__':'__main__','__file__':executed_name}
status='EXECUTION_FAILURE';error=None;original_report=None
try:
 emit('begin_original_source_byte_exact_execution')
 sys.settrace(trace)
 exec(compile(source_bytes,executed_name,'exec'),namespace)
 sys.settrace(None)
 checks=dict(namespace['checks'])
 original_report=json.loads((O/'rr_null_physical_axis_tps_source_reader.json').read_text())
 assert list(checks)==expected and len(checks)==18 and all(checks.values())
 assert original_report['count']==18 and original_report['checks']==checks
 assert hashlib.sha256(executed.read_bytes()).hexdigest()==ORIGINAL_SHA
 for name in NAMES:assert hashlib.sha256((O/name).read_bytes()).hexdigest()==input_info[name]['sha256']
 status='PASS_ORIGINAL_READER_BYTE_EXACT_FULL_RATIONAL_EXECUTION'
except BaseException as exc:
 sys.settrace(None);error=''.join(traceback.format_exception(exc));print(error,flush=True)
 if isinstance(exc,AssertionError):status='ASSERTION_FAILURE_REVIEW_REQUIRED'
finally:
 finished=True;peak_rss=max(peak_rss,rss())
 report=emit('original_execution_finished',status,{'original_source_bytes_unchanged':True,'original_18_exact_checks_executed':status.startswith('PASS'),'original_final_assertion_executed':status.startswith('PASS'),'original_exact_check_count':18,'original_assert_statement_count':1,'execution_text_changed':False,'source_output':original_report,'error':error,'resource_control':'Original Python/native alarms disabled; original RLIMIT_AS delegated to outer RSS monitor','rss_limit_mib':args.rss_limit_mib,'series_truncation_used':False,'full_rational_function_arithmetic':True,'fresh_run_complete':status.startswith('PASS')})
 (O/'original_reader_run.json').write_text(json.dumps(report,indent=2),encoding='utf8')
if not status.startswith('PASS'):raise SystemExit(1)
