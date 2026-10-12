"""Run the complete dependency DAG against a newly reconstructed state source.

Usage: sage -python run_fresh_rebuild.py /path/to/new_generic_source.sobj
The staged scripts contain historical algorithms, not historical results.
Each dependent certificate is recomputed in this directory during this run.
"""
from pathlib import Path
import subprocess,sys,time,json,hashlib,pickle
from sage.all import load
HERE=Path(__file__).absolute().parent
source=Path(sys.argv[1]).absolute()
raw=source.read_bytes()
try: data=pickle.loads(raw)
except Exception:data=load(str(source))
if 'basis' not in data and 'G5' in data:data['basis']=data['G5']
for key in ['Fraw','C3','Q4','basis']:assert key in data,key
dest=HERE/'rebuild'/'round3'/'generic_degree5_sparse_prepared.sobj';dest.parent.mkdir(exist_ok=True,parents=True)
dest.write_bytes(pickle.dumps(data,protocol=4))
jobs=[]
def job(round_,name,*args):jobs.append((HERE/'rebuild'/round_/(name+'.py'),list(args)))
for name in '''homothety_first_jet_reduction homothety_first_jet_minor_atlas
homothety_second_jet_reduction homothety_invariant_linear_chart
homothety_fiber_leading_compression homothety_observation_boundary_rank
homothety_phase_fiber_reader homothety_onezero_fiber_leading
homothety_axis_infinity_leading homothety_axis_jet_bound_reader
homothety_axis_plane_normal_sources homothety_axis_plane_forcing_schur
homothety_axis_plane_quadratic_schur homothety_middle_null_cubic_face
homothety_middle_null_cubic_solution homothety_middle_null_cubic_unprojected
homothety_middle_null_quartic_face'''.split():job('round12',name)
for name in ['homothety_middle_null_endpoint_cubic','homothety_middle_endpoint_uniform_quartic']:
 job('round12',name);job('round12',name,'--pzero')
for name in '''homothety_middle_null_endpoint_reader homothety_middle_endpoint_whole_reader
homothety_middle_resonant_whole_coefficient_reader homothety_middle_resonant_independent_reader
homothety_observation_transport_reader'''.split():job('round12',name)
job('round23','dilation_pure_short_source_reader')
for mode in ['general','pzero','qzero','resonance']:job('round23','dilation_state_column_linear_module',mode)
for name in ['dilation_state_homothety_boundary_carrier','dilation_state_boundary_source_reader','dilation_global_axis_pressure_sources']:job('round23',name)
jobs.append((HERE/'carrier_printed_check.py',[]))
jobs.append((HERE/'boundary_printed_check.py',[]))
start=time.time();records=[]
logdir=HERE/'logs';logdir.mkdir(exist_ok=True)
for i,(script,args) in enumerate(jobs):
 tag=str(i+1).zfill(2)+'_'+script.stem+('_'+'_'.join(a.strip('-') for a in args) if args else '')
 print('START',tag,flush=True);t=time.time()
 proc=subprocess.run([sys.executable,str(script),*args],capture_output=True,text=True)
 logfile=logdir/(tag+'.log');logfile.write_text(proc.stdout+'\nSTDERR:\n'+proc.stderr,encoding='utf-8')
 record={'order':i+1,'script':str(script.relative_to(HERE)),'args':args,'exit_code':proc.returncode,
         'seconds':time.time()-t,'log':str(logfile.relative_to(HERE)),
         'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest()}
 records.append(record)
 report={'source_path':str(source),'source_sha256':hashlib.sha256(raw).hexdigest(),
         'copied_state_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
         'historical_result_files_used':False,'jobs':records,'seconds':time.time()-start,
         'status':'PASS' if i==len(jobs)-1 and proc.returncode==0 else ('RUNNING' if proc.returncode==0 else 'FAIL')}
 (HERE/'fresh_rebuild_run.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print('END',tag,'exit',proc.returncode,'seconds',round(record['seconds'],2),flush=True)
 if proc.returncode:
  print(proc.stdout[-3000:]+'\n'+proc.stderr[-3000:],flush=True)
  raise SystemExit(proc.returncode)
print('COMPLETE',len(records),'jobs',round(time.time()-start,2),'seconds',flush=True)
