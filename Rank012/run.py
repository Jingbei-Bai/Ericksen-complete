"""Reconstruct and execute the calculations in this archive with SageMath."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, os, re, shutil, signal
import subprocess, sys, time, uuid
from datetime import datetime, timezone

HERE=Path(__file__).resolve().parent

def sha(path):
    with Path(path).open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()

def write(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.new')
    temp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    temp.replace(path)

def member(root,name):
    name=str(name).replace('\\','/')
    if Path(name).is_absolute() or re.match(r'^[A-Za-z]:',name):
        raise ValueError('Expected a relative package path: '+name)
    path=(root/name).resolve()
    if not path.is_relative_to(root.resolve()): raise ValueError('Path escapes the package: '+name)
    return path

def utc(): return datetime.now(timezone.utc).isoformat()

def work_tree_hashes(root):
    return {p.relative_to(root).as_posix():sha(p) for p in sorted(root.rglob('*'))
            if p.is_file() and p.relative_to(root).parts[0]!='run'}

def save_checkpoint(root,status):
    path=root/'run/work_tree.json'
    write(path,{'files':work_tree_hashes(root)})
    status['work_tree_checkpoint_sha256']=sha(path)

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result

def prepare(root,manifest,plan):
    if root.exists(): raise ValueError('Choose a new work directory, or use --resume for this run.')
    root.mkdir(parents=True)
    for row in manifest['files']:
        source=member(HERE,row['archive_path'])
        if not source.is_file() or sha(source)!=row['sha256']:
            raise ValueError('Package file differs from its manifest: '+row['archive_path'])
        target=member(root,row['work_path'])
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target)
    write(root/'execution_jobs.json',plan)
    write(root/'MANIFEST.json',{'files':[dict(path=r['work_path'],sha256=r['sha256'],kind=r['kind']) for r in manifest['files']]})

def check_receipt(root,job,result):
    """Check full source execution and the required completed coefficient rows."""
    wrapper=job.get('wrapper_member')
    stem=Path(job['source_member']).stem
    if wrapper:
        out=member(root,wrapper).parent
        extra=job.get('arguments',[])[1:]
        if job.get('execution_receipt'):
            path=member(root,job['execution_receipt'])
        elif wrapper.endswith(('rank3_fibres/run_one.py','rank3_ordinary/run_one.py')):
            path=out/'job_reports'/(stem+('_'+extra[0] if extra else '')+'.json')
        elif wrapper.endswith('isometries/replay_one.py'):
            path=out/'replay'/Path(job['arguments'][0]).parent/(stem+'_'+('_'.join(extra) if extra else 'default')+'_fresh_run.json')
        elif wrapper.endswith('rank3/replay_source_reader.py'):
            if stem=='independent_subcritical_single_RR_yw_minus_closed_reader' and not extra: extra=['H2']
            path=out/'source_replays'/(stem+('_'+extra[0] if extra else ''))/'replay_record.json'
        elif wrapper.endswith('rank3_coverage/run_generator.py'):
            path=out/(stem+('_'+'_'.join(extra) if extra else '')+'_run.json')
        else: raise ValueError('Unknown source execution wrapper: '+wrapper)
        data=json.loads(path.read_text())
        executed=data.get('executed_source_sha256',data.get('executed_sha256'))
        if data.get('status')!='PASS' or executed!=result['source_sha256'] or data.get('source_text_transformations',[]):
            raise ValueError('The full source execution receipt did not match '+job['id'])
        result['execution_text_receipt']={'verified':True,'path':str(path),'sha256':sha(path),'executed_source_sha256':executed}
    else:
        result['execution_text_receipt']={'verified':True,'method':'Direct source execution'}
    prefix='research/ericksen_20261007/full_recompute_20261009/'
    cover=member(root,prefix+'isometries/rank3_coverage/replay/round10')
    expectations={'independent_delta_subcritical_H_initials':('rows',18),
        'independent_ratio_p_boundary_global_pair_syzygies':('records',3),
        'independent_subcritical_two_RR_inactive_pair_generic_boundary':('rows',31)}
    if stem in expectations:
        path=cover/(stem+'.json');data=json.loads(path.read_text());key,count=expectations[stem]
        actual=len(data[key]) if isinstance(data[key],list) else data[key]
        if actual!=count: raise ValueError('Incomplete coefficient rows: '+stem)
        result['required_row_count']=count
    if stem=='independent_subcritical_two_RR_inactive_pair_free_boundary':
        if re.findall(r'^SAVED\s+(\d+)\s',Path(result['stdout_log']).read_text(errors='replace'),re.M)!=['31']:
            raise ValueError('Expected all 31 cyclic boundary rows')
        result['required_row_count']=31
    if stem=='independent_ratio_incoming12_all_third_range':
        raw=member(root,prefix+'rank3_ordinary/rebuild/round10/independent_ratio_raw_Q0_normal_two.sobj')
        data=json.loads(Path(result['execution_text_receipt']['path']).read_text())
        matching=[r for r in data.get('imports',[]) if Path(r['path']).resolve()==raw.resolve()]
        if not matching or not all(r['sha256']==sha(raw) and r.get('regenerated') is True for r in matching):
            raise ValueError('The complete incoming equations require freshly generated Q0 data')
        mode='raw' if job.get('environment',{}).get('INCOMING_KEEP_CUBIC')=='1' else 'range'
        report=json.loads((raw.parent/('independent_ratio_incoming12_all_third_'+mode+'.json')).read_text())
        if report['rows']!=37 or set(report['constant_support'])!={str(i) for i in range(37)}:
            raise ValueError('Expected all 37 incoming equation rows')
        result['required_row_count']=37

def execute(root,jobs,env,args,plan_hash,manifest_hash):
    records=root/'run';ledger=records/'results.json';records.mkdir(exist_ok=True)
    if args.resume:
        status=json.loads(ledger.read_text())
        if status['plan_sha256']!=plan_hash or status['manifest_sha256']!=manifest_hash:
            raise ValueError('Resume requires the same sources and task selection')
        checkpoint=records/'work_tree.json'
        if status['status']=='RUNNING' or not checkpoint.is_file() or sha(checkpoint)!=status.get('work_tree_checkpoint_sha256'):
            raise ValueError('Resume requires a completed or cleanly interrupted run; choose a new work directory after an abrupt termination.')
        if json.loads(checkpoint.read_text())['files']!=work_tree_hashes(root):
            raise ValueError('Generated outputs or execution receipts changed since the checkpoint')
        for result in status['results'].values():
            if result.get('supervisor_pid') and result['status']=='RUNNING':
                process=Path('/proc')/str(result['supervisor_pid'])/'cmdline'
                if process.exists() and b'external_resource_supervisor.py' in process.read_bytes():
                    raise ValueError('A calculation from this run is still running')
            if result['status']=='EXECUTION_PASS':
                if sha(result['stdout_log'])!=result['stdout_log_sha256']:
                    raise ValueError('An earlier execution log changed')
            else: result['retry_requested']=True
        for row in json.loads((HERE/'MANIFEST.json').read_text())['files']:
            if sha(member(root,row['work_path']))!=row['sha256']:
                raise ValueError('A declared source or input changed in the work directory')
    else:
        status={'run_id':str(uuid.uuid4()),'status':'RUNNING','started_utc':utc(),
            'plan_sha256':plan_hash,'manifest_sha256':manifest_hash,'registered_jobs':len(jobs),
            'inherited_execution_results_used':False,'generated_outputs_initially_absent':True,'results':{}}
    runner=load_module('calculation_runner',root/'run_calculations.py')
    active={};stop=False
    def request_stop(signum,frame):
        nonlocal stop
        stop=True
    signal.signal(signal.SIGINT,request_stop);signal.signal(signal.SIGTERM,request_stop)
    while True:
        for ident,entry in list(active.items()):
            code=entry['process'].poll()
            if code is None: continue
            entry['log'].close();result=status['results'][ident]
            supervision=json.loads(Path(result['supervisor_record']).read_text())
            result.update(status='EXECUTION_PASS' if code==0 and supervision.get('child_exit_code')==0 else 'EXECUTION_FAILURE',
                exit_code=code,finished_utc=utc(),seconds=time.monotonic()-entry['start'],
                peak_process_group_rss_bytes=supervision['peak_process_group_rss_bytes'],
                stdout_log_sha256=sha(result['stdout_log']))
            if supervision['status'].startswith('RESOURCE_ABORT'): result['status']='RESOURCE_ABORT'
            if result['status']=='EXECUTION_PASS':
                try:
                    if sha(member(root,result['source_member']))!=result['source_sha256']:
                        raise ValueError('Source changed during execution')
                    check_receipt(root,entry['job'],result)
                except Exception as error: result.update(status='EXECUTION_FAILURE',error=str(error))
            print(result['status'],ident,round(result['seconds'],2),'s',flush=True)
            del active[ident]
        if stop:
            for entry in active.values():
                report=json.loads(Path(entry['record']).read_text())
                try: os.killpg(report['leader']['pid'],signal.SIGTERM)
                except ProcessLookupError: pass
                entry['process'].wait();entry['log'].close()
                status['results'][entry['job']['id']]['status']='INTERRUPTED'
            status.update(status='INTERRUPTED',finished_utc=utc());save_checkpoint(root,status);write(ledger,status);return 130
        passed={k for k,v in status['results'].items() if v['status']=='EXECUTION_PASS'}
        attempted={k for k,v in status['results'].items() if not v.get('retry_requested',False)}
        namespaces={entry['namespace'] for entry in active.values()}
        ready=[]
        for job in jobs:
            namespace=job.get('wrapper_member') or job['id']
            if job['id'] in attempted or not set(job.get('depends_on',[]))<=passed or namespace in namespaces: continue
            if len(ready)>=max(0,args.workers-len(active)): break
            ready.append(job);namespaces.add(namespace)
        for job in ready:
            ident=job['id'];number=jobs.index(job)+1
            directory=records/'jobs'/f'{number:03d}'
            prior=status['results'].get(ident);attempt=1 if prior is None else prior.get('attempt',1)+1
            directory=directory/f'attempt{attempt}';directory.mkdir(parents=True)
            result={'id':ident,'run_id':status['run_id'],'status':'RUNNING','started_utc':utc(),
                'attempt':attempt,'stdout_log':str(directory/'stdout.log')}
            if job.get('action','run')=='copy_generated':
                src,dst=member(root,job['from']),member(root,job['to']);dst.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(src,dst);Path(result['stdout_log']).write_text('Copied freshly generated predecessor output.\n')
                result.update(status='EXECUTION_PASS',exit_code=0,finished_utc=utc(),seconds=0,
                    stdout_log_sha256=sha(result['stdout_log']),output_sha256=sha(dst))
                status['results'][ident]=result;continue
            command=runner.command(job,root)
            substitutions={'@WORKROOT@':str(root),'@LEDGER@':str(ledger),'@RUNROOT@':str(records)}
            for token,value in substitutions.items(): command=[a.replace(token,value) for a in command]
            source=member(root,job['source_member'])
            result.update(source_member=job['source_member'],source_sha256=sha(source),
                command=command,supervisor_record=str(directory/'supervision.json'))
            log=(directory/'supervisor.log').open('wb')
            supervisor=[sys.executable,str(root/'external_resource_supervisor.py'),
                '--record',result['supervisor_record'],'--log',result['stdout_log'],
                '--rss-limit-mib',str(args.rss_limit_mib),'--min-available-mib','2048','--',*command]
            job_environment=job.get('environment',{})
            if not isinstance(job_environment,dict) or any(not isinstance(k,str) or not isinstance(v,str) for k,v in job_environment.items()):
                raise ValueError('Task environment must map variable names to strings')
            result['task_environment']=job_environment
            proc=subprocess.Popen(supervisor,env=dict(env,**job_environment),cwd=root,stdout=log,stderr=subprocess.STDOUT)
            result['supervisor_pid']=proc.pid
            active[ident]={'job':job,'process':proc,'log':log,'start':time.monotonic(),
                'namespace':job.get('wrapper_member') or ident,'record':result['supervisor_record']}
            status['results'][ident]=result
            print('START',number,'/',len(jobs),ident,flush=True)
        status['completed_actual_passes']=sum(r['status']=='EXECUTION_PASS' for r in status['results'].values())
        status['last_updated_utc']=utc();write(ledger,status)
        if not active and not ready:
            status.update(status='PASS' if len(passed)==len(jobs) else 'FAIL',finished_utc=utc(),
                blocked_jobs=[j['id'] for j in jobs if j['id'] not in status['results']])
            save_checkpoint(root,status)
            write(ledger,status)
            print(status['status'],status['completed_actual_passes'],'/',len(jobs),flush=True)
            return 0 if status['status']=='PASS' else 1
        time.sleep(.5)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir',type=Path,default=Path('work'))
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--rss-limit-mib',type=int,default=8192)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--list',action='store_true')
    parser.add_argument('--task',action='append',default=[])
    parser.add_argument('--proof',action='append',default=[])
    parser.add_argument('--list-proofs',action='store_true')
    args=parser.parse_args()
    proofmap=json.loads((HERE/'PROOF_MAP.json').read_text())
    if args.list_proofs:
        for n in proofmap['nodes']:
            print(n['node']+' : '+', '.join(n['paper_labels']))
        return 0
    if args.proof:
        selected={n['node']:n for n in proofmap['nodes']}
        for ident in args.proof:
            if ident not in selected: parser.error('Unknown mathematical part: '+ident)
            args.task += [d['task_id'] for d in selected[ident]['direct_tasks']]
        args.task=list(dict.fromkeys(args.task))
        if not args.task: parser.error('Selected parts have analytic proofs and no symbolic calculation')
    binder=load_module('current_paper_index',HERE/'paper_index.py')
    binder.validate(HERE)
    if args.workers<1: parser.error('--workers must be positive')
    runner=load_module('package_runner',HERE/'runtime/run_calculations.py')
    jobs=runner.topological_jobs(runner.read_jobs(HERE/'tasks.json'),[],args.task)
    if args.list:
        for i,j in enumerate(jobs,1): print(f'{i:03d} {j["id"]} <- '+', '.join(j.get('depends_on',[])))
        return 0
    if sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE'): raise ValueError('Run without Python optimization; calculations use assertions.')
    if not Path('/proc/meminfo').exists(): raise ValueError('Use the fixed Linux or WSL SageMath environment.')
    root=args.work_dir.resolve();manifest=json.loads((HERE/'MANIFEST.json').read_text())
    plan={'jobs':jobs};plan_hash=hashlib.sha256(json.dumps(plan,sort_keys=True).encode()).hexdigest()
    manifest_hash=sha(HERE/'MANIFEST.json')
    if not args.resume: prepare(root,manifest,plan)
    environment=root/'run/environment.json'
    subprocess.run([sys.executable,str(HERE/'check_environment.py'),'--output',str(environment)],check=True)
    env=dict(os.environ,ERICKSEN_BUNDLE_ROOT=str(root),ERICKSEN_FORCE_RECOMPUTE='1',
        PYTHONPATH=str(root/'runtime_support')+os.pathsep+os.environ.get('PYTHONPATH',''),
        PYTHONUNBUFFERED='1',PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    return execute(root,jobs,env,args,plan_hash,manifest_hash)

if __name__=='__main__':
    try: raise SystemExit(main())
    except (ValueError,OSError,subprocess.CalledProcessError) as error:
        print('ERROR:',error,file=sys.stderr);raise SystemExit(1)
