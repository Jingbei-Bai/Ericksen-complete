"""External process-group RSS supervision; never shares the CAS GIL.

Only the independently created child process group or an explicitly identified
group whose leader PID equals --pid is supervised. A resource stop is recorded
from actual measurements and never classified as a mathematical assertion fail.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,os,signal,subprocess,sys,time
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--pid',type=int)
p.add_argument('--record',type=Path,required=True)
p.add_argument('--log',type=Path)
p.add_argument('--rss-limit-mib',type=int,required=True)
p.add_argument('--min-available-mib',type=int,default=2048)
p.add_argument('command',nargs=argparse.REMAINDER)
a=p.parse_args();a.record.parent.mkdir(parents=True,exist_ok=True)
command=a.command[1:] if a.command and a.command[0]=='--' else a.command
assert bool(a.pid)!=bool(command),'Use either --pid or a child command'
def info(pid):
 st=Path('/proc')/str(pid)
 text=(st/'stat').read_text();parts=text[text.rfind(')')+2:].split()
 pages=int((st/'statm').read_text().split()[1]);cmd=(st/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
 return {'pid':pid,'pgrp':int(parts[2]),'start_ticks':int(parts[19]),'rss_bytes':pages*os.sysconf('SC_PAGE_SIZE'),'command':cmd}
def available():
 for line in Path('/proc/meminfo').read_text().splitlines():
  if line.startswith('MemAvailable:'):return int(line.split()[1])*1024
 raise RuntimeError('MemAvailable missing')
def group_rss(group):
 members=[]
 for directory in Path('/proc').iterdir():
  if directory.name.isdigit():
   try:
    item=info(int(directory.name))
    if item['pgrp']==group:members.append(item)
   except (OSError,ValueError,IndexError):pass
 return sum(item['rss_bytes'] for item in members),members
child=None;loghandle=None
if command:
 assert a.log is not None,'A raw stdout/stderr log is required for new execution'
 a.log.parent.mkdir(parents=True,exist_ok=True);loghandle=a.log.open('wb')
 child=subprocess.Popen(command,start_new_session=True,stdout=loghandle,stderr=subprocess.STDOUT)
 pid=child.pid
else:pid=a.pid
first=info(pid);group=first['pgrp']
assert group==pid and group!=os.getpgrp(),first
st=time.monotonic();peak=0;minimum=available();cause=None;last=None
record={'schema':1,'status':'RUNNING','supervisor_pid':os.getpid(),'supervisor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'leader':first,'command':command or first['command'],'raw_stdout_log':str(a.log) if a.log else None,'rss_limit_bytes':a.rss_limit_mib*1024**2,'minimum_memavailable_bytes':a.min_available_mib*1024**2,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
def save():
 record.update({'seconds_supervised':time.monotonic()-st,'peak_process_group_rss_bytes':peak,'lowest_memavailable_bytes':minimum,'last_measurement':last})
 tmp=a.record.with_suffix(a.record.suffix+'.tmp');tmp.write_text(json.dumps(record,indent=2),encoding='utf8');tmp.replace(a.record)
while True:
 if child is not None and child.poll() is not None:break
 try:
  current=info(pid)
  if current['start_ticks']!=first['start_ticks']:break
 except (OSError,ValueError,IndexError):break
 total,members=group_rss(group);mem=available();peak=max(peak,total);minimum=min(minimum,mem)
 last={'process_group_rss_bytes':total,'memavailable_bytes':mem,'members':members,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 if total>a.rss_limit_mib*1024**2:cause='PROCESS_GROUP_RSS_EXCEEDED_ACTUAL_BUDGET'
 elif mem<a.min_available_mib*1024**2:cause='MEMAVAILABLE_BELOW_ACTUAL_RESERVE'
 if cause:
  record['status']='RESOURCE_ABORT_NOT_MATHEMATICAL_FAILURE';record['actual_abort_reason']=cause;save()
  os.killpg(group,signal.SIGTERM);record['termination_signals']=['SIGTERM']
  deadline=time.monotonic()+3
  while time.monotonic()<deadline:
   if child is not None and child.poll() is not None:break
   try:info(pid)
   except OSError:break
   time.sleep(.1)
  try:
   info(pid);os.killpg(group,signal.SIGKILL);record['termination_signals'].append('SIGKILL')
  except OSError:pass
  break
 save();time.sleep(.5)
if child is not None:
 code=child.wait();record['child_exit_code']=code
else:
 record['child_exit_code']=None
if cause is None:record['status']='ACTUAL_CHILD_EXECUTION_PASS' if child is not None and code==0 else 'PROCESS_FINISHED_WITHOUT_RESOURCE_STOP'
record['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
if a.log and a.log.exists():record['raw_stdout_log_sha256']=hashlib.sha256(a.log.read_bytes()).hexdigest()
save()
if loghandle:loghandle.close()
print(json.dumps(record,indent=2),flush=True)
if cause or (child is not None and code):raise SystemExit(1)
