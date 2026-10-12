"""Match printed N+- formulas to corrected, freshly sourced Laurent norms.

No archived third row, stream, or infinity-probe output is read here.
"""
from pathlib import Path
from sage.all import *
import json,hashlib
h=Path(__file__).absolute().parent
ordinary=h.parent.parent/'rank3_ordinary'
f=ordinary/'rebuild/round10/independent_incoming_laurent_source_reader.sobj'
report=ordinary/'job_reports/independent_incoming_laurent_source_reader.json'
assert json.loads(report.read_text())['status']=='PASS'
d=load(str(f));checks=[]
def ck(n,v):assert v,n;checks.append(n)
for item in d['results']:
 K=item['leading_coefficient'].parent();A,B=K.gens()
 cm=2352*A**3*B**2+961*A**2*B**3-672*A**3*B-3970*A**2*B**2-1984*A*B**3+48*A**3+577*A**2*B+1664*A*B**2+1024*B**3
 cp=1200*A**3*B**2+1089*A**2*B**3+480*A**3*B-2946*A**2*B**2-2112*A*B**3+48*A**3-575*A**2*B+1792*A*B**2+1024*B**3
 nm=QQ(9216)/49*A*A*(B-1)**4*(B-A)**14/((A-1)**18*B**17)*cm*(7*A*B-A-6*B)**2*(11*A*B-3*A-8*B)**4/((5*A*B-A-4*B)**6*(2*A*B-A-B)**2)
 np=QQ(9216)/49*A*A*(B-1)**4*(B-A)**12/((A-1)**18*B**17)*cp*(5*A*B+A-6*B)**2*(5*A*B+3*A-8*B)**4/(3*A*B+A-4*B)**6
 sign=int(item['ray_sign']);cycle=int(item['cycle']);suffix=str((sign,cycle))
 ck('fresh corrected source norm equals literal printed N'+suffix,item['leading_coefficient']==(nm if sign<0 else np))
 ck('literal printed Laurent exponent and absolute remainder '+suffix,item['norm']['valuation_lower_bound']==-30 and item['norm']['precision_absolute']>=-29)
 ck('all printed leading coefficients nonzero in symbolic material field '+suffix,bool(item['leading_coefficient']))
result=dict(status='PASS',checks=checks,fresh_source=str(f),fresh_source_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),fresh_job_report_sha256=hashlib.sha256(report.read_bytes()).hexdigest(),latex_equations=['appr3tau:negative','appr3tau:positive','appr3tau:laurent'],scope='Literal printed rational coefficients and Laurent orders agree with the corrected freshly reconstructed source norm. Entire material factor/inverse audit remains the178-check fresh Laurent reader. This checker does not load archived third rows or archived infinity probes.')
(h/'nullwall_printed_norms.json').write_text(json.dumps(result,indent=2))
print('FRESH NULLWALL PRINTED NORMS PASS',len(checks),flush=True)
