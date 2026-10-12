"""Construct Macaulay inputs from the direction matrices and normalize exact rows."""
from pathlib import Path
from sage.all import *
import argparse, hashlib, json, shutil

def sha(path):
    with path.open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('phase',choices=['sources','precondition'])
parser.add_argument('--bundle-root',type=Path,required=True)
parser.add_argument('--work-root',type=Path,required=True)
parser.add_argument('--ledger',type=Path,required=True)
args=parser.parse_args()
root=args.bundle_root.resolve();work=args.work_root.resolve();work.mkdir(parents=True,exist_ok=True)
original=root/'research/ericksen_20261007/round10'
fresh=root/'research/ericksen_20261007/full_recompute_20261009/rank3'
ledger=json.loads(args.ledger.read_text())
job=ledger['results']['rank3_sources:fresh_rank_matrix_certificates']
assert job['status']=='EXECUTION_PASS' and job['exit_code']==0 and job['run_id']==ledger['run_id']
assert job['source_sha256']==sha(fresh/'fresh_rank_matrix_certificates.py')
assert json.loads((fresh/'fresh_rank_matrix_certificates.json').read_text())['status']=='PASS'
checks={}
def ck(name,value):
    assert bool(value),name
    checks[name]=True
def convert(ring,p): return ring({tuple(map(int,e)):QQ(c) for e,c in p.dict().items()})

if args.phase=='sources':
    source=load(str(fresh/'fresh_rank_matrix_certificates.sobj'))
    reference=load(str(original/'general_curve_canonical_ranktwo_linear.sobj'))
    ring=reference['ring']
    parameters=[convert(ring,p) for p in source['lowrank_parameter_sources']]
    ck('28_coefficient_polynomials',len(parameters)==len(reference['parameter_sources'])==28)
    for i,(a,b) in enumerate(zip(parameters,reference['parameter_sources'])): ck('coefficient_'+str(i),a==b)
    save(dict(reference,parameter_sources=parameters),str(work/'general_curve_canonical_ranktwo_linear.sobj'))
    reference=load(str(original/'general_curve_zero_det_coefficients_input.sobj'))
    determinants=list(load(str(fresh/'fresh_FG_sources.sobj'))['D'])
    ring=reference['ring'];field=ring.fraction_field();items=[];sources=[]
    ck('three_determinant_coefficients',len(reference['items'])==len(determinants)==3)
    for i,(raw,item) in enumerate(zip(determinants,reference['items'])):
        coefficient=convert(ring,raw);unit=field(coefficient)/field(item['coefficient'])
        ck('coefficient_unit_'+str(i),unit in QQ and unit!=0)
        coefficient=ring(coefficient/QQ(unit))
        residual=ring(field(coefficient)/field(item['protected_factor']))
        ck('coefficient_normalization_'+str(i),coefficient==item['coefficient'])
        ck('protected_factor_identity_'+str(i),residual==item['source']==reference['sources'][i])
        items.append(dict(item,coefficient=coefficient,source=residual));sources.append(residual)
    save(dict(reference,items=items,sources=sources),str(work/'general_curve_zero_det_coefficients_input.sobj'))
    primitive=reference['source_file']
    assert Path(primitive).name==primitive
    shutil.copyfile(original/primitive,work/primitive)
elif args.phase=='precondition':
    data=load(str(work/'general_curve_ranktwo_base_macaulay_input.sobj'))
    matrix=data['matrix'];target=data['target_vector'];field=matrix.base_ring();rational_ring=field.ring()
    integer_ring=PolynomialRing(ZZ,'H');rows=data['selected_rows'];n=matrix.ncols()
    integer_rows=[];integer_rhs=[];scales=[];contents=[]
    assert len(rows)==n and len(set(rows))==n
    for i in rows:
        values=[rational_ring(c) for c in matrix.row(i)]+[rational_ring(target[i])]
        scale=lcm([c.denominator() for v in values for c in v.list()] or [ZZ.one()])
        normalized=[integer_ring(scale*v) for v in values];content=gcd(normalized)
        ck('nonzero_row_factors_'+str(i),scale!=0 and content!=0)
        normalized=[v//content for v in normalized]
        integer_rows.append(normalized[:-1]);integer_rhs.append(normalized[-1]);scales.append(scale);contents.append(content)
    column_contents=[]
    for j in range(n):
        content=gcd([integer_rows[i][j] for i in range(n)])
        ck('constant_column_unit_'+str(j),bool(content) and content.degree()==0)
        column_contents.append(content)
        for i in range(n): integer_rows[i][j]//=content
    ck('all_matrix_entries',all(field(integer_rows[j][k])*field(contents[j])*field(column_contents[k])==scales[j]*matrix[i,k] for j,i in enumerate(rows) for k in range(n)))
    ck('all_target_entries',all(field(integer_rhs[j])*field(contents[j])==scales[j]*target[i] for j,i in enumerate(rows)))
    save(dict(input_file='general_curve_ranktwo_base_macaulay_input.sobj',integer_matrix_rows=integer_rows,
        integer_rhs=integer_rhs,row_contents=contents,column_contents=column_contents,rational_scales=scales,
        source_identity=True),str(work/'general_curve_ranktwo_base_flint_preconditioned.sobj'))
(work/(args.phase+'_identity.json')).write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n')
print('PASS',args.phase,len(checks),'exact identities',flush=True)
