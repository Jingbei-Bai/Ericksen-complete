"""Independent full-source/degree/ZZ-node audit of two exact Cramer certificates."""
from pathlib import Path
import sys,time,json,hashlib
from sage.all import *
out=Path(__file__).resolve().parent;kind=sys.argv[1];assert kind in ('ranktwo','zero_det_B')
prefix='general_curve_'+kind;start=time.monotonic();checks={}
inp=load(str(out/(prefix+'_interpolation_input.sobj')))
mac=load(str(out/inp['input_file']))
rep=json.loads((out/(prefix+'_integer_interpolation.json')).read_text())
deg=json.loads((out/(prefix+'_interpolation_degree_certificate.json')).read_text())
QH=PolynomialRing(QQ,'readerH');hh=QH.gen();F=QH.fraction_field()
ZH=PolynomialRing(ZZ,'readerH');CH=PolynomialRing(QQ,names=('readerH','readerZ','readerS'));h,z,s=CH.gens()
def cp(p):
    if hasattr(p,'dict'):return CH({tuple(map(int,e)):QQ(c) for e,c in p.dict().items()})
    if hasattr(p.numerator(),'dict'):return CH(cp(p.numerator())/cp(p.denominator()))
    return CH(QQ(p))
def uf(p):
    nn=p.numerator();dd=p.denominator()
    return F(QH(nn.list() if hasattr(nn,'list') else nn))/QH(dd.list() if hasattr(dd,'list') else dd)
def uq(p):return QH(list(p))
def uz(p):return ZH(list(p))
def layers(p):
    ans={}
    for e,c in p.dict().items():
        eh,ez,es=map(int,e);key=(ez,es);ans[key]=ans.get(key,QH.zero())+QQ(c)*hh**eh
    return ans
def twodict(p):return {tuple(map(int,e)):uf(c) for e,c in p.dict().items()}
if kind=='ranktwo':
    src=load(str(out/'general_curve_canonical_ranktwo_linear.sobj'))
    original=[cp(f) for f in src['parameter_sources']]
    target=(h-z)**2*(z-1)**2*s**2*(s-1)**2
else:
    src=load(str(out/'general_curve_zero_det_coefficients_input.sobj'))
    canon=load(str(out/src['source_file']))
    checks['canonical_primitive_source']=src['primitive']==canon['primitive']
    checks['canonical_clear_content']=src['determinant']*src['denominator']==src['primitive']*src['content']
    guard=h*z*(h-z)*(z-1)*(h+1)*(1+h-z)*(s-1)*(h*s-z+1)*(h*s-z+s)*s*(h*s-z)
    for i,it in enumerate(src['items']):
        checks[f'primitive_coefficient_{i}']=cp(it['coefficient'])==cp(src['primitive'][i])
        checks[f'coefficient_stripping_{i}']=cp(it['coefficient'])==cp(it['protected_factor'])*cp(it['source'])
        checks[f'only_guard_stripped_{i}']=all(not guard%f for f,e in cp(it['protected_factor']).factor())
    original=[cp(f) for f in src['sources']]
    B=h*z*(h*h-1)*s*s+2*(z*z-h*h)*s-h*z*(z*z-1)
    target=s*(s-1)*B
source_layers=[layers(p) for p in original];target_layer=layers(target)
checks['all_original_necessary_sources']=all(twodict(p)=={e:F(v) for e,v in layer.items()} for p,layer in zip(mac['sources'],source_layers)) and len(mac['sources'])==len(original)
checks['target_exact']=twodict(mac['target'])=={e:F(v) for e,v in target_layer.items()}
mons=[tuple(map(int,e)) for e in mac['monomials']]
cols=[]
for sid,shift in mac['provenance']:
    u,v=map(int,shift);assert u>=0 and v>=0
    cols.append({(e[0]+u,e[1]+v):F(c) for e,c in source_layers[int(sid)].items()})
checks['complete_coefficient_support']=set(target_layer).union(*(set(c) for c in cols))<=set(mons) and len(set(mons))==len(mons)
A=[[col.get(e,F.zero()) for col in cols] for e in mons];b=[F(target_layer.get(e,QH.zero())) for e in mons]
checks['all_original_macaulay_entries']=all(A[i][j]==uf(mac['matrix'][i,j]) for i in range(len(A)) for j in range(len(cols)))
checks['all_original_target_entries']=b==[uf(v) for v in mac['target_vector']]
n=len(cols);nrows=len(A);selected=list(map(int,mac['selected_rows']))
assert len(selected)==n and len(set(selected))==n
Bmat=[[uz(v) for v in row] for row in inp['square_matrix'].rows()];brhs=[uz(v) for v in inp['square_rhs']]
ALL=[[uz(v) for v in row] for row in inp['all_rows'].rows()];rhs=[uz(v) for v in inp['all_rhs']]
contents=[uf(v) for v in inp['all_row_contents']];qscales=[QQ(v) for v in inp['all_row_rational_scales']]
column_contents=[uf(v) for v in inp['column_contents']]
checks['nonzero_normalizing_factors']=all(contents) and all(qscales) and all(column_contents)
checks['only_constant_column_contents']=all(v.numerator().degree()==0 and v.denominator().degree()==0 for v in column_contents)
checks['all_original_row_normalizations']=all(F(ALL[i][j])*contents[i]*column_contents[j]==qscales[i]*A[i][j] for i in range(nrows) for j in range(n))
checks['all_original_rhs_normalizations']=all(F(rhs[i])*contents[i]==qscales[i]*b[i] for i in range(nrows))
if kind=='ranktwo':
    pre=load(str(out/'general_curve_ranktwo_base_flint_preconditioned.sobj'))
    rscale=[uf(v) for v in pre['row_contents']];qs=[QQ(v) for v in pre['rational_scales']]
    checks['square_matrix_source_normalization']=all(F(Bmat[i][j])*rscale[i]*column_contents[j]==qs[i]*A[selected[i]][j] for i in range(n) for j in range(n))
    checks['square_rhs_source_normalization']=all(F(brhs[i])*rscale[i]==qs[i]*b[selected[i]] for i in range(n))
else:
    checks['square_matrix_source_normalization']=Bmat==[ALL[i] for i in selected]
    checks['square_rhs_source_normalization']=brhs==[rhs[i] for i in selected]
u=list(map(int,deg['selected_row_potentials']));v=list(map(int,deg['column_potentials']));ul=int(deg['remaining_row_potential']);degree=int(deg['degree_bound'])
checks['degree_potential_sum']=degree==sum(u)+ul+sum(v)
checks['EVERY_square_augmented_degree']=all(not f or f.degree()<=u[i]+v[j] for i,row in enumerate(Bmat) for j,f in enumerate(row+[brhs[i]]))
checks['EVERY_remaining_augmented_degree']=all(not f or f.degree()<=ul+v[j] for i,row in enumerate(ALL) for j,f in enumerate(row+[rhs[i]]))
checks['enough_distinct_nodes']=len(set(rep['nodes']))==len(rep['nodes'])==degree+1 and int(rep['completed_nodes'])==degree+1
assert all(checks.values()),[k for k,v in checks.items() if not v]
print(kind,'ALL source/normalization/degree checks passed',time.monotonic()-start,flush=True)
FF=GF(32003);nodes=[];blockchecks=[];blocks=out/(prefix+'_integer_nodes')
for blk in rep['blocks']:
    path=blocks/blk['file'];assert hashlib.sha256(path.read_bytes()).hexdigest()==blk['sha256']
    data=load(str(path));assert len(data)==int(blk['count'])
    for nd in data:
        hv=ZZ(nd['H']);nums=vector(ZZ,nd['numerators']);den=ZZ(nd['denominator']);assert den and len(nums)==n
        square=matrix(ZZ,[[p(hv) for p in row] for row in Bmat]);bb=vector(ZZ,[p(hv) for p in brhs])
        md=square.change_ring(FF).det();assert md and int(md)==int(nd['det_mod'])
        assert not any(square*nums-den*bb)
        allmat=matrix(ZZ,[[p(hv) for p in row] for row in ALL]);rb=vector(ZZ,[p(hv) for p in rhs])
        assert not any(allmat*nums-den*rb)
        nodes.append(int(hv))
    blockchecks.append(blk['file'])
    progress={'status':'EXACT_INTEGER_WITNESSES_REPLAYING','kind':kind,'nodes':len(nodes),'required':degree+1,'seconds':time.monotonic()-start}
    (out/('independent_'+kind+'_interpolation_progress.json')).write_text(json.dumps(progress));print(progress,flush=True)
checks['all_blocks_and_nodes_match_manifest']=nodes==list(map(int,rep['nodes'])) and len(nodes)==degree+1
checks['ALL_nodes_square_nonsingular_and_ALL_original_integer_rows']=True
assert all(checks.values())
scope='Exact all-parameter Cramer source identity, including every original coefficient row, all source-preserving normalizations, rigorous augmented degree bound, and all stored characteristic-zero integer witnesses. New H denominators have only finitely many zeros; no pointwise exclusion there is claimed.'
result={'checks_passed':len(checks),'checks':checks,'nodes':len(nodes),'degree_bound':degree,'matrix_shape':[nrows,n],'seconds':time.monotonic()-start,'scope':scope}
save({'checks':checks,'source_layers':source_layers,'target_layer':target_layer,'selected_rows':selected,'verified_nodes':nodes,'verified_blocks':blockchecks,'scope':scope},str(out/('independent_'+kind+'_interpolation_reader.sobj')))
(out/('independent_'+kind+'_interpolation_reader.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

