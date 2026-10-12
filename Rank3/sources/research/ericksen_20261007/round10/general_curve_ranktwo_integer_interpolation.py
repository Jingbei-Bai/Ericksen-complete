"""Deterministic characteristic-zero identity certificate by exact interpolation.

There are D+1 distinct integer nodes, with a checked integer degree bound D.
Each square matrix is nonsingular (nonzero modular determinant), its solution
is computed over QQ, and ALL original normalized source rows are checked over
ZZ with a common integer denominator. This is not modular interpolation.
"""
from sage.all import *
from pathlib import Path
import json,time,hashlib
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_ranktwo_interpolation_input.sobj'));degree=json.loads((out/'general_curve_ranktwo_interpolation_degree_certificate.json').read_text());D=int(degree['degree_bound']);need=D+1;B=d['square_matrix'];c=d['square_rhs'];A=d['all_rows'];b=d['all_rhs'];ff=GF(32003);blocks=out/'general_curve_ranktwo_integer_nodes';blocks.mkdir(exist_ok=True)
assert B.nrows()==B.ncols()==185 and A.nrows()==276
# Recheck the elementary integer degree inequalities before any evaluations.
u=degree['selected_row_potentials'];v=degree['column_potentials'];ul=degree['remaining_row_potential'];assert D==sum(u)+ul+sum(v)
for i in range(B.nrows()):
    for j,f in enumerate(list(B.row(i))+[c[i]]):assert not f or f.degree()<=u[i]+v[j]
for i in range(A.nrows()):
    for j,f in enumerate(list(A.row(i))+[b[i]]):assert not f or f.degree()<=ul+v[j]
nodes=[];skipped=[];block=[];manifest=[];h=2;print('EXACT_NODE_LOOP',need,'degree',D,flush=True)
def checkpoint(final=False):
    global block
    if block:
        name='block_'+str(len(manifest)).zfill(4)+'.sobj';tmp=blocks/(name+'.next.sobj');dest=blocks/name;save(block,str(tmp));tmp.replace(dest);manifest.append({'file':name,'count':len(block),'first_H':block[0]['H'],'last_H':block[-1]['H'],'sha256':hashlib.sha256(dest.read_bytes()).hexdigest()});block=[]
    rep={'status':'EXACT_ALL_PARAMETER_IDENTITY' if final else 'EXACT_NODES_RUNNING','degree_bound':D,'required_nodes':need,'completed_nodes':len(nodes),'nodes':nodes,'skipped_modular_singular_nodes':skipped,'prime_for_nonsingularity':32003,'blocks':manifest,'seconds':time.monotonic()-start,'all_rows':276,'square_size':185,'scope':'For every augmented determinant, D+1 distinct exact-Q zero evaluations and the integer degree bound imply the polynomial is identically zero. Nonzero modular determinants prove nonsingularity only; every residual is checked in characteristic zero.'};tmp=out/'general_curve_ranktwo_integer_interpolation_progress.next.json';tmp.write_text(json.dumps(rep,indent=2)+'\n');tmp.replace(out/'general_curve_ranktwo_integer_interpolation_progress.json')
    if final:(out/'general_curve_ranktwo_integer_interpolation.json').write_text(json.dumps(rep,indent=2)+'\n')
    print('EXACT_NODES',len(nodes),time.monotonic()-start,flush=True)
while len(nodes)<need:
    M=matrix(ZZ,[[f(h) for f in row] for row in B.rows()]);detmod=M.change_ring(ff).det()
    if not detmod:skipped.append(h);h+=1;continue
    rhs=vector(ZZ,[f(h) for f in c]);x=M.solve_right(rhs);den=lcm([q.denominator() for q in x]);nums=vector(ZZ,[den*q for q in x]);assert den and M*nums==den*rhs
    AM=matrix(ZZ,[[f(h) for f in row] for row in A.rows()]);ab=vector(ZZ,[f(h) for f in b]);assert AM*nums==den*ab
    block.append({'H':h,'det_mod':int(detmod),'numerators':list(nums),'denominator':den});nodes.append(h);h+=1
    if len(block)==50:checkpoint()
assert len(set(nodes))==need;checkpoint(final=True)
