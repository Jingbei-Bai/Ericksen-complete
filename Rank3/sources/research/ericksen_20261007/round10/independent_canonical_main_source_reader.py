"""Independent source-to-four-certificate replay for the main canonical chart."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
old=load(str(p/'independent_ratio_canonical_amplitude_normal.sobj'))
new=load(str(p/'independent_second_prolongation_reader.sobj'))
short=load(str(p/'independent_ratio_canonical_amplitude_short_rows.sobj'))
chart=load(str(p/'independent_ratio_canonical_amplitude_tau_chart.sobj'))
src=load(str(p/'independent_ratio_canonical_second_residual_degree_source.sobj'))
R=old['matrix'].base_ring();K=R.base_ring();a,r=K.gens();P=K.ring();q=1+r;m=1+a*r;t=R.gens();AF=R.fraction_field();D=t[2]*sum(t)
def cv(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
def guarded(f,extra=()):
    f=P(f)
    for g in [P(a),P(r),P(q),P(m),P(a*q+1)]+list(extra):
        while f and f.total_degree()>0:
            quo,rem=f.quo_rem(g)
            if rem:break
            f=quo
    return bool(f) and f.total_degree()==0
M=old['matrix']; sr=matrix(R,[[cv(v) for v in row] for row in short['rows'].rows()])
for ind,(i,j,c) in enumerate(short['pairs']):
    raw=M[j,c]*M.row(i)-M[i,c]*M.row(j)
    e=next(k for k in range(6) if sr[ind,k]); scalar=AF(raw[e])/AF(sr[ind,e]);ck('short_scale_material_'+str(ind),scalar in K)
    scalar=K(scalar);ck('short_scale_protected_'+str(ind),guarded(scalar.numerator()) and guarded(scalar.denominator()))
    for k in range(6):ck('original_short_row_%d_%d'%(ind,k),raw[k]==scalar*sr[ind,k])
selected=matrix(R,[M.row(0),sr.row(0),sr.row(1),M.row(1),M.row(3)])
minor=selected.matrix_from_columns([0,1,2,4,5]).det()
for f,pow in minor.factor():
    if f.is_constant():
        c=K(f);ck('non_tau_material_minor_'+str(len(checks)),guarded(c.numerator()) and guarded(c.denominator()))
    else:ck('non_tau_amplitude_minor_'+str(len(checks)),f/f.lc() in (t[2],sum(t)))
ck('non_tau_minor_nonzero',minor!=0)
lift=vector(AF,[AF(R(f.numerator()))/AF(R(f.denominator())) for f in chart['normal_lift_for_tau_one']])
ck('tau_lift_coordinate',lift[3]==1)
for i,v in enumerate(selected.change_ring(AF)*lift):ck('five_lift_rows_'+str(i),v==0)
for j,f in enumerate(lift):
    for fac,pow in f.denominator().factor():
        if fac.is_constant():
            c=K(fac);ck('lift_material_'+str(len(checks)),guarded(c.numerator()) and guarded(c.denominator()))
        else:ck('lift_amplitude_'+str(len(checks)),fac/fac.lc() in (t[2],sum(t)))
allrows=list(M.rows())+[[cv(v) for v in row] for row in new['normal_matrix'].rows()]
polys=[];indices=[]
for i,row in enumerate(allrows):
    f=(vector(AF,row)*lift)*D;ck('full_residual_polynomial_'+str(i),f.denominator() in K)
    if f:polys.append(R(f));indices.append(i)
ck('source_count_84',len(polys)==84)
for j,f in enumerate(polys):ck('original_84_residual_'+str(j),f==cv(src['source_polynomials'][j]))
mons=[tuple(map(int,e)) for e in src['monomials']];A=src['coefficient_matrix']
for j,f in enumerate(polys):
    dd={tuple(map(int,e)):c for e,c in f.dict().items()}
    for i,e in enumerate(mons):ck('actual_coefficient_%d_%d'%(i,j),K(A[i,j])==dd.get(e,K.zero()))
ck('target_D_source',D==cv(src['common_amplitude_denominator']))
rhs=vector(P,[P(D.monomial_coefficient(R.monomial(*e))) for e in mons])
print('ORIGINAL 91 NORMAL ->84 SOURCES ->110x84 PASS',len(checks),round(time.monotonic()-st,2),flush=True)
names=['independent_ratio_canonical_second_primitive_certificate.sobj','independent_ratio_canonical_second_exchange_certificate.sobj','independent_ratio_canonical_second_exchange2_certificate.sobj','independent_ratio_canonical_second_exchange_index31_certificate.sobj']
ds=load(str(p/'independent_ratio_canonical_second_three_denominators.sobj'))
four=load(str(p/'independent_ratio_canonical_second_exchange_index31_denominator.sobj'))
bare=[P(f) for f in ds['polynomials']]+[P(four['bare_denominator'])]
denominators=[];units=None;matrix0=None;records=[]
for ci,name in enumerate(names):
    c=load(str(p/name));NN=vector(P,[P(v) for v in c['numerators']]);den=P(c['denominator']);MM=matrix(P,c['normalized_matrix']);uu=[K(v) for v in c['column_units']]
    ck('same_rhs_'+str(ci),vector(P,c['rhs'])==rhs)
    if ci==0:
        units=uu;matrix0=MM
        for j,u in enumerate(uu):
            ck('column_material_num_'+str(j),guarded(u.numerator()));ck('column_material_den_'+str(j),guarded(u.denominator()))
            for i in range(A.nrows()):ck('original_column_%d_%d'%(i,j),K(MM[i,j])*u==A[i,j])
    else:
        ck('same_material_column_units_'+str(ci),uu==units)
        for i in range(MM.nrows()):ck('same_polynomial_matrix_%d_%d'%(ci,i),list(MM.row(i))==list(matrix0.row(i)))
    for i in range(MM.nrows()):
        val=sum((MM[i,j]*NN[j] for j in range(MM.ncols()) if MM[i,j] and NN[j]),P.zero())
        ck('full_original_identity_%d_%d'%(ci,i),val==den*rhs[i])
    ratio=K(den)/K(bare[ci]);ck('bare_den_material_numerator_'+str(ci),guarded(ratio.numerator()));ck('bare_den_material_denominator_'+str(ci),guarded(ratio.denominator()))
    ck('target_only_t2Lambda_'+str(ci),cv(c['target'])==D)
    denominators.append(den);records.append({'source':name,'denominator_degree':int(den.total_degree()),'bare_degree':int(bare[ci].total_degree()),'bare_prefactor':str(ratio.factor())})
    print('FOUR CERTIFICATE',ci,'ALL110 SOURCE ROWS PASS',round(time.monotonic()-st,2),flush=True)
inp=load(str(p/'independent_ratio_canonical_second_resultant_inputs.sobj'));ext=load(str(p/'independent_ratio_canonical_second_extra_resultant_index31_inputs.sobj'))
B=inp['ring'];rr,v=B.gens();conv=list(inp['polynomials'])+[B(ext['polynomials'][1])];recs=list(inp['records'])+[{'r_shift':ext['r_shift'],'rational_scale':ext['rational_scale']}]
ck('extra_shared_first_polynomial',B(ext['polynomials'][0])==B(conv[0]))
for i,(f,g,rec) in enumerate(zip(bare,conv,recs)):
    shift=int(rec['r_shift']);sc=QQ(rec['rational_scale']);rebuild=B.zero()
    for e,c in f.dict().items():rebuild+=sc*c*rr**(int(e[1])-int(e[0])+shift)*v**int(e[0])
    ck('all_parameter_v_equals_ar_'+str(i),rebuild==g)
    ck('converted_integral_'+str(i),all(QQ(c).denominator()==1 for c in g.coefficients()))
    ck('shift_only_r_'+str(i),shift>=0 and sc!=0)
ck('positive_aq1_identity',a*q+1==a*(1+r)+1)
payload={'checks':checks,'converted_ring':B,'converted_polynomials':conv,'material_denominators':denominators,'bare_denominators':bare,'records':records,'target':D,'non_tau_five_minor':minor,'source_indices':indices,'scope':'Four exact original-source identities for D=t2*(t0+t1+t2). Only protected ordered-material factors occur outside the four bare denominators; no all-amplitude exclusion yet until resultant coverage is independently audited.'}
save(payload,str(p/'independent_canonical_main_source_reader.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'original_rows':110,'certificates':4,'sources':84,'records':records,'non_tau_five_minor':str(minor.factor()),'target':str(D),'scope':payload['scope']}
(p/'independent_canonical_main_source_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
