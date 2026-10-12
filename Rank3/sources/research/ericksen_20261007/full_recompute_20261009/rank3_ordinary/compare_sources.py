"""Bind the complete incoming sources to the reconstructed state equations."""
from sage.all import *
from pathlib import Path
import hashlib,json,time
from job_cache import valid_report
OUT=Path(__file__).absolute().parent
NEW=OUT/'rebuild/round10'; RESEARCH=OUT.parent.parent
checks=[]; st=time.monotonic()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def ck(name,value):
    assert value,name
    checks.append(name)
def bind_full_source_inputs():
    fresh=OUT.parent/'rank3'
    state_paths=[fresh/'fresh_global_polynomials.sobj',fresh/'fresh_joint_fields.sobj',fresh/'fresh_global_source_DAG.sobj']
    fg,fj,nh=[load(str(p)) for p in state_paths]
    R=fg['ring']; K=R.base_ring(); x=R.gens()
    filenames={
        'gauge':'independent_ratio_p_boundary_second_gauge.sobj',
        'H':'independent_ratio_global_R_derivative_degreefall.sobj',
        'E':'independent_ratio_even_derivative_cyclic_direct.sobj',
        'Q':'rank3_joint_fiber_leading.sobj'}
    paths={role:{} for role in filenames}; receipts=[]
    Q0=NEW/'independent_ratio_raw_Q0_normal_two.sobj'
    for stem in ('independent_ratio_incoming12_generic_normal',
                 'independent_ratio_incoming12_generic_normal_two',
                 'independent_ratio_incoming12_all_third_range'):
        path=OUT/'job_reports'/(stem+'.json'); ck(stem+' current receipt',valid_report(path))
        rec=json.loads(path.read_text()); source=RESEARCH/'round10'/(stem+'.py')
        ck(stem+' complete source',rec['executed_source_sha256']==sha(source))
        for role,name in filenames.items():
            matches=[r for r in rec['imports'] if Path(r['path']).name==name]
            ck(stem+' actual input '+role,bool(matches))
            for item in matches:
                pp=Path(item['path']).resolve()
                ck(stem+' input bytes '+role,pp.is_relative_to(RESEARCH.resolve()) and pp.is_file() and sha(pp)==item['sha256'])
                paths[role][str(pp)]=item['sha256']
        if stem.endswith('all_third_range'):
            matches=[r for r in rec['imports'] if Path(r['path']).resolve()==Q0.resolve()]
            ck('complete fresh Q0 source',bool(matches) and all(r.get('regenerated') is True and r['sha256']==sha(Q0) for r in matches))
        receipts.append({'source':str(source),'source_sha256':sha(source),'receipt':str(path),'receipt_sha256':sha(path)})
    def coefficient(c):
        if c in QQ:return K(c)
        C=c.parent(); mp=C.ring().hom(list(K.gens()),K)
        return mp(c.numerator())/mp(c.denominator())
    def polynomial(value):
        dd=value if isinstance(value,dict) else value.dict()
        return R({tuple(map(int,e)):coefficient(c) for e,c in dd.items()})
    def layer(f,n):return R({tuple(e):c for e,c in f.dict().items() if sum(e)==n})
    supports=(((0,3,1),(1,4,1)),((2,2,1),(3,5,-1)),((4,0,1),(5,1,-1)))
    source_records=None
    for role,objects in paths.items():
        for path in objects:
            d=load(path)
            if role=='gauge':
                ck('all six R',len(d['six_full_sources'])==6)
                for i in range(6):ck('actual R'+str(i),polynomial(d['six_full_sources'][i])==fg['sources'][i])
                for i in range(3):
                    for j in range(9):ck('actual gauged field'+str((i,j)),polynomial(d['first_gauged_fields'][i,j])==fg['fields'][i,j])
                    for j in range(2):ck('actual base field'+str((i,j)),polynomial(d['base_fields'][i,j])==fg['base_fields'][i,j])
                    ck('actual gauge unit'+str(i),coefficient(d['source_unit_coefficients'][i])==nh['units'][i])
            elif role=='Q':
                ck('all21 Q',len(d['quartic_source_rows'])==len(fj['quartic_source_rows'])==21)
                for i in range(21):ck('actual Q'+str(i),polynomial(d['quartic_source_rows'][i])==fj['quartic_source_rows'][i])
            elif role=='H':
                ck('all18 H',len(d['records'])==len(nh['records'])==18)
                opool=[polynomial(dd) for dd in d['polynomial_pool']]; npool=[R(dd) for dd in nh['polynomial_pool']]
                for i,(o,n) in enumerate(zip(d['records'],nh['records'])):
                    od=o['full_signed_product_DAG']; nd=n['full_signed_product_DAG']
                    ck('full H term count'+str(i),len(od)==len(nd))
                    ck('actual H index'+str(i),(o['direction'],o['source_position'])==(n['direction'],n['source_position']))
                    for j,((os,ol,orr),(ns,nl,nr)) in enumerate(zip(od,nd)):
                        ck('actual H summand'+str((i,j)),os==ns and opool[ol]==npool[nl] and opool[orr]==npool[nr])
                source_records=[{'kind':'Q','index':i} for i in (0,1,5,6)]+[{'kind':'R','index':i} for i in range(6)]
                source_records += [{'kind':'H','direction':d['records'][i]['direction'],'index':d['records'][i]['source_position']} for i in range(18)]
            elif role=='E':
                ck('all9 E',d['supports']==supports and len(d['items'])==9)
                for row,it in enumerate(d['items']):
                    i,qi=int(it['direction']),int(it['Q_source_index'])
                    derivative=layer(fj['quartic_source_rows'][qi],4).derivative(x[6+i]); co=[{},{}]
                    for ee,c in derivative.dict().items():
                        e=list(map(int,ee)); k=next(k for k,(pos,_,_) in enumerate(supports[i]) if e[pos]);e[supports[i][k][0]]-=1;co[k][tuple(e)]=c
                    cofactors=[R(v) for v in co]
                    ck('complete E derivative'+str(row),derivative==sum(cofactors[k]*x[pos] for k,(pos,_,_) in enumerate(supports[i])))
                    for k in range(2):ck('actual E cofactor'+str((row,k)),polynomial(it['global_support_coefficient_dicts'][k])==cofactors[k])
                for i,qi in enumerate((1,5,6)):ck('actual E full Q'+str(qi),polynomial(d['original_full_Q_coefficients'][i])==fj['quartic_source_rows'][qi])
                e_records=[{'kind':'E','direction':int(it['direction']),'index':int(it['Q_source_index'])} for it in d['items']]
    ck('source ordering',source_records is not None)
    return {'source_inputs':[{'role':role,'path':p,'sha256':h} for role,items in paths.items() for p,h in items.items()],
            'state_sources':[{'path':str(p),'sha256':sha(p)} for p in state_paths],
            'producer_receipts':receipts,'complete_Q0_sha256':sha(Q0)},source_records+e_records
binding,records=bind_full_source_inputs()
first=load(str(NEW/'independent_ratio_incoming12_generic_normal.sobj'))
second=load(str(NEW/'independent_ratio_incoming12_generic_normal_two.sobj'))
null=load(str(NEW/'independent_ratio_incoming12_generic_nullwall_all.sobj'))
third=load(str(NEW/'independent_ratio_incoming12_all_third_raw.sobj'))
ck('all37 first and second sources',first['matrix'].nrows()==second['coefficient_matrix'].nrows()==37)
ck('complete source order',first['source_records']==second['source_records']==records and len(records)==37)
ck('all37 nullwall sources',null['complete'] and len(null['items'])==37 and {int(it['source_row']) for it in null['items']}==set(range(37)))
ck('all fifteen coefficients',all(len(it['rhs'])==15 for it in null['items']))
ck('all37 third sources',len(third['items'])==37 and {int(it['source_row']) for it in third['items']}==set(range(37)))
ck('full raw cubic amplitude algebra',third['basis']==['1','l','y','l*y','y^2','l*y^2','y^3','l*y^3'] and all(len(it['constant'])==8 and all(len(v)==8 for v in it['unknown_coefficients']+it['third_range_coefficients']) for it in third['items']))
ck('all third range and kernel coefficients',all(len(it['unknown_coefficients'])==5 and len(it['third_range_coefficients'])==6 for it in third['items']))
direct_path=OUT.parent/'rank3_fibres/ordinary_H_direct_coefficients.sobj'
direct=load(str(direct_path)); current={int(it['source_row']):it for it in null['items']}
required={(11,1),(11,4),(13,1),(13,4),(23,1),(25,1)}
entries={(int(it['row']),int(it['column'])):it['source_value'] for it in direct['entries']}
ck('eight full H coefficients',set(entries)=={(r,c) for r in (11,13,23,25) for c in (1,4)} and required<=set(entries))
ck('same amplitude algebra',direct['relation']==null['relation'])
for (row,col),value in entries.items():ck('direct full H coefficient'+str((row,col)),value==current[row]['rhs'][col])
ck('direct H source input',any(r['role']=='H' and r['path']==direct['H_input']['path'] and r['sha256']==direct['H_input']['sha256'] for r in binding['source_inputs']))
ck('direct current nullwall input',direct['nullwall_sha256']==sha(NEW/'independent_ratio_incoming12_generic_nullwall_all.sobj'))
report={'status':'PASS','count':len(checks),'checks':checks,'source_binding':binding,
        'complete_source_rows':37,'direct_H_coefficients':8,'covered_H_positions':sorted(required),
        'seconds':time.monotonic()-st,'coefficient_object':str(direct_path),'coefficient_object_sha256':sha(direct_path)}
(OUT/'source_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',len(checks),'complete incoming source identities',flush=True)
