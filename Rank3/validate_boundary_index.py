"""Check valuation, coordinate-support and source references in the boundary index."""
from pathlib import Path
import argparse, hashlib, itertools, json

def validate(package):
    package=Path(package)
    index=json.loads((package/'BOUNDARY_INDEX.json').read_text(encoding='utf-8'))
    proof=json.loads((package/'PROOF_MAP.json').read_text(encoding='utf-8'))
    nodes={n['node']:n for n in proof['nodes']}
    jobs=json.loads((package/'tasks.json').read_text(encoding='utf-8'))['jobs']
    byid={j['id']:j for j in jobs}
    rows=json.loads((package/'MANIFEST.json').read_text(encoding='utf-8'))['files']
    bywork={r['work_path']:r for r in rows}
    assert len(index['valuation_regions'])==5
    assert [r['condition'] for r in index['valuation_regions']]==[
        's_0 = 0','0 < s_0 < n/2','s_0 = n/2','n/2 < s_0 < n','n <= s_0 <= infinity']
    po=['12','13','21','23','31','32']
    expected=set()
    for r in itertools.product([0,1],repeat=3):
        if not any(r):continue
        for p in itertools.product([0,1],repeat=6):
            if any(p) and p[0]*r[0]==p[5]*r[2] and p[1]*r[0]==p[3]*r[1] and p[2]*r[1]==p[4]*r[2]:
                expected.add((r,p))
    actual=set();counts={1:0,2:0,3:0}
    for item in index['intermediate_support_patterns']:
        r=tuple(int(i in item['R_support']) for i in [1,2,3])
        p=tuple(int(i in item['p_support']) for i in po)
        assert (r,p) not in actual
        actual.add((r,p));counts[sum(r)]+=1
        assert item['rank_D']==sum(bool(p[i] or p[j]) for i,j in [(4,5),(2,3),(0,1)])>0
        assert set(item['regions'])=={'supercritical','critical','subcritical'}
        for route in item['regions'].values():
            assert route['nodes'] and route['argument']
            assert all(n in nodes for n in route['nodes'])
    assert actual==expected and counts=={1:45,2:21,3:7}
    catalog={n['node']:n for n in index['proof_nodes']}
    for name,node in catalog.items():
        reference=nodes[name]
        assert node['direct_tasks']==reference['direct_tasks']
        assert node['necessary_dependency_task_ids']==reference['necessary_dependency_task_ids']
        for t in node['direct_tasks']:
            j=byid[t['task_id']]
            assert t['source']==j.get('source_member')
            assert t['input_views']==j.get('input_views',{})
            if j.get('action','run')=='run':
                row=bywork[t['source']]
                assert hashlib.sha256((package/row['archive_path']).read_bytes()).hexdigest()==t['source_sha256']
    for locus in index['exceptional_loci']:
        assert locus['identity'] and locus['result'] and all(n in nodes for n in locus['nodes'])
    visiting=set();done=set()
    def visit(name):
        assert name in byid and name not in visiting
        if name in done:return
        visiting.add(name)
        for dep in byid[name].get('depends_on',[]):visit(dep)
        visiting.remove(name);done.add(name)
    for j in jobs:visit(j['id'])
    return {'status':'INDEX_PASS','valuation_regions':5,'support_patterns':73,
            'intermediate_region_support_pairs':219,'exceptional_loci':len(index['exceptional_loci']),
            'proof_nodes':len(catalog),'registered_tasks':len(jobs),'dependency_graph_acyclic':True}

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    args=a.parse_args();print(json.dumps(validate(args.package),indent=2))
