"""Validate manuscript labels, exact formula excerpts and computation entries."""
from pathlib import Path
import argparse, hashlib, json, re

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def normalize(text):
    text=re.sub(r'(?<!\\)%[^\n]*','',text)
    return re.sub(r'\s+','',text)

def anchor(text,label):
    hits=list(re.finditer(r'\\label\s*\{'+re.escape(label)+r'\}',text))
    if len(hits)!=1: raise ValueError('Expected one manuscript label: '+label)
    pos=hits[0].start()
    mathenv={'equation','equation*','align','align*','gather','gather*',
             'lemma','proposition','theorem','corollary','definition'}
    stack=[]
    for tok in re.finditer(r'\\(begin|end)\{([^}]+)\}',text[:pos]):
        if tok[1]=='begin': stack.append((tok[2],tok.start()))
        else:
            if stack and stack[-1][0]==tok[2]: stack.pop()
    frames=[s for s in stack if s[0] in mathenv]
    if frames:
        env,start=frames[-1]
        close=re.search(r'\\end\{'+re.escape(env)+r'\}',text[pos:])
        if close is None: raise ValueError('Unclosed environment: '+label)
        end=pos+close.end()
    else:
        # A section label binds its whole mathematical subsection, including
        # its unlabelled coefficient displays, rather than a stale line number.
        commands=list(re.finditer(r'\\(?:subsubsection|subsection|section|paragraph)\*?\{',text[:pos]))
        start=commands[-1].start() if commands else pos
        boundary=re.search(r'\\(?:subsubsection|subsection|section)\*?\{',text[pos+hits[0].end()-pos:])
        end=hits[0].end()+boundary.start() if boundary else len(text)
        if end-start>60000: end=hits[0].end()
    value=text[start:end].strip()
    return {'label':label,'line':text.count('\n',0,pos)+1,'text':value,
            'normalized_sha256':hashlib.sha256(normalize(value).encode()).hexdigest()}

def validate(package,paper=None,supplement=None):
    package=Path(package)
    mapping=json.loads((package/'PROOF_MAP.json').read_text())
    tasks=json.loads((package/'tasks.json').read_text())['jobs']
    byid={j['id']:j for j in tasks}
    if len(byid)!=len(tasks): raise ValueError('Duplicate task identifiers')
    manifest=json.loads((package/'MANIFEST.json').read_text())
    bywork={r['work_path']:r for r in manifest['files']}
    if mapping['tasks_sha256']!=digest(package/'tasks.json'): raise ValueError('Task index is outdated')
    roots=set()
    for node in mapping['nodes']:
        for entry in node['direct_tasks']:
            ident=entry['task_id']
            if ident not in byid: raise ValueError('Unregistered indexed task: '+ident)
            j=byid[ident];roots.add(ident)
            for field,actual in [('source',j.get('source_member')),('wrapper',j.get('wrapper_member')),
                                 ('arguments',j.get('arguments',[])),('input_views',j.get('input_views',{})),
                                 ('environment',j.get('environment',{}))]:
                if entry.get(field)!=actual: raise ValueError('Indexed task interface differs: '+ident+' '+field)
            if j.get('action','run')=='run':
                row=bywork[j['source_member']]
                if entry['source_sha256']!=row['sha256'] or digest(package/row['archive_path'])!=row['sha256']:
                    raise ValueError('Indexed source differs: '+ident)
    closure=set(roots);todo=list(roots)
    while todo:
        for dep in byid[todo.pop()].get('depends_on',[]):
            if dep not in byid: raise ValueError('Unregistered dependency: '+dep)
            if dep not in closure: closure.add(dep);todo.append(dep)
    if closure!=set(byid): raise ValueError('Tasks without a proof use: '+repr(sorted(set(byid)-closure)))
    for item in mapping['formula_inputs']:
        path=package/item['path']
        if digest(path)!=item['sha256']: raise ValueError('Formula input differs: '+item['path'])
        value=path.read_text(encoding='utf-8-sig')
        for block in item.get('exact_blocks',[]):
            if normalize(block) not in normalize(value): raise ValueError('Formula block is absent: '+item['path'])
    if paper:
        text=Path(paper).read_text(encoding='utf-8-sig')
        for node in mapping['nodes']:
            for item in node['manuscript_anchors']:
                if anchor(text,item['label'])['normalized_sha256']!=item['normalized_sha256']:
                    raise ValueError('Manuscript mathematics differs: '+item['label'])
        for item in mapping['formula_inputs']:
            if item.get('manuscript')=='main':
                for block in item.get('exact_blocks',[]):
                    if normalize(block) not in normalize(text): raise ValueError('Printed formula differs from main manuscript')
    if supplement:
        tex=Path(supplement).read_text(encoding='utf-8-sig')
        for item in mapping['formula_inputs']:
            if item.get('manuscript')=='supplement':
                for block in item.get('exact_blocks',[]):
                    if normalize(block) not in normalize(tex): raise ValueError('Printed matrix supplement differs')
    return {'status':'INDEX_PASS','proof_nodes':len(mapping['nodes']),'indexed_roots':len(roots),
            'registered_tasks':len(tasks),'dependency_closure':len(closure),
            'paper_sha256':mapping['paper']['sha256'],'formula_inputs':len(mapping['formula_inputs'])}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--package',type=Path,default=Path(__file__).resolve().parent)
    p.add_argument('--paper',type=Path)
    p.add_argument('--supplement',type=Path)
    p.add_argument('--output',type=Path)
    a=p.parse_args();result=validate(a.package,a.paper,a.supplement)
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
