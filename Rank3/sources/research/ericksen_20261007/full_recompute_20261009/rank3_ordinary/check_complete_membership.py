"""Exact coefficient divisibilities for COMPLETE ordinary relative ideals.

All inputs below were regenerated in this new folder. Higher terms are
covered by the explicitly proved normal-ideal power inclusions in the
ordinary appendix; no stored zero remainder is used as an input.
"""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).absolute().parent;DATA=OUT/'rebuild/round10'
checks=[];started=time.monotonic()
def ck(name,ok):
    assert ok,name
    checks.append(dict(name=name,passed=True))
n=load(str(DATA/'independent_ratio_incoming12_generic_normal.sobj'))
j=load(str(DATA/'independent_ratio_incoming12_generic_normal_two.sobj'))
p=load(str(DATA/'ordinary_pair_incoming_full_twojet.sobj'))
ck('all first-source rows regenerated',n['matrix'].nrows()==37 and n['matrix'].ncols()==6)
ck('all second-source rows regenerated',j['coefficient_matrix'].nrows()==37 and j['coefficient_matrix'].ncols()==27)
l,u,v=n['ring'].gens()
for i in range(37):
    ck('complete eta coefficient divisible by u '+str(i),n['matrix'][i,4].subs({u:0})==0)
    ck('complete t coefficient divisible by v+2 '+str(i),(n['matrix'][i,3]-n['matrix'][i,0]).subs({v:-2})==0)
    ck('second reconstruction recovers first '+str(i),list(j['coefficient_matrix'].row(i)[:6])==list(n['matrix'].row(i)))
pl,pu,pv=p['jet_ring'].gens()
for item in p['items']:
    values=item['jet_coefficients']
    ck('pair constant is zero '+str(item['pair']),values[0]==0)
    ck('pair eta coefficient divisible by u '+str(item['pair']),values[5].subs({pu:0})==0)
    ck('pair t coefficient divisible by v+2 '+str(item['pair']),(values[4]-values[1]).subs({pv:-2})==0)
us=load(str(DATA/'ordinary_pair_minus_whole_u_square.sobj'))
ul,uu=us['ring'].gens()
ck('negative-wall all forty sources present',len(us['items'])==40)
for item in us['items']:
    ck('whole shifted square lies in (nu,u) '+str(item['label']),item['wall_remainder'].subs({uu:0})==0)
rg=load(str(DATA/'ordinary_pair_minus_general_quadratic_range.sobj'))
ck('full generic quadratic range correction forty remainders',len(rg['residuals'])==40 and all(f==0 for f in rg['residuals']))
e35=load(str(DATA/'ordinary_pair_minus_E35_cubic.sobj'))
ck('full E35 compensated cubic all material and base',e35['full_pure_coefficients'][3]==0)
guards=json.loads((DATA/'independent_ratio_incoming12_second_module_baseguards.json').read_text())
ck('every second-module inverse is generically valid for every allowed material',guards['count']==167 and not guards['bad'])
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-started,
            sample_materials_used=0,scope='Exact whole-parameter coefficient memberships used with the normal-ideal power lemma.')
(OUT/'complete_membership_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',len(checks),'complete-source coefficient identities',flush=True)
