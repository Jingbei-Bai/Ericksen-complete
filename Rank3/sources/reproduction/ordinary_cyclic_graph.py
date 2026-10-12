"""Recover the cyclic rational graph from the complete ordinary R equations."""
from pathlib import Path
import argparse, json, time
from sage.all import QQ, PolynomialRing, load, save, matrix, vector, prod

parser = argparse.ArgumentParser()
parser.add_argument('--work-root', type=Path, required=True)
root = parser.parse_args().work_root.resolve()
started = time.monotonic()
rank = root / 'research/ericksen_20261007/full_recompute_20261009/rank3'
fresh = load(str(rank / 'fresh_global_polynomials.sobj'))
joint = load(str(rank / 'fresh_joint_fields.sobj'))
R = fresh['ring']
oldK = R.base_ring()
K = PolynomialRing(QQ, names=('cyc_a', 'cyc_r')).fraction_field()
a, r = K.gens()
q = 1 + r
arguments = (a, a*q/(1+a*r), 1/r, -q/r)
phi = oldK.ring().hom(arguments, K)
T = PolynomialRing(K, names=R.variable_names()[:6])
F = T.fraction_field()
checks = []
cache = {}

def ck(name, value):
    assert value, name
    checks.append(name)

def coefficient(c):
    c = oldK(c)
    if c not in cache:
        cache[c] = phi(c.numerator()) / phi(c.denominator())
    return cache[c]

def pure(f):
    ck('coefficient polynomial depends only on p',
       all(not any(e[6:]) for e in f.dict()))
    return T({tuple(map(int, e[:6])): coefficient(c)
              for e, c in f.dict().items()})

C = matrix(T, 3, 3, [pure(v) for v in fresh['C'].list()])
L = matrix(T, 3, 3, [pure(v) for v in fresh['L'].list()])
offset = vector(T, [pure(v) for v in fresh['offset']])
ck('both graph chart determinants nonzero', C.det() != 0 and L.det() != 0)
rr = -(C*L).change_ring(F).inverse() * offset.change_ring(F)
ck('complete graph system', (C*L).change_ring(F)*rr + offset == 0)
den = rr[0].denominator()
for value in rr[1:]:
    den = den.lcm(value.denominator())
ck('common graph denominator nonzero', den != 0)
numerators = []
for j, value in enumerate(rr):
    product = value * F(den)
    ck('exact common denominator '+str(j), product.denominator().is_constant())
    numerators.append(T(product))

def pull(f, degree, rr_degree):
    result = T.zero()
    for e, c in f.dict().items():
        if sum(e) != degree:
            continue
        re = sum(e[6:])
        ck('graph substitution RR degree bound', re <= rr_degree)
        result += T({tuple(map(int, e[:6])): coefficient(c)}) * den**(rr_degree-re) * prod(
            numerators[j]**int(e[6+j]) for j in range(3))
    return result

for j, source in enumerate(fresh['sources']):
    ck('all six complete ordinary R tops '+str(j), pull(source, 4, 1) == 0)
ck('joint and ordinary coefficient fields agree', joint['ring'].base_ring() == oldK)
items = []
for index in (0, 1, 5, 6):
    source = joint['quartic_source_rows'][index]
    degree = max(map(sum, source.dict()))
    rr_degree = max(sum(e[6:]) for e in source.dict() if sum(e) == degree)
    numerator = pull(source, degree, rr_degree)
    items.append(dict(source=index, numerator=numerator, denominator_power=rr_degree))
target = root / 'research/ericksen_20261007/round10/independent_ratio_p_graph_original_sources_cyclic.sobj'
target.parent.mkdir(parents=True, exist_ok=True)
payload = dict(coefficient_field=K, p_ring=T, graph_denominator=den,
               graph_numerators=numerators, items=items,
               original_chart_determinants=(L.det(), C.det()), cyclic_base_only=True)
save(payload, str(target), compress=False)
restored = load(str(target))
ck('exact graph and quartic readback', restored == payload)
report = dict(status='PASS', count=len(checks), checks=checks,
              seconds=time.monotonic()-started)
(target.parent/'fresh_ordinary_cyclic_graph.json').write_text(json.dumps(report, indent=2)+'\n')
print('PASS cyclic graph', len(checks), round(report['seconds'], 2), flush=True)
