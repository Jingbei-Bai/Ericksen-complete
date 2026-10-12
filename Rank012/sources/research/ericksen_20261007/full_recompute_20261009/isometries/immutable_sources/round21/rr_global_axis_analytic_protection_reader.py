"""Independent pointwise-parameter protection of the normalized axis system.

Generic e-valuation alone is insufficient.  This checks all coefficient poles
and the leading denominator in e against the stated whole-material guard.
"""
from sage.all import *
from pathlib import Path
import json, signal, time

signal.alarm(120)
start = time.time()
O = Path(__file__).resolve().parent
cp = load(str(O / 'rr_null_physical_axis_reduction_checkpoint_plus.sobj'))
mat = cp['reduced_matrix']
FE = mat.base_ring()
E = FE.ring()
e = E.gen()
F = E.base_ring()
B = F.ring()
s, t = B.gens()
guard = (s*t*(s*s-1)*(t*t-1)*(s*s-t*t)*(s*s+t*t)
         *(s*s+t*t-2)*((s*s-1)**2+(t*t-1)**2))
normal = diagonal_matrix(FE, [e**4, e**7, e**(-2)]) * mat
cache = {}

def supported(poly):
    poly = B(poly)
    if not poly:
        return False
    if poly in cache:
        return cache[poly]
    z = poly
    while z.total_degree() > 0:
        g = z.gcd(guard)
        if g.total_degree() == 0:
            cache[poly] = False
            return False
        z = z.quo_rem(g)[0]
    cache[poly] = True
    return True

checks = {}
for i in range(3):
    for j in range(6):
        f = normal[i,j]
        if not f:
            checks['entry_%s_%s_zero' % (i,j)] = True
            continue
        n, d = f.numerator(), f.denominator()
        vd = d.valuation()
        # F[e]'s monic rational representation can create fake parameter
        # poles.  Clear both numerator and denominator together, then remove
        # their common B-content before testing the actual e-leading term.
        coeffs = list(n) + list(d)
        common_den = B.one()
        for v in coeffs:
            common_den = common_den.lcm(F(v).denominator())
        cleared = [B(F(v)*common_den) for v in coeffs]
        content = gcd(cleared)
        lead = B(F(d[vd])*common_den).quo_rem(content)[0]
        ok = supported(lead)
        ok = ok and n.valuation() >= vd
        checks['entry_%s_%s_uniform_analytic' % (i,j)] = bool(ok)
assert all(checks.values()), (checks, [str(p.factor()) for p,b in cache.items() if not b])
report = {
    'status': 'PASS_WHOLE_MATERIAL_AXIS_ANALYTIC_PROTECTION',
    'checks': checks, 'count': len(checks),
    'guard': str(guard.factor()),
    'guard_domain': '|t| > s > 1; all listed real factors are nonzero',
    'scope': 'Each normalized full rational coefficient is analytic near e=0 at every allowed real material point, not merely regular over the generic parameter field.',
    'distinct_polynomial_support_checks': len(cache),
    'seconds': time.time()-start,
}
(O / 'rr_global_axis_analytic_protection_reader.json').write_text(
    json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps(report, indent=2), flush=True)
