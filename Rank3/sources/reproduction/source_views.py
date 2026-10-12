"""Build exact legacy source schemas solely from complete fresh global sources.

This is a new mathematical schema adapter, not a byte-exact execution claim
for any historical generator. All polynomial coefficients and supports are
retained. The separate homogeneous layers are the same auxiliary layers used
by the original schemas; the full H identity is never replaced by a layer.

Producer: rank3_sources:fresh_global_sources.
Inputs: fresh_global_source_DAG.sobj and fresh_global_polynomials.sobj.
Outputs: three complete, producer-bound original-name SOBJ/JSON pairs.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

from sage.all import PolynomialRing, identity_matrix, load, matrix, save, vector

PREFIX = Path('research/ericksen_20261007')
FRESH = PREFIX / 'full_recompute_20261009/rank3'
PRODUCER_ID = 'rank3_sources:fresh_global_sources'
NAMES = {
    'global': 'independent_ratio_global_R_derivative_degreefall',
    'second': 'independent_ratio_p_boundary_second_gauge',
    'support': 'independent_ratio_p_boundary_support_gauges',
}
REQUIRED_KEYS = {
    'global': ('coefficient_field', 'polynomial_pool', 'records', 'C_refs',
               'D_refs', 'L_refs', 'Z_source_refs', 'full_base_field_refs',
               'full_first_gauged_field_refs', 'full_source_refs', 'sigma',
               'source_order', 'source_top4_refs'),
    'second': ('coefficient_field', 'first_gauged_fields', 'base_fields',
               'six_full_sources', 'source_Q0', 'source_permutation',
               'source_unit_coefficients'),
    'support': ('ring', 'coefficient_field', 'D', 'items', 'Zvec',
                'Z_state_jacobian', 'Z_base_jacobian'),
}
CONSTRUCTION_REFERENCES = {
    'global': 'round10/independent_ratio_global_R_derivative_degreefall.py:15-115',
    'second': 'round10/independent_ratio_p_boundary_second_gauge.py:22-130; '
              'round10/independent_ratio_p_boundary_graph_cover.py:12-17',
    'support': 'round10/independent_ratio_p_boundary_support_gauges.py:16-76; '
               'round10/independent_ratio_p_boundary_second_principal_tangency.py:22-26',
}
CONSTRUCTION_SOURCE_SHA256 = {
    'round10/independent_ratio_global_R_derivative_degreefall.py':
        '29b08324275aac8e7ee2c971b884cba6a260fadeec68f7a49ecf7130c41f7534',
    'round10/independent_ratio_p_boundary_second_gauge.py':
        '73e6f2c2356f0af2a42774f4ccba574d601882e4f6d99a6a15c5c21659444ded',
    'round10/independent_ratio_p_boundary_graph_cover.py':
        '630d9e2e612fae1c52c12013e9cec819af27e6075686dda5f0685686c32d5288',
    'round10/independent_ratio_p_boundary_support_gauges.py':
        'f885ecdd5379b5f4f582dd00d83d8ae28cab452d3871a7e70732a2a3ab067c75',
    'round10/independent_ratio_p_boundary_second_principal_tangency.py':
        'a6476e2806c8bb524593c3f7887e9123b1fdebdddb11859cc7e2f1b25f2cab1b',
    'full_recompute_20261009/rank3/fresh_global_sources.py':
        '81d0e80003401db934e3076f3458f6d95b9a6054d8a64b27d1cb66766539d362',
}


def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


class FreshViews:
    def __init__(self, root):
        self.root = root
        self.start = time.monotonic()
        self.paths = [root / FRESH / 'fresh_global_source_DAG.sobj',
                      root / FRESH / 'fresh_global_polynomials.sobj']
        self.inputs = [{'path': str(p), 'sha256': file_sha(p)} for p in self.paths]
        self.dag = load(str(self.paths[0]))
        self.fresh = load(str(self.paths[1]))
        self.native = self.fresh['ring']
        self.K = self.native.base_ring()
        self.a, self.b, self.h, self.z = self.K.gens()
        self.r = (self.a - self.b) / (self.a * (self.b - 1))
        self.q = 1 + self.r
        self.H = self.r * self.h
        self.Z = self.r * self.z / self.q
        self.alpha = (self.Z - self.H) / (self.Z * (self.H + 1))
        self.beta = self.H * (self.Z - 1) / (self.Z * (self.H + 1))
        self.units = tuple(self.K(c) for c in self.dag['units'])
        self.sigma = tuple(int(c) for c in self.dag['sigma'])
        self.order = tuple(int(c) for c in self.dag['source_order'])
        self.checks = []
        self.ck('fresh coefficient field agrees in both inputs',
                self.dag['coefficient_field'] == self.K)
        self.ck('nine ordered state variables', self.native.ngens() == 9)
        self.ck('six original source order', self.order == (10, 11, 13, 15, 18, 19))
        self.ck('three canonical source permutations', self.sigma == (2, 1, 0))
        self.ck('six complete R sources', len(self.fresh['sources']) == 6)
        self.ck('all eighteen complete derivative records', len(self.dag['records']) == 18)
        self.ck('fresh units match the complete original formulas', self.units == (
            2 * (self.H - self.Z) / (self.H * self.Z),
            -2 * (self.Z - 1) / self.Z,
            self.r * (self.H + 1) / (self.H * self.q)))

    def ck(self, label, value):
        if not bool(value):
            raise AssertionError(label)
        self.checks.append(label)

    def convert(self, ring, poly):
        # Complete coefficient/exponent transport. There is no support filter.
        return ring({tuple(int(e) for e in ex): self.K(c)
                     for ex, c in poly.dict().items() if c})

    def data(self, ring):
        cp = lambda f: self.convert(ring, f)
        sources = [cp(f) for f in self.fresh['sources']]
        hat = matrix(ring, 3, 9, [cp(f) for f in self.fresh['fields'].list()])
        base = matrix(ring, 3, 2, [cp(f) for f in self.fresh['base_fields'].list()])
        D = matrix(ring, 6, 3, [cp(f) for f in self.fresh['D'].list()])
        C = matrix(ring, 3, 3, [cp(f) for f in self.fresh['C'].list()])
        L = matrix(ring, 3, 3, [cp(f) for f in self.fresh['L'].list()])
        offset = vector(ring, [cp(f) for f in self.fresh['offset']])
        Zvec = C * L * vector(ring, ring.gens()[6:]) + offset
        self.ck('complete factored source rows in ' + str(ring.variable_names()),
                vector(ring, [self.layer(ring, f, 4) for f in sources]) == D * Zvec)
        self.ck('all six full source degree bounds', all(self.degree(f) <= 4 for f in sources))
        for i in range(3):
            for j in range(9):
                self.ck('full first gauge degree %d %d' % (i, j),
                        self.degree(hat[i, j]) <= (3 if j == 6 + i else 2))
            for j in range(2):
                self.ck('full base degree %d %d' % (i, j), self.degree(base[i, j]) <= 1)
            self.ck('complete cubic diagonal %d' % i,
                    self.layer(ring, hat[i, 6 + i], 3) == self.units[i] * Zvec[self.sigma[i]])
        return dict(ring=ring, sources=sources, hat=hat, base=base,
                    D=D, C=C, L=L, offset=offset, Zvec=Zvec,
                    Q0=cp(self.fresh['Q0']))

    @staticmethod
    def layer(ring, f, degree):
        return ring({tuple(int(v) for v in e): c for e, c in f.dict().items()
                     if sum(e) == degree})

    @staticmethod
    def degree(f):
        return max((sum(e) for e in f.dict()), default=-1)

    @staticmethod
    def hist(f):
        result = {}
        for e in f.dict():
            degree = sum(int(v) for v in e)
            result[degree] = result.get(degree, 0) + 1
        return result

    @staticmethod
    def bd(ring, f, parameter):
        answer = {}
        for e, c in f.dict().items():
            derivative = c.derivative(parameter)
            if derivative:
                answer[tuple(int(v) for v in e)] = derivative
        return ring(answer)

    def global_view(self):
        R = PolynomialRing(self.K, names=tuple('glob_p%d' % i for i in range(6)) +
                           tuple('glob_rr%d' % i for i in range(3)))
        d = self.data(R)
        xx = R.gens()
        src, hat, base, D, C, L, Zvec = (d[k] for k in
            ('sources', 'hat', 'base', 'D', 'C', 'L', 'Zvec'))
        tops = [self.layer(R, f, 4) for f in src]
        lower3 = [self.layer(R, f, 3) for f in src]
        partials = [[f.derivative(v) for v in xx] for f in src]
        basepartials = [[self.bd(R, f, p) for p in (self.h, self.z)] for f in src]
        hat2 = [[self.layer(R, hat[i, j], 2) for j in range(9)] for i in range(3)]
        base1 = [[self.layer(R, base[i, j], 1) for j in range(2)] for i in range(3)]
        top_partials = [[f.derivative(v) for v in xx] for f in tops]
        top_basepartials = [[self.bd(R, f, p) for p in (self.h, self.z)] for f in tops]
        ecoef = []
        for i in range(3):
            rows = []
            for k in range(3):
                v = [R.zero() for _ in range(6)]
                if i == 0:
                    v[3], v[4] = C[k, 0], C[k, 1]
                elif i == 1:
                    v[5], v[2] = self.alpha / self.beta * C[k, 1], C[k, 2]
                else:
                    v[1], v[0] = self.alpha * C[k, 0], -self.beta * C[k, 2]
                self.ck('exact full E leading module %d %d' % (i, k),
                        sum((v[j] * D[j, self.sigma[i]] for j in range(6)), R.zero()) ==
                        (C * L.column(i))[k])
                for j in range(6):
                    if v[j]:
                        self.ck('E selects only its source column %d %d %d' % (i, k, j),
                                all(D[j, ell] == 0 for ell in range(3) if ell != self.sigma[i]))
                rows.append(v)
            ecoef.append(rows)
        pool, refs = [], {}

        def ref(f):
            f = R(f)
            if f not in refs:
                refs[f] = len(pool)
                pool.append({tuple(int(v) for v in e): c for e, c in f.dict().items()})
            return refs[f]

        def rr(v):
            return [ref(f) for f in v]

        oldpool = [R({tuple(int(v) for v in e): self.K(c) for e, c in dd.items()})
                   for dd in self.dag['polynomial_pool']]
        records = []
        for i in range(3):
            for j in range(6):
                corr = [sum((D[j, k] * ecoef[i][k][ell] for k in range(3)), R.zero())
                        for ell in range(6)]
                full = [(1, hat[i, v], partials[j][v]) for v in range(9) if partials[j][v]]
                full += [(1, base[i, v], basepartials[j][v]) for v in range(2) if basepartials[j][v]]
                full += [(-1, self.units[i] * corr[ell], src[ell]) for ell in range(6) if corr[ell]]
                saved = self.dag['records'][6 * i + j]
                self.ck('fresh source record coordinates %d %d' % (i, j),
                        (int(saved['direction']), int(saved['source_position'])) == (i, j))
                expected = [(int(sign), oldpool[int(left)], oldpool[int(right)])
                            for sign, left, right in saved['full_signed_product_DAG']]
                # Factor-level equality verifies the complete identity without
                # expanding and recancelling any degree-six polynomial product.
                self.ck('complete signed factor DAG matches fresh producer %d %d' % (i, j),
                        full == expected)
                top = [(1, hat2[i][v], top_partials[j][v]) for v in range(9)
                       if hat2[i][v] and top_partials[j][v]]
                top += [(1, base1[i][v], top_basepartials[j][v]) for v in range(2)
                        if base1[i][v] and top_basepartials[j][v]]
                top += [(1, self.units[i] * Zvec[self.sigma[i]], lower3[j].derivative(xx[6 + i]))]
                top += [(-1, self.units[i] * corr[ell], lower3[ell]) for ell in range(6) if corr[ell]]
                for _, left, right in top:
                    self.ck('complete auxiliary degree-five term %d %d' % (i, j),
                            not left or not right or self.degree(left) + self.degree(right) == 5)
                self.ck('full degree-six cancellation %d %d' % (i, j),
                        tops[j].derivative(xx[6 + i]) == sum(
                            (D[j, k] * (C * L.column(i))[k] for k in range(3)), R.zero()))
                records.append(dict(direction=i, source_position=j,
                    original_source_index=self.order[j], degree_bound=5,
                    full_signed_product_DAG=[(sign, ref(a), ref(b)) for sign, a, b in full],
                    ordinary_degree5_signed_product_DAG=[(sign, ref(a), ref(b)) for sign, a, b in top],
                    correction_coefficients=rr(corr)))
        payload = dict(format='portable polynomial pool; each ref is an integer index; '
                              'DAG triples=(sign,left_ref,right_ref)',
            coefficient_field=self.K, original_parameters=(self.a, self.b, self.h, self.z),
            parameter_images=(self.a, self.r, self.H, self.Z), variable_names=R.variable_names(),
            D_refs=[rr(row) for row in D.rows()], C_refs=[rr(row) for row in C.rows()],
            L_refs=[rr(row) for row in L.rows()], Z_source_refs=rr(Zvec),
            units=self.units, sigma=self.sigma, source_order=self.order,
            full_source_refs=rr(src), source_top4_refs=rr(tops), source_layer3_refs=rr(lower3),
            full_first_gauged_field_refs=[rr(row) for row in hat.rows()],
            full_base_field_refs=[rr(row) for row in base.rows()],
            source_state_partial_refs=[rr(row) for row in partials],
            source_base_partial_refs=[rr(row) for row in basepartials],
            E_source_coefficient_refs=[[rr(row) for row in rows] for rows in ecoef],
            records=records, polynomial_pool=pool,
            identity='H_ij=hatFi(Rj)-u_i sum_k D_jk E_ik(R), with full h,z derivatives.',
            initial_identity='Full ordinary degree-six cancellation; degree-five DAG includes every term.',
            scope='Complete original 6R/18H identity schema from fresh necessary sources; no state-pivot division.')
        return payload, {'sources': 6, 'records': 18, 'state_partials_shape': [6, 9],
                         'base_partials_shape': [6, 2], 'E_coefficients_shape': [3, 3, 6],
                         'D_shape': [6, 3], 'C_shape': [3, 3], 'L_shape': [3, 3],
                         'pool_entries': len(pool), 'variable_names': list(R.variable_names())}

    def second_view(self):
        R = PolynomialRing(self.K, names=('sg_p12', 'sg_p13', 'sg_p21', 'sg_p23',
                                         'sg_p31', 'sg_p32', 'sg_rr1', 'sg_rr2', 'sg_rr3'))
        d = self.data(R)
        xx = R.gens()
        W = matrix(R, 3, 6)
        # Literal original graph_cover numerator, in ordered p12,p13,p21,p23,p31,p32.
        W[0, 0], W[0, 1] = self.alpha * xx[0] * xx[3], self.alpha * xx[1] * xx[2]
        W[1, 2], W[1, 5] = -self.alpha * xx[1] * xx[5], -self.alpha * xx[0] * xx[4]
        W[2, 3], W[2, 4] = self.alpha * xx[3] * xx[4], -self.alpha * xx[2] * xx[5]
        Delta = d['L'].det()
        self.ck('complete Delta has ordinary degree three', self.degree(Delta) == 3)
        self.ck('complete Delta D inverse without individual coordinate division',
                W * d['D'] == Delta * identity_matrix(R, 3))
        src = vector(R, d['sources'])
        WR = W * src
        hat, base = d['hat'], d['base']
        numerator = Delta * hat
        for i in range(3):
            numerator[i, 6 + i] -= self.units[i] * WR[self.sigma[i]]
        base_numerator = Delta * base
        records = []
        for i in range(3):
            self.ck('full source gauge leading WR %d' % i,
                    self.layer(R, WR[self.sigma[i]], 6) == Delta * d['Zvec'][self.sigma[i]])
            for j in range(9):
                self.ck('complete gauge degree bound %d %d' % (i, j),
                        self.degree(numerator[i, j]) <= 5)
                self.ck('whole gauge source identity %d %d' % (i, j),
                        numerator[i, j] == Delta * hat[i, j] - (
                            self.units[i] * WR[self.sigma[i]] if j == 6 + i else 0))
            for j in range(2):
                self.ck('complete gauge base degree bound %d %d' % (i, j),
                        self.degree(base_numerator[i, j]) <= 4)
            records.append(dict(direction=i,
                fiber_histograms=[self.hist(numerator[i, j]) for j in range(9)],
                base_histograms=[self.hist(base_numerator[i, j]) for j in range(2)]))
        dp = [sum((numerator[i, j] * Delta.derivative(xx[j]) for j in range(6)), R.zero())
              for i in range(3)]
        db = [sum((base_numerator[i, j] * self.bd(R, Delta, p)
                   for j, p in enumerate((self.h, self.z))), R.zero()) for i in range(3)]
        total = [dp[i] + db[i] for i in range(3)]
        for i in range(3):
            self.ck('full Delta total derivative degree bound %d' % i, self.degree(total[i]) <= 7)
            self.ck('Delta independent of RR variables %d' % i, Delta.derivative(xx[6 + i]) == 0)
        payload = dict(ring=R, coefficient_field=self.K, parameters=(self.a, self.b, self.h, self.z),
            variables=xx, base_parameters=(self.h, self.z), new_parameter_images=(self.a, self.r, self.H, self.Z),
            source_Q0=d['Q0'], six_full_sources=src, left_inverse_numerator=W, Delta=Delta,
            source_combination_numerators=WR, source_unit_coefficients=self.units,
            source_permutation=self.sigma, first_gauged_fields=hat, base_fields=base,
            second_gauged_field_numerators=numerator, second_gauged_base_numerators=base_numerator,
            principal_ordinary_field_numerators=matrix(R, 3, 9, [self.layer(R, f, 5) for f in numerator.list()]),
            principal_ordinary_base_numerators=matrix(R, 3, 2, [self.layer(R, f, 4) for f in base_numerator.list()]),
            Delta_total_numerator_derivatives=total, Delta_p_derivative_parts=dp,
            Delta_base_derivative_parts=db, records=records,
            source_files=[str(p) for p in self.paths],
            scope='Complete Delta-rescaled gauge, only necessary source multiples; no detC division.',
            quotient_rule='Fi(f/Delta^k)=(Delta*Pi(f)-k*f*Pi(Delta))/Delta^(k+2); '
                          'Pi includes both h,z coefficient derivatives.')
        return payload, {'first_fields_shape': [3, 9], 'base_fields_shape': [3, 2],
                         'second_fields_shape': [3, 9], 'left_inverse_shape': [3, 6],
                         'six_sources': 6, 'variable_names': list(R.variable_names())}

    def support_view(self):
        # The original support constructor uses exactly the second-gauge ring.
        R = PolynomialRing(self.K, names=('sg_p12', 'sg_p13', 'sg_p21', 'sg_p23',
                                         'sg_p31', 'sg_p32', 'sg_rr1', 'sg_rr2', 'sg_rr3'))
        d = self.data(R)
        xx = R.gens()
        Zvec = d['Zvec']
        J = matrix(R, 3, 9, [f.derivative(v) for f in Zvec for v in xx])
        Jbase = matrix(R, 3, 2, [self.bd(R, f, p) for f in Zvec for p in (self.h, self.z)])
        choices = (((0, 4, 1), (1, 5, -1)), ((2, 2, 1), (5, 3, -1)),
                   ((3, 0, 1), (4, 1, 1)))
        src = vector(R, d['sources'])
        items, records = [], []
        for chart, sel in enumerate(itertools.product(*choices)):
            m = R.one()
            for _, var, _ in sel:
                m *= xx[var]
            W = matrix(R, 3, 6)
            for col, (source, var, sign) in enumerate(sel):
                # Exact polynomial division by a selected coordinate, with every
                # exponent preserved. No fraction-field conversion is needed.
                coefficients = {}
                for exponent, coefficient in m.dict().items():
                    exponent = list(int(v) for v in exponent)
                    self.ck('selected coordinate divides every monomial %d %d' % (chart, col),
                            exponent[var] > 0)
                    exponent[var] -= 1
                    coefficients[tuple(exponent)] = coefficient
                quotient = R(coefficients)
                self.ck('exact selected monomial division %d %d' % (chart, col),
                        quotient * xx[var] == m)
                W[col, source] = sign * quotient
            self.ck('whole support inverse %d' % chart, W * d['D'] == m * identity_matrix(R, 3))
            WR = W * src
            full = m * d['hat']
            for i in range(3):
                full[i, 6 + i] -= self.units[i] * WR[self.sigma[i]]
            basefull = m * d['base']
            for i, f in enumerate(full.list()):
                self.ck('full support degree bound %d %d' % (chart, i), self.degree(f) <= 5)
            for i, f in enumerate(basefull.list()):
                self.ck('full support base degree bound %d %d' % (chart, i), self.degree(f) <= 4)
            self.ck('whole support source numerator top %d' % chart,
                    vector(R, [self.layer(R, f, 6) for f in WR]) == m * Zvec)
            principal = matrix(R, 3, 9, [self.layer(R, f, 5) for f in full.list()])
            basepr = matrix(R, 3, 2, [self.layer(R, f, 4) for f in basefull.list()])
            rec = dict(chart=chart, chosen_rows=[s[0] for s in sel],
                       nonzero_p_indices=[s[1] for s in sel], denominator=str(m),
                       principal_terms=[len(f.dict()) for f in principal.list()],
                       seconds=time.monotonic() - self.start)
            items.append(dict(chart=chart, selection=sel, denominator=m,
                left_inverse_numerator=W, source_numerators=WR,
                full_field_numerators=full, full_base_numerators=basefull,
                principal_field_numerators=principal, principal_base_numerators=basepr, record=rec))
            records.append(rec)
        self.ck('all eight distinct support selections', len(items) == 8 and
                len({tuple(it['selection']) for it in items}) == 8)
        payload = dict(ring=R, coefficient_field=self.K, Zvec=Zvec,
            Z_state_jacobian=J, Z_base_jacobian=Jbase, D=d['D'], items=items,
            tangency_DAG='T_ij=Pi*Z_state_jacobian_j+Bi*Z_base_jacobian_j',
            necessity='On a chosen support chart full leading sources give Z=0; '
                      'Pi((W R)j) has leading Pi(m) Zj+m Tij. '
                      'Coefficient h,z derivatives are taken before specialization.',
            coverage='Exactly rank D=3: choose one nonzero variable in each '
                     '(p31,p32),(p21,p23),(p12,p13). No detL,detC or RR division.',
            scope='All eight complete support source gauges, including detL=0; inactive columns excluded.')
        return payload, {'charts': 8, 'D_shape': [6, 3], 'W_shape_each_chart': [3, 6],
                         'fields_shape_each_chart': [3, 9], 'base_shape_each_chart': [3, 2],
                         'Z_state_jacobian_shape': [3, 9], 'Z_base_jacobian_shape': [3, 2],
                         'selections_in_original_product_order': [it['selection'] for it in items],
                         'variable_names': list(R.variable_names())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--profiles', nargs='+', choices=tuple(NAMES), default=list(NAMES))
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    root, output = args.bundle_root.absolute(), args.output_root.absolute()
    output.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    views = FreshViews(root)
    result = dict(schema='ericksen_full_fresh_source_views/v1', producer_dependency=PRODUCER_ID,
                  input_receipts=views.inputs, mathematical_source=file_sha(Path(__file__)),
                  source_views_are_new_adapters=True, inherited_historical_arrays_used=False,
                  profiles={})
    for profile in args.profiles:
        before = len(views.checks)
        payload, dimensions = getattr(views, profile + '_view')()
        views.ck('complete consumed schema ' + profile,
                 all(k in payload for k in REQUIRED_KEYS[profile]))
        payload['checks'] = (len(views.checks[before:]) if profile == 'support'
                             else list(views.checks[before:]))
        payload['provenance'] = dict(producer_dependency=PRODUCER_ID, input_receipts=views.inputs,
            full_polynomial_support_retained=True, profiles=profile,
            construction_reference=CONSTRUCTION_REFERENCES[profile],
            construction_source_sha256=CONSTRUCTION_SOURCE_SHA256,
            historical_array_comparison_performed=False)
        target = output / (NAMES[profile] + '.sobj')
        save(payload, str(target), compress=False)
        report = dict(status='PASS_FULL_FRESH_SOURCE_SCHEMA_VIEW', profile=profile,
            full_schema_keys=sorted(payload), required_consumer_keys=list(REQUIRED_KEYS[profile]),
            dimensions=dimensions, checks=len(views.checks[before:]),
            check_labels=list(views.checks[before:]),
            input_receipts=views.inputs, output_sobj_sha256=file_sha(target),
            producer_dependency=PRODUCER_ID, construction_reference=CONSTRUCTION_REFERENCES[profile],
            construction_source_sha256=CONSTRUCTION_SOURCE_SHA256,
            source_views_are_new_adapters=True, historical_array_comparison_performed=False)
        if profile == 'global':
            report.update(full_sources=18, degree_bound=5,
                full_DAG_product_counts=[len(v['full_signed_product_DAG']) for v in payload['records']],
                top_DAG_product_counts=[len(v['ordinary_degree5_signed_product_DAG']) for v in payload['records']])
        elif profile == 'second':
            report.update(records=payload['records'],
                Delta_derivative_terms=[len(f.dict()) for f in payload['Delta_total_numerator_derivatives']],
                Delta_base_part_terms=[len(f.dict()) for f in payload['Delta_base_derivative_parts']])
        else:
            report['records'] = [item['record'] for item in payload['items']]
        report.update(seconds=time.monotonic() - start, scope=payload['scope'])
        report_path = output / (NAMES[profile] + '.json')
        report_path.write_text(json.dumps(report, indent=2) + '\n')
        result['profiles'][profile] = dict(sobj=str(target), sobj_sha256=file_sha(target),
            report=str(report_path), report_sha256=file_sha(report_path), dimensions=dimensions)
        print('COMPLETE FRESH VIEW', profile, target, flush=True)
    result.update(status='PASS_FULL_FRESH_SOURCE_SCHEMA_VIEWS', seconds=time.monotonic() - start,
                  initial_contract_checks=views.checks[:7])
    record = args.record or output / 'source_views_generation.json'
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(json.dumps(result, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
