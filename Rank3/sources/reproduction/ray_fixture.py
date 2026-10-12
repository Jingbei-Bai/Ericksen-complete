"""Restore exact ray coefficient data and its precision metadata.

The whole Laurent reader computes its own norm. Its comparison input uses only
the leading coefficient and the stored absolute precision. The metadata object
therefore intentionally provides no series arithmetic or coefficient access.
"""
from sage.all import QQ, PolynomialRing

FIXTURE_KIND = 'ericksen_exact_ray_leading_coefficient_and_precision/v1'


class PrecisionMetadata:
    __slots__ = ('_precision',)

    def __init__(self, precision):
        self._precision = int(precision)

    def precision_absolute(self):
        return self._precision


def restore_ray_fixture(payload):
    if not isinstance(payload, dict) or payload.get('fixture_kind') != FIXTURE_KIND:
        return payload
    assert set(payload) == {'fixture_kind', 'variable_count', 'leading_coefficient', 'absolute_precision'}
    assert payload['variable_count'] == 2
    ring = PolynomialRing(QQ, names=('fixture_a', 'fixture_b'))
    field = ring.fraction_field()

    def polynomial(records):
        terms = {}
        for monomial, numerator, denominator in records:
            exponent = tuple(int(e) for e in monomial)
            assert len(exponent) == 2 and all(e >= 0 for e in exponent)
            assert int(denominator) != 0 and exponent not in terms
            coefficient = QQ(int(numerator)) / QQ(int(denominator))
            assert coefficient
            terms[exponent] = coefficient
        return ring(terms)

    coefficient = payload['leading_coefficient']
    assert set(coefficient) == {'numerator', 'denominator'}
    numerator = polynomial(coefficient['numerator'])
    denominator = polynomial(coefficient['denominator'])
    assert denominator
    leading = field(numerator) / field(denominator)
    assert leading.numerator() * denominator == numerator * leading.denominator()
    return {'material_ring': ring,
            'leading_coefficient': leading,
            'source_norm': PrecisionMetadata(payload['absolute_precision'])}
