"""A complete, exact F(e) backend with finite polynomial expression DAGs.

This module never takes a jet, quotient, tolerance, or precision argument.
Every polynomial leaf retains all its coefficients. Polynomial addition and
multiplication retain complete expressions. Valuation/zero queries inspect
exact F coefficients up to a finite degree upper bound, stopping only at the
first nonzero coefficient or after all possible coefficients have vanished.

The interface is deliberately limited to the original reader324 operations;
it is not a general replacement for every Sage fraction-field method.
"""
from sage.all import *
from sage.rings.fraction_field import FractionField_generic
from sage.structure.element import FieldElement, parent
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.structure.richcmp import op_EQ, op_NE
from collections import Counter


class FullPolynomial:
    _intern={}

    def __new__(cls, field, kind, children=(), coefficients=None, shift=0):
        if kind=='leaf':
            data=tuple(sorted((int(k),field(c)) for k,c in coefficients.items() if c))
            key=(id(field),kind,data)
        else:
            data=None
            if kind in ('add','mul'):
                left,right=children
                if kind=='add':
                    if left.upper<0:return right
                    if right.upper<0:return left
                    if (left.kind=='neg' and left.children[0] is right) or (right.kind=='neg' and right.children[0] is left):
                        return cls(field,'leaf',coefficients={})
                else:
                    if left.upper<0 or right.upper<0:return cls(field,'leaf',coefficients={})
                    if left.kind=='leaf' and left.upper==0 and left.coefficient(0)==1:return right
                    if right.kind=='leaf' and right.upper==0 and right.coefficient(0)==1:return left
                if left.kind=='leaf' and left.upper==0 and right.kind=='leaf' and right.upper==0:
                    value=left.coefficient(0)+right.coefficient(0) if kind=='add' else left.coefficient(0)*right.coefficient(0)
                    return cls(field,'leaf',coefficients={0:value})
                children=tuple(sorted(children,key=id))
            elif kind=='neg':
                child=children[0]
                if child.kind=='neg':return child.children[0]
                if child.kind=='leaf':return cls(field,'leaf',coefficients={k:-c for k,c in child.cache.items()})
            elif kind=='shift':
                child=children[0]
                if not shift:return child
                assert shift>=0 and child.lower>=shift
                if child.kind=='leaf':return cls(field,'leaf',coefficients={k-shift:c for k,c in child.cache.items()})
            key=(id(field),kind,tuple(id(c) for c in children),int(shift))
        if key in cls._intern:return cls._intern[key]
        node=object.__new__(cls)
        node._initial_spec=(field,kind,children,data,shift)
        cls._intern[key]=node
        return node

    def __init__(self, field, kind, children=(), coefficients=None, shift=0):
        if hasattr(self,'initialized'):return
        field,kind,children,data,shift=self._initial_spec
        self.initialized=True
        self.field, self.kind, self.children = field, kind, children
        self.shift = shift
        self.cache = {}
        self.cached_valuation = None
        self.cached_factors=None
        if kind == 'leaf':
            self.cache = dict(data)
            self.lower = min(self.cache, default=0)
            self.upper = max(self.cache, default=-1)
        elif kind == 'add':
            self.lower = min(c.lower for c in children)
            self.upper = max(c.upper for c in children)
        elif kind == 'mul':
            self.lower = sum(c.lower for c in children)
            self.upper = sum(c.upper for c in children)
        elif kind == 'neg':
            self.lower, self.upper = children[0].lower, children[0].upper
        elif kind == 'shift':
            assert shift >= 0 and children[0].lower >= shift
            self.lower, self.upper = children[0].lower-shift, children[0].upper-shift
        else:
            raise ValueError(kind)

    def coefficient(self, k):
        k = int(k)
        if k < self.lower or k > self.upper:
            return self.field.zero()
        if k in self.cache:
            return self.cache[k]
        if self.kind == 'leaf':
            value = self.field.zero()
        elif self.kind == 'add':
            value = self.children[0].coefficient(k)+self.children[1].coefficient(k)
        elif self.kind == 'neg':
            value = -self.children[0].coefficient(k)
        elif self.kind == 'shift':
            value = self.children[0].coefficient(k+self.shift)
        else:
            left, right = self.children
            first = max(left.lower, k-right.upper)
            last = min(left.upper, k-right.lower)
            value = sum((left.coefficient(j)*right.coefficient(k-j)
                         for j in range(first, last+1)), self.field.zero())
        self.cache[k] = value
        return value

    def valuation(self):
        if self.cached_valuation is None:
            value = Infinity
            for k in range(self.lower, self.upper+1):
                if self.coefficient(k):
                    value = Integer(k)
                    break
            # Commit only after all queries that determine the answer succeed.
            self.cached_valuation = value
        return self.cached_valuation

    def is_zero(self):
        return self.valuation() == Infinity

    def native(self, ring):
        """Full expansion, used only by small implementation contracts."""
        return ring({k:self.coefficient(k) for k in range(self.lower,self.upper+1)
                     if self.coefficient(k)})

    def __repr__(self):
        return 'FullPolynomial(%s, degree_bound=%s)' % (self.kind,self.upper)


def add(left, right):
    if left.upper < 0:
        return right
    if right.upper < 0:
        return left
    return FullPolynomial(left.field,'add',(left,right))


def multiply(left, right):
    if left.upper < 0 or right.upper < 0:
        return FullPolynomial(left.field,'leaf',coefficients={})
    return FullPolynomial(left.field,'mul',(left,right))


def factor_data(polynomial):
    """Only explicit product factors; no polynomial factorization inference."""
    if polynomial.cached_factors is None:
        F=polynomial.field
        if polynomial.kind=='mul':
            c1,f1=factor_data(polynomial.children[0]);c2,f2=factor_data(polynomial.children[1])
            scalar=c1*c2;factors=Counter(dict(f1))+Counter(dict(f2))
        elif polynomial.kind=='neg':
            scalar,stored=factor_data(polynomial.children[0]);scalar=-scalar;factors=Counter(dict(stored))
        elif polynomial.kind=='leaf' and polynomial.upper<=0:
            scalar=polynomial.coefficient(0);factors=Counter()
        else:
            scalar=F.one();factors=Counter({polynomial:1})
        polynomial.cached_factors=(scalar,tuple(factors.items()))
    scalar,stored=polynomial.cached_factors
    return scalar,Counter(dict(stored))


def factor_product(F,scalar,factors):
    result=FullPolynomial(F,'leaf',coefficients={0:scalar})
    for factor,count in sorted(factors.items(),key=lambda item:id(item[0])):
        for unused in range(count):result=multiply(result,factor)
    return result


def fraction_sum_numerator(left,right,subtract=False):
    c1,d1=factor_data(left.d);c2,d2=factor_data(right.d)
    common=d1 & d2
    q1=factor_product(left.parent().coefficient_field,c1,d1-common)
    q2=factor_product(left.parent().coefficient_field,c2,d2-common)
    second=multiply(right.n,q1)
    if subtract:second=FullPolynomial(second.field,'neg',(second,))
    numerator=add(multiply(left.n,q2),second)
    denominator=multiply(left.d,q2)
    return numerator,denominator


class CompleteFraction(FieldElement):
    def __init__(self, domain, numerator, denominator=None, **unused):
        FieldElement.__init__(self,domain)
        n = domain.polynomial(numerator)
        d = domain.polynomial(1 if denominator is None else denominator)
        if d.is_zero():
            raise ZeroDivisionError('Zero complete polynomial denominator')
        # These cancellations are solely between structurally identical,
        # interned complete factors. The exact denominator is nonzero.
        cn,nf=factor_data(n);cd,df=factor_data(d);shared=nf & df
        n=factor_product(domain.coefficient_field,cn/cd,nf-shared)
        d=factor_product(domain.coefficient_field,domain.coefficient_field.one(),df-shared)
        # This only removes a guaranteed structural common monomial. All
        # coefficients and every remaining degree are retained.
        common = min(n.lower,d.lower) if n.upper >= 0 else 0
        if common:
            n=FullPolynomial(domain.coefficient_field,'shift',(n,),shift=common)
            d=FullPolynomial(domain.coefficient_field,'shift',(d,),shift=common)
        self.n,self.d=n,d

    def numerator(self):
        return self.n

    def denominator(self):
        return self.d

    def __bool__(self):
        return not self.n.is_zero()

    def _add_(self, other):
        other=self.parent()(other)
        numerator,denominator=fraction_sum_numerator(self,other)
        return CompleteFraction(self.parent(),numerator,denominator)

    def _neg_(self):
        return CompleteFraction(self.parent(),
             FullPolynomial(self.parent().coefficient_field,'neg',(self.n,)),self.d)

    def _mul_(self, other):
        other=self.parent()(other)
        return CompleteFraction(self.parent(),multiply(self.n,other.n),multiply(self.d,other.d))

    def _div_(self, other):
        other=self.parent()(other)
        if not other:
            raise ZeroDivisionError('Division by zero exact rational function')
        return CompleteFraction(self.parent(),multiply(self.n,other.d),multiply(self.d,other.n))

    def __invert__(self):
        if not self:
            raise ZeroDivisionError('Inverse of zero exact rational function')
        return CompleteFraction(self.parent(),self.d,self.n)

    def _richcmp_(self, other, operation):
        other=self.parent()(other)
        if operation not in (op_EQ,op_NE):
            raise TypeError('This rational-function field is not ordered')
        difference,unused=fraction_sum_numerator(self,other,subtract=True)
        equal=difference.is_zero()
        return equal if operation==op_EQ else not equal

    def _repr_(self):
        return 'CompleteFraction(%r / %r)' % (self.n,self.d)

    def native(self):
        E=self.parent().old_ring
        return self.parent().old_field(self.n.native(E))/self.parent().old_field(self.d.native(E))


class CompleteFractionField(FractionField_generic):
    def __init__(self, E, old_field):
        self.old_ring=E
        self.old_field=old_field
        self.coefficient_field=E.base_ring()
        F=self.coefficient_field
        # A distinct carrier ring prevents Sage's generic fraction-parent
        # comparison from identifying this representation with old_field.
        # Its generator is identified exactly with the old e; no new algebraic
        # variable is introduced into values.
        self.carrier_ring=PolynomialRing(F,'e_complete_reader324_v2')
        FractionField_generic.__init__(self,self.carrier_ring,element_class=CompleteFraction)
        for source in (F,E,old_field,self.carrier_ring):
            self.register_coercion(SetMorphism(Hom(source,self,Rings()),self._element_constructor_))

    def polynomial(self, value):
        if isinstance(value,FullPolynomial):
            assert value.field is self.coefficient_field
            return value
        z=self.old_ring(value.dict()) if parent(value) is self.carrier_ring else self.old_ring(value)
        return FullPolynomial(self.coefficient_field,'leaf',coefficients=z.dict())

    def _element_constructor_(self, value, denominator=None, **unused):
        if isinstance(value,CompleteFraction) and value.parent() is self:
            return value if denominator is None else value/self(denominator)
        if denominator is not None:
            return self(value)/self(denominator)
        if parent(value) is self.old_field:
            return CompleteFraction(self,value.numerator(),value.denominator())
        return CompleteFraction(self,value)

    def _repr_(self):
        return 'Complete exact DAG representation of %s' % self.old_field

    def gen(self, i=0):
        assert i==0
        return self(self.old_ring.gen())


def small_contract():
    K=NumberField(polygen(QQ)**2+1,'ii')
    P=PolynomialRing(K,names=('s','t'));F=P.fraction_field();s,t=F.gens()
    E=PolynomialRing(F,'e');e=E.gen();FE=E.fraction_field();L=CompleteFractionField(E,FE)
    a=FE((s+t)*e**7+(s-t)*e**9)/FE((s*t+1)*e**2+e**5)
    b=FE(e**4+(s/t)*e**13)/FE(1+s*e**3)
    A,B=L(a),L(b)
    tests={}
    for key,new,old in [('sum',A+B,a+b),('difference',A-B,a-b),
                        ('product',A*B,a*b),('quotient',A/B,a/b),
                        ('inverse',~A,1/a),('positive_power',A**3,a**3),
                        ('negative_power',A**(-2),a**(-2))]:
        tests[key]=new.native()==old
    tests['valuation_complete_high_degree']=A.n.valuation()-A.d.valuation()==a.numerator().valuation()-a.denominator().valuation()
    tests['exact_zero_full_bound']=not (A-A)
    tests['nonzero_above_six']=bool(L(e**37)) and L(e**37).n.valuation()==37
    tests['nonzero_beyond_any_fixed_contract_order']=bool(L(e**123)) and L(e**123).n.valuation()==123
    tests['zero_infinity']=L.zero().n.valuation()==Infinity
    tests['complex_coefficient_retained']=L(F(K.gen())).native()==FE(K.gen())
    S=PolynomialRing(L,names=('a','b','l'));aa,bb,ll=S.gens()
    R=PolynomialRing(F,names=('a','b'));ra,rb=R.gens();ev=R.hom([aa,bb],S)
    tests['polynomial_hom']=ev(ra*rb*s)==aa*bb*L(s)
    tests['polynomial_arithmetic']=(ll-A*aa-B*bb)**2==ll**2+A*A*aa**2+B*B*bb**2-2*A*aa*ll-2*B*bb*ll+2*A*B*aa*bb
    M=matrix(FE,[[a,b],[1/a,a+b]])
    T=diagonal_matrix(L,[e**3,FE(e**(-2))])*M
    tests['matrix_parent']=T.base_ring() is L
    tests['matrix_entry']=T[1,1].native()==FE(e**(-2))*(a+b)
    tests['field_equality']=A==L(a) and A!=B
    tests['common_denominator_multiplicity']=(A+A).d is A.d
    tests['inverse_explicit_factor_cancellation']=(A*(~A)).n.upper==0 and (A*(~A)).d.upper==0
    tests['nested_rational_arithmetic']=((A+B+A-B)/(A+1)).native()==(a+b+a-b)/(a+1)
    tests['complete_carrier_inclusion']=L(L.carrier_ring.gen()).native()==FE(e)
    assert all(tests.values()),tests
    return {'status':'PASS_SMALL_COMPLETE_DAG_CONTRACT','checks':tests,'count':len(tests),
            'full_reader_executed':False,'series_truncation_used':False,
            'coefficient_queries_have_no_precision_argument':True}


if __name__=='__main__':
    import argparse,hashlib,json,time
    from pathlib import Path
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    options=parser.parse_args()
    st=time.time();report=small_contract();report['seconds']=time.time()-st
    report['backend_module_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if options.output:
        options.output.write_text(json.dumps(report,indent=2),encoding='utf8')
    print(json.dumps(report,indent=2),flush=True)
