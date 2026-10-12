"""Required actual wall first-variation source, without abandoned projection tail.

This is exactly the raw coefficient stage of
independent_inactive_RRtorus_wall_linear_variation.py. The original tail has
no stored result and no dependency in the maintained proof.
"""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).absolute().parent;start=time.monotonic()
d=load(str(p/'independent_inactive_RRtorus_normal_block_reader.sobj'))
raw=load(str(p/'independent_inactive_RRtorus_infinity_normal.sobj'))
Q=d['ring'];l,rho=Q.gens();K=Q.base_ring();a,r,H,Z=K.gens()
m=1+a*r;q=1+r;F=Q.fraction_field();old=raw['coefficient_field']
ph=old.ring().hom([a,a*q/m,H/r,q*Z/r],K);cache={}
def cf(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def cvp(f):
 if not hasattr(f,'dict'):return Q(cf(f))
 return Q({tuple(map(int,ex[:2])):cf(c) for ex,c in f.dict().items()})
def cv(f):return F(cvp(f.numerator()))/cvp(f.denominator())
ids=[0,1,2,3,5,7,11,13,17,19,23,25]
N=matrix(F,[[cv(raw['boundary_matrix'][i,j]) for j in range(6)] for i in ids])
kap=raw['ring'].gens()[2]
def atzero(f):
 numerator=f.numerator().subs({kap:0});denominator=f.denominator().subs({kap:0})
 assert denominator,'Actual transported derivative must be regular'
 return F(cvp(numerator))/cvp(denominator)
Nk=matrix(F,[[atzero(raw['transported_matrix'][i,j].derivative(kap)) for j in range(6)] for i in ids])
save({'ring':Q,'source_rows':ids,'normal':N,'kappa_derivative':Nk},str(p/'independent_inactive_RRtorus_wall_variation_raw.sobj'),compress=False)
print('ACTUAL RAW WALL SOURCE REBUILT',N.dimensions(),Nk.dimensions(),time.monotonic()-start,flush=True)
