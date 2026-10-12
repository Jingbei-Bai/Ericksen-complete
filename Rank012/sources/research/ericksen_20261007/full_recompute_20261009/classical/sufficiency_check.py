"""Exact, all-parameter checks for the six classical-family sufficiency proof.

Run with an ordinary SymPy Python.  No numerical material, spectrum, or
parameter specialization is used as evidence.  The positivity and nonvanishing
assumptions identify regular coordinate charts; signed determinants retain the
orientation restrictions in paper/main.tex.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).absolute().parent
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--main-tex', type=Path)
args = parser.parse_args()
INTRO = args.main_tex.resolve() if args.main_tex else next(
    (base / 'paper/main.tex' for base in HERE.parents
     if (base / 'paper/main.tex').is_file()), None)
assert INTRO is not None, 'Keep paper/main.tex in the proof bundle.'
checks = []


def check(name, expression):
    if isinstance(expression, s.MatrixBase):
        residual = expression.applyfunc(s.simplify)
        ok = all(x == 0 for x in residual)
    else:
        residual = s.simplify(expression)
        ok = residual == 0
    checks.append({"name": name, "passed": bool(ok), "residual": str(residual)})
    if not ok:
        raise AssertionError((name, residual))


r, R, r0 = s.symbols("r R r_star", positive=True)
a, b, c, A, C, D, E, K = s.symbols("a b c A C D E K", real=True)
B0 = s.symbols("B_0", real=True)
alpha = s.symbols("alpha", positive=True)
beta, k, gamma = s.symbols("beta k gamma", real=True)
epsilon = s.symbols("epsilon", real=True)
u, v = s.symbols("W_1 W_2", real=True)
X, Y, Z, Theta, theta, x, y, z = s.symbols("X Y Z Theta theta x y z", real=True)
d, x0 = s.symbols("d x_0", real=True)

# Coordinate Jacobians -> physical gradients in oriented orthonormal bases.
F1 = s.Matrix([[a/r, 0, 0], [0, b*r, 0], [0, c, 1/(a*b)]])
r_of_X = s.sqrt(2*a*X+d)
map1 = s.Matrix([r_of_X, b*Y, Z/(a*b)+c*Y])
physical1 = s.diag(1, r_of_X, 1) * map1.jacobian([X, Y, Z])
check("F1: chart differentiation", physical1 - F1.subs(r, r_of_X))
check("F1: signed Jacobian 1", F1.det()-1)

F2 = s.Matrix([[a*R, 0, 0], [0, b/R, 0], [0, c/R, 1/(a*b)]])
map2 = s.Matrix([a*R**2/2+x0, b*Theta, Z/(a*b)+c*Theta])
physical2 = map2.jacobian([R, Theta, Z]) * s.diag(1, 1/R, 1)
check("F2: chart differentiation", physical2-F2)
check("F2: signed Jacobian 1", F2.det()-1)

F3 = s.Matrix([[A*R/r, 0, 0], [0, C*r/R, D*r], [0, E/R, K]])
r_of_R = s.sqrt(A*R**2+B0)
map3 = s.Matrix([r_of_R, C*Theta+D*Z, E*Theta+K*Z])
physical3 = s.diag(1, r_of_R, 1)*map3.jacobian([R, Theta, Z])*s.diag(1, 1/R, 1)
check("F3: chart differentiation retaining D and E", physical3-F3.subs(r, r_of_R))
check("F3: signed Jacobian A(CK-DE)", F3.det()-A*(C*K-D*E))
B3_expected = s.Matrix([
    [A**2*R**2/r**2, 0, 0],
    [0, r**2*(C**2/R**2+D**2), r*(C*E/R**2+D*K)],
    [0, r*(C*E/R**2+D*K), E**2/R**2+K**2],
])
check("F3: full B block", F3*F3.T-B3_expected)

for label, F in [("F1", F1), ("F2", F2), ("F3", F3)]:
    B = F*F.T
    Bi = B.inv()
    for tensor_label, tensor in [("B", B), ("B inverse", Bi)]:
        check(label+": "+tensor_label+" first-to-second component", tensor[0, 1])
        check(label+": "+tensor_label+" first-to-third component", tensor[0, 2])
        # All entries are even rational functions of R (F1 does not use R),
        # hence functions of R^2, the current radial/Cartesian coordinate.
        check(label+": "+tensor_label+" even in R", tensor.subs(R, -R)-tensor)
    check(label+": inverse product", B*Bi-s.eye(3))

# F2: arbitrary stress coefficients dependent on x still close.
Sxx, Syy, Szz, Syz = [s.Function(name)(x) for name in ["Sxx", "Syy", "Szz", "Syz"]]
S2 = s.Matrix([[Sxx, 0, 0], [0, Syy, Syz], [0, Syz, Szz]])
div2 = s.Matrix([sum(s.diff(S2[i, j], q) for j, q in enumerate([x,y,z])) for i in range(3)])
check("F2: arbitrary Cartesian first-coordinate stress divergence",
      div2-s.Matrix([s.diff(Sxx,x), 0, 0]))
check("F2: Bxx and inverse Bxx in explicit pressure",
      s.Matrix([(F2*F2.T)[0,0]-a**2*R**2,
                (F2*F2.T).inv()[0,0]-1/(a**2*R**2)]))

# Derive cylindrical divergence from Euclidean derivatives in a rotating frame.
# The covariant derivative operators on tensor components are:
# D_r S=S', D_theta S=(J S-S J)/r, D_z S=0.
Sr, St, Sz, Srt, Srz, Stz = [
    s.Function(name)(r) for name in ["Srr","Stt","Szz","Srt","Srz","Stz"]
]
Scyl = s.Matrix([[Sr,Srt,Srz], [Srt,St,Stz], [Srz,Stz,Sz]])
J = s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
div_cyl_derived = s.diff(Scyl,r)[:,0] + ((J*Scyl-Scyl*J)/r)[:,1]
div_cyl_expected = s.Matrix([
    s.diff(Sr,r)+(Sr-St)/r, s.diff(Srt,r)+2*Srt/r, s.diff(Srz,r)+Srz/r
])
check("General cylindrical divergence derived using rotating orthonormal frame",
      div_cyl_derived-div_cyl_expected)
t = s.symbols("t", positive=True)
pcyl = Sr+s.Integral((Sr.subs(r,t)-St.subs(r,t))/t, (t,r0,r))
check("General radial-block pressure, arbitrary coefficient functions",
      s.diff(pcyl,r)-(s.diff(Sr,r)+(Sr-St)/r))

# F4 keeps det Q=epsilon, including epsilon=-1.  Reduce epsilon^2=1.
rprime = epsilon*R**2/r**2
signed4 = epsilon*rprime*r**2/R**2
check("F4: signed Jacobian modulo epsilon^2=1", s.rem(s.expand(signed4-1), epsilon**2-1, epsilon))
br, bt = R**4/r**4, r**2/R**2
check("F4: radial squared stretch modulo epsilon^2=1",
      s.rem(s.expand(rprime**2-br), epsilon**2-1, epsilon))
check("F4: strain determinant 1", br*bt**2-1)
# Euclidean spherical identity: div(f I+g n tensor n)=
# (f'+g'+2g/r)n; n dot grad n=0 and div n=2/r.
xx = s.Matrix([x,y,z])
rad = s.sqrt(xx.dot(xx))
nn = xx/rad
check("Spherical radial direction has divergence 2/r",
      sum(s.diff(nn[i], xx[i]) for i in range(3))-2/rad)
check("Spherical radial direction is constant along radial rays",
      nn.jacobian(xx)*nn)
psph = Sr+2*s.Integral((Sr.subs(r,t)-St.subs(r,t))/t, (t,r0,r))
check("F4: spherical pressure, arbitrary coefficient functions",
      s.diff(psph,r)-(s.diff(Sr,r)+2*(Sr-St)/r))

# F5: no assumption of a particular spectrum, energy, or shear.
F5 = s.Matrix([[alpha,0,0], [alpha*beta,alpha*k,0], [0,0,gamma]])
map5 = s.Matrix([alpha*R, beta*s.log(R)+k*Theta, gamma*Z])
physical5 = s.diag(1,alpha*R,1)*map5.jacobian([R,Theta,Z])*s.diag(1,1/R,1)
check("F5: chart differentiation", physical5-F5)
check("F5: signed Jacobian alpha^2 k gamma", F5.det()-alpha**2*k*gamma)
B5 = F5*F5.T
B5i = s.Matrix([
    [(beta**2+k**2)/(alpha**2*k**2), -beta/(alpha**2*k**2), 0],
    [-beta/(alpha**2*k**2), 1/(alpha**2*k**2), 0],
    [0,0,gamma**-2]
])
check("F5: complete inverse", B5*B5i-s.eye(3))
check("F5: I1 constant formula", s.trace(B5)-(alpha**2*(1+beta**2+k**2)+gamma**2))
check("F5: I2 constant formula", s.trace(B5i)-((1+beta**2+k**2)/(alpha**2*k**2)+gamma**-2))
M = alpha**2*u+v/(alpha**2*k**2)
S5 = 2*u*B5-2*v*B5i
check("F5: radial pressure coefficient", S5[0,0]-S5[1,1]-2*M*(1-beta**2-k**2))
check("F5: angular pressure coefficient", S5[0,1]-2*beta*M)
div5 = s.Matrix([(S5[0,0]-S5[1,1])/r, 2*S5[0,1]/r, 0])
p5 = 2*M*((1-beta**2-k**2)*s.log(r/r0)+2*beta*theta)
gradient5 = s.Matrix([s.diff(p5,r), s.diff(p5,theta)/r, s.diff(p5,z)])
check("F5: full local pressure gradient for arbitrary W1 W2", div5-gradient5)
check("F5: angular part has zero local curl",
      s.diff(r*div5[1],r)/r)

# The displayed removable F3 axis: C=+/-1, K=1/(AC), B0=E=0.
# Use C^2=1 only; both orientation choices are retained.
Faxis = s.Matrix([[s.sqrt(A),0,0], [0,C*s.sqrt(A),D*r], [0,0,1/(A*C)]])
Baxis = Faxis*Faxis.T
def sign_reduce(expr):
    num, den = s.fraction(s.cancel(expr))
    return s.rem(num,C**2-1,C)
for name, expr in [
    ("F3 axis Brr", Baxis[0,0]-A),
    ("F3 axis Btt", Baxis[1,1]-A-D**2*r**2),
    ("F3 axis inverse Brr", Baxis.inv()[0,0]-1/A),
    ("F3 axis inverse Btt", Baxis.inv()[1,1]-1/A),
    ("F3 axis I1", s.trace(Baxis)-(2*A+A**-2+D**2*r**2)),
    ("F3 axis I2", s.trace(Baxis.inv())-(2/A+A**2+A*D**2*r**2)),
]:
    check(name+" modulo C^2=1", sign_reduce(expr))

result = {
    "status": "PASS",
    "method": "exact SymPy rational identities and formal derivatives",
    "checks_count": len(checks),
    "checks": checks,
    "source": str(INTRO),
    "source_sha256": hashlib.sha256(INTRO.read_bytes()).hexdigest(),
    "scope": {
        "F0": "Any L in SL(3,R): constant Cartesian stress, constant pressure.",
        "F1_F2": "All real a,b,c with ab != 0 on the stated positive-radius charts.",
        "F3": "All real A,C,D,E,K with A(CK-DE)=1, arbitrary B0, R,r>0; no loss of D or E and no A>0 assumption on annuli.",
        "F4": "epsilon=det Q=+/-1; R,r>0, r^3=epsilon R^3+B0.",
        "F5": "alpha>0, real beta,k,gamma with alpha^2*k*gamma=1; no simple-spectrum restriction.",
        "material": "W arbitrary and sufficiently smooth; radial-block checks allow arbitrary radial W1,W2, and F5 uses arbitrary constants W1,W2 at its invariant pair.",
        "regular_continuations": "Closedness extends across coordinate singularities by continuity for C3 regular maps; local pressure follows from the local Poincare lemma.",
        "global_pressure": "Not asserted: F5 pressure generally contains an angular term."
    },
    "formulas": {
        "signed_jacobians": ["det L=1", "1", "1", "A(CK-DE)=1", "(det Q)*r'*(r/R)^2=epsilon^2=1", "alpha^2*k*gamma=1"],
        "F5_pressure": "2*(alpha^2*W1+W2/(alpha^2*k^2))*((1-beta^2-k^2)*log(r/r_star)+2*beta*(theta-theta_star))"
    }
}
(HERE/"sufficiency_check.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
print(json.dumps({"status": result["status"],"checks_count": len(checks),"output": str(HERE/"sufficiency_check.json")}))
