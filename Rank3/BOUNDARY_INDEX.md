# Rank 3 boundary index

The valuation partition is Lemma `r3:boundarycoverage`. The 73 rows below are the coordinate support patterns allowed by the leading incidence equations in the three intermediate regions. Each row has a separate route for the supercritical, critical and subcritical regions. Coefficient rank and parameter loci further refine these support patterns.

`rank_D` refers to the six-by-three incidence coefficient matrix in Appendix C.1. It does not refer to a Schur matrix or a linear normal matrix.

## Valuation regions

| Condition | Chart | Calculation nodes | Mathematical argument |
|---|---|---|---|
| s_0 = 0 | original weighted chart | BC02, BC10, BC11, BC12, BC13 | Cone normal ideal and the four-point intersection ideal contain tau. |
| 0 < s_0 < n/2 | W_T | BC21, BC30, BC33 | Complete leading coefficient identities; transverse normal equations on full derivative support. |
| s_0 = n/2 | W_S intersect W_T | BC14, BC15, BC30, BC31, BC32, BC33 | Critical amplitude and direction matrices, their rank strata, original-row Schur reductions and residual equations. |
| n/2 < s_0 < n | W_S | BC06, BC16, BC17, BC18, C9_ENDPOINT, BC19, BC20, BC30, BC32, BC33 | Relative local ideals on inactive supports; two horizontal derivations and generic support exclusion on the other supports. |
| n <= s_0 <= infinity | W_O | BC22, BC23, BC24, BC25, BC26, BC27, BC29 | Local radical containment in F_in, followed by algebraicity and differential exclusion of its dominant component. |

## Intermediate supports

| R support | p support | rank D | Canonical support | Supercritical | Critical | Subcritical |
|---|---|---:|---|---|---|---|
| 3 | 23 | 1 | 23 (mask 8) | BC21 | BC14 | BC16 |
| 3 | 21 | 1 | 21 (mask 4) | BC21 | BC14 | BC17 |
| 3 | 21,23 | 1 | 21,23 (mask 12) | BC21 | BC14 | BC16 |
| 3 | 13 | 1 | 13 (mask 2) | BC21 | BC14 | BC16 |
| 3 | 13,23 | 2 | 13,23 (mask 10) | BC21 | BC14 | BC18 |
| 3 | 13,21 | 2 | 13,21 (mask 6) | BC21 | BC14 | BC16 |
| 3 | 13,21,23 | 2 | 13,21,23 (mask 14) | BC21 | BC14 | BC16 |
| 3 | 12 | 1 | 12 (mask 1) | BC21 | BC14 | BC17 |
| 3 | 12,23 | 2 | 12,23 (mask 9) | BC21 | BC14 | BC16 |
| 3 | 12,21 | 2 | 12,21 (mask 5) | BC21 | BC14 | BC16 |
| 3 | 12,21,23 | 2 | 12,21,23 (mask 13) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 3 | 12,13 | 1 | 12,13 (mask 3) | BC21 | BC14 | BC16 |
| 3 | 12,13,23 | 2 | 12,13,23 (mask 11) | BC21 | BC14 | BC16 |
| 3 | 12,13,21 | 2 | 12,13,21 (mask 7) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 3 | 12,13,21,23 | 2 | 12,13,21,23 (mask 15) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 2 | 32 | 1 | 13 (mask 2) | BC21 | BC14 | BC16 |
| 2 | 31 | 1 | 12 (mask 1) | BC21 | BC14 | BC17 |
| 2 | 31,32 | 1 | 12,13 (mask 3) | BC21 | BC14 | BC16 |
| 2 | 13 | 1 | 21 (mask 4) | BC21 | BC14 | BC17 |
| 2 | 13,32 | 2 | 13,21 (mask 6) | BC21 | BC14 | BC16 |
| 2 | 13,31 | 2 | 12,21 (mask 5) | BC21 | BC14 | BC16 |
| 2 | 13,31,32 | 2 | 12,13,21 (mask 7) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 2 | 12 | 1 | 23 (mask 8) | BC21 | BC14 | BC16 |
| 2 | 12,32 | 2 | 13,23 (mask 10) | BC21 | BC14 | BC18 |
| 2 | 12,31 | 2 | 12,23 (mask 9) | BC21 | BC14 | BC16 |
| 2 | 12,31,32 | 2 | 12,13,23 (mask 11) | BC21 | BC14 | BC16 |
| 2 | 12,13 | 1 | 21,23 (mask 12) | BC21 | BC14 | BC16 |
| 2 | 12,13,32 | 2 | 13,21,23 (mask 14) | BC21 | BC14 | BC16 |
| 2 | 12,13,31 | 2 | 12,21,23 (mask 13) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 2 | 12,13,31,32 | 2 | 12,13,21,23 (mask 15) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 2,3 | 21,31 | 2 | 21,31 | BC21 | BC15 | BC19, BC32 |
| 2,3 | 13 | 1 | 13 | BC21 | BC15, BC32 | BC19, BC32 |
| 2,3 | 13,21,31 | 3 | 13,21,31 | BC21 | BC15 | BC06, BC32 |
| 2,3 | 12 | 1 | 12 | BC21 | BC15, BC32 | BC19, BC32 |
| 2,3 | 12,21,31 | 3 | 12,21,31 | BC21 | BC15 | BC06, BC32 |
| 2,3 | 12,13 | 1 | 12,13 | BC21 | BC15, BC32 | BC19, BC32 |
| 2,3 | 12,13,21,31 | 3 | 12,13,21,31 | BC21 | BC15 | BC06, BC32 |
| 1 | 32 | 1 | 21 (mask 4) | BC21 | BC14 | BC17 |
| 1 | 31 | 1 | 23 (mask 8) | BC21 | BC14 | BC16 |
| 1 | 31,32 | 1 | 21,23 (mask 12) | BC21 | BC14 | BC16 |
| 1 | 23 | 1 | 12 (mask 1) | BC21 | BC14 | BC17 |
| 1 | 23,32 | 2 | 12,21 (mask 5) | BC21 | BC14 | BC16 |
| 1 | 23,31 | 2 | 12,23 (mask 9) | BC21 | BC14 | BC16 |
| 1 | 23,31,32 | 2 | 12,21,23 (mask 13) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 1 | 21 | 1 | 13 (mask 2) | BC21 | BC14 | BC16 |
| 1 | 21,32 | 2 | 13,21 (mask 6) | BC21 | BC14 | BC16 |
| 1 | 21,31 | 2 | 13,23 (mask 10) | BC21 | BC14 | BC18 |
| 1 | 21,31,32 | 2 | 13,21,23 (mask 14) | BC21 | BC14 | BC16 |
| 1 | 21,23 | 1 | 12,13 (mask 3) | BC21 | BC14 | BC16 |
| 1 | 21,23,32 | 2 | 12,13,21 (mask 7) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 1 | 21,23,31 | 2 | 12,13,23 (mask 11) | BC21 | BC14 | BC16 |
| 1 | 21,23,31,32 | 2 | 12,13,21,23 (mask 15) | BC21 | BC14 | BC18, C9_ENDPOINT |
| 1,3 | 23 | 1 | 12 | BC21 | BC15, BC32 | BC19, BC32 |
| 1,3 | 21 | 1 | 13 | BC21 | BC15, BC32 | BC19, BC32 |
| 1,3 | 21,23 | 1 | 12,13 | BC21 | BC15, BC32 | BC19, BC32 |
| 1,3 | 12,32 | 2 | 21,31 | BC21 | BC15 | BC19, BC32 |
| 1,3 | 12,23,32 | 3 | 12,21,31 | BC21 | BC15 | BC06, BC32 |
| 1,3 | 12,21,32 | 3 | 13,21,31 | BC21 | BC15 | BC06, BC32 |
| 1,3 | 12,21,23,32 | 3 | 12,13,21,31 | BC21 | BC15 | BC06, BC32 |
| 1,2 | 32 | 1 | 13 | BC21 | BC15, BC32 | BC19, BC32 |
| 1,2 | 31 | 1 | 12 | BC21 | BC15, BC32 | BC19, BC32 |
| 1,2 | 31,32 | 1 | 12,13 | BC21 | BC15, BC32 | BC19, BC32 |
| 1,2 | 13,23 | 2 | 21,31 | BC21 | BC15 | BC19, BC32 |
| 1,2 | 13,23,32 | 3 | 13,21,31 | BC21 | BC15 | BC06, BC32 |
| 1,2 | 13,23,31 | 3 | 12,21,31 | BC21 | BC15 | BC06, BC32 |
| 1,2 | 13,23,31,32 | 3 | 12,13,21,31 | BC21 | BC15 | BC06, BC32 |
| 1,2,3 | 21,31 | 2 | 21,31 | BC21, BC33 | BC21, BC33 | BC20, BC33 |
| 1,2,3 | 13,23 | 2 | 13,23 | BC21, BC33 | BC21, BC33 | BC20, BC33 |
| 1,2,3 | 13,21,23,31 | 3 | 13,21,23,31 | BC30 | BC30 | BC02, BC06, BC30, BC32 |
| 1,2,3 | 12,32 | 2 | 12,32 | BC21, BC33 | BC21, BC33 | BC20, BC33 |
| 1,2,3 | 12,21,31,32 | 3 | 12,21,31,32 | BC30 | BC30 | BC02, BC06, BC30, BC32 |
| 1,2,3 | 12,13,23,32 | 3 | 12,13,23,32 | BC30 | BC30 | BC02, BC06, BC30, BC32 |
| 1,2,3 | 12,13,21,23,31,32 | 3 | 12,13,21,23,31,32 | BC21 | BC14, BC31 | BC02, BC06, BC30, BC32 |

## Exceptional parameter and rank loci

### rank D = 0

The disjoint column supports are (31,32), (21,23), (12,13). Rank zero forces all six p entries to vanish.

Empty on the intermediate projective p chart.

Nodes: BC02. Manuscript labels: `appr3:sources`, `r3:boundarycoverage`.

### vanishing linear normal matrix on a pure axis

The residual quadrics have only the origin as a common zero; their multiplication matrix gives I_perp^4 subset I.

tau^4 in I; rank D is not the rank of this normal matrix.

Nodes: BC17. Manuscript labels: `r3short:pure-axis`, `r3short:initial`, `r3:normalexclusion`.

### D_9 = 0 on a two-coordinate axis

29 D_15 + 19 D_9 = 960 q; at D_9=0, rho^2+16q=415q/29.

Use the D_15 pivot; the ideal (p_bar_32 p_bar_13,p_bar_13^2) gives p_bar_13=0 and delta=0 on reduced arcs.

Nodes: BC19, BC32. Manuscript labels: `r3:boundarysources`.

### D_- = 0 or D_+ = 0

D_+ - D_- = 2k is a unit.

The walls are disjoint. The positive wall has residual determinant -3 ell/256; the negative wall uses (rho E,E^4).

Nodes: BC19, BC20, BC33. Manuscript labels: `r3:boundarysources`.

### canonical paired locus and k+q = 0

J=(p_bar_23 W_ell,p_bar_23 kappa,delta W_ell,delta W_rho,delta kappa,p_bar_23^2); det C_J=-5q rho/(1152 kappa_0^5).

kappa delta and tau belong to I. The three exact quotient-ring identities retain all cross terms.

Nodes: BC20, BC33. Manuscript labels: `r3:boundarysources`.

### simultaneously vanishing canonical elimination minors at material resonance

They require -kappa_0(1+kappa_0)=3/4, whereas this expression is at most 1/4.

Empty for real positive distinct principal values.

Nodes: BC20, BC33. Manuscript labels: `r3:boundaryreduction`, `r3:boundarysources`.

### critical missing row at its rank drop

The source and even equations require 1+mq=0; mq=(b/a)m^2>0.

Empty for positive distinct principal values.

Nodes: BC15, BC32. Manuscript labels: `r3:boundarysources`.

### critical full support, rank drop, point (1,1), exceptional direction or material wall

Use the full source matrix, its amplitude rank ideal, the rank-drop ideal and the exact Schur and resultant equations.

The direction and material cases are separated before inverting a minor.

Nodes: BC31. Manuscript labels: `appr3:sources`, `r3:boundarysources`.

### conic, w^2=k

The complete Q_6 derivative rows restore the two missing normal directions; determinant -6144 k^2 w x.

Full normal rank on this locus.

Nodes: BC18. Manuscript labels: `appr3whole:identities`.

### conic, w^2=-k

The conic equation forces x=z=0.

Redirect to single-coordinate support mask 10; use its soft quartic ideal.

Nodes: BC18. Manuscript labels: `appr3whole:identities`, `r3:boundarysources`.

### conic, w=0

The relative implicit graph gives delta in (w), kappa in (v_end,(p_bar_12+p_bar_21)^2); the intersection with (w) is w times this reference ideal.

kappa delta and tau in I; the transpose gives mask 13.

Nodes: BC18, C9_ENDPOINT. Manuscript labels: `appr3whole:endpoint`, `appr3whole:identities`.

### loss of a derivative coordinate or of an amplitude pair

The special-point nonzero coordinate set changes to another of the 73 incidence support patterns.

Use the indexed proper support; the vanishing coordinate remains a local variable on an approaching arc.

Nodes: BC14, BC15, BC16, BC19, BC21, BC30, BC31, BC32, BC33. Manuscript labels: `r3:boundarycoverage`, `r3:boundarysources`.

### ordinary intersections: u=0, v=-2, v=2 eta, Theta=0, extra tau wall and k+q=0

The local normal ideals have the radical of F_in, with the resonant quadratic normal change retained.

Containment followed by the generic algebraic differential exclusion.

Nodes: BC22, BC23, BC24, BC25, BC26, BC27, BC29. Manuscript labels: `r3:boundarysources`, `appr3:ordinaryalgebraicity`, `appr3ord:complete`, `appr3tau:wall`.

## Exact source interfaces

The source paths, SHA256 digests, arguments and prerequisite task identifiers are listed in `BOUNDARY_INDEX.json` and `PROOF_MAP.json`. Run a node and its prerequisites with:

```sh
python run.py --proof BC30 --work-dir work-bc30
```

Use a new directory for each independent execution. Nodes without a registered computational task state mathematical arguments in the manuscript.
