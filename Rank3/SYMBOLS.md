# Mathematical symbols

The Sage names are local to the coefficient rings indicated below.

| Computational interface | Sage names | Manuscript meaning |
|---|---|---|
| Rank-three spectral coefficient rings | `r`, `q`, `m` | \(\kappa_0,\ 1+\kappa_0,\ 1+a\kappa_0\) |
| Joint fields before the ratio change | `h`, `z` | \(h_2/h_1,\ h_3/h_1\) |
| Ratio-curve rings with generators `(a,r,H,Z,S,N)` | `H`, `Z` | \(\eta=r h,\quad\zeta=r z/q\) |
| The same ratio-curve rings | `S`, `N` | \(\varsigma=d\zeta/d\eta,\quad d\varsigma/d\eta\) |
| Moment elimination | `F`, `G`, `B` | \(\mathcal P_1,\mathcal P_2,\mathcal P_0\) |
| Boundary coefficient rings | `rho` | Local boundary coordinate \(\varrho=\overline R_2\) on \(\overline R_3=1\) |
| State matrices | Python indices `0,1,2` | Principal labels \(1,2,3\) |

The parameters of boundary charts, including a blowup coordinate named `sigma`, are local chart coordinates. They are distinct from the curve slope `S`. The source polynomial parents specify their generators explicitly.
