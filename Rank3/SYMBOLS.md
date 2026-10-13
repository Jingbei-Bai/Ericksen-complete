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
| Boundary coordinates | `p12`, ..., `p32` on a homogeneous chart | \(\overline p_{ij}=\tau\widehat p_{ij}\) |
| Weighted homogeneous derivative vector | `v1`, `v2`, `v3` | \(\boldsymbol v=\tau^2\mathsf R\) |
| Ordinary chart | ordinary `R` vector | \(\boldsymbol R=\tau\mathsf R\) |
| Intermediate chart with selected entry `j` | normalized `R` vector | \(\overline{\boldsymbol R}=\boldsymbol v/(\boldsymbol v)_j\) |
| Ratio incidence coefficients | `alpha`, ratio `beta` | \(\alpha=(\zeta-\eta)/(\zeta(\eta+1)),\quad\beta=\eta(\zeta-1)/(\zeta(\eta+1))\) |
| T chart in `independent_critical_inactive_fullRR_beta_polynomial_reader.py` | the fourth field generator `bd_beta` | \(\beta_{\mathrm T}=\tau/\sigma^2\), with \(\kappa=\sigma\beta_{\mathrm T}\) |
| Paired and axis incidence charts | `l` or `ell` | \(\ell=\overline p_{32}\) |
| Canonical paired chart | amplitude `s` | \(\overline p_{23}\) |
| Two-coordinate axis chart | amplitude `s` | \(\overline p_{13}\) |

The parameters of boundary charts, including a blowup coordinate named `sigma`, are local chart coordinates. They are distinct from the curve slope `S`. The source polynomial parents specify their generators explicitly.

The valuation orders are \(n=\operatorname{ord}\tau\) and
\(s_0=\min_j\operatorname{ord}(\boldsymbol v)_j\). The material
constant \(\kappa_0\) and the blowup coordinate \(\kappa\) have the
definitions in Appendix C.1 and C.2. `rank_D` in `BOUNDARY_INDEX.json`
is the rank of \(D(\overline p)\); a source or normal matrix has its
own row and column labels in the source data.

The relative conic calculations use the following names in `reproduction/conic_endpoint/`:

| Sage name | Manuscript meaning |
|---|---|
| `k`, `r`, `m`, `q` | \(a\kappa_0,\ \kappa_0,\ 1+a\kappa_0,\ 1+\kappa_0\) |
| `eta` or `H`, `zeta` or `Z` | \(\eta,\zeta\) |
| `w`, `t`, `kappa`, `delta` | \(\overline p_{23},\overline{\boldsymbol R}_2,\kappa,\delta\) |
| `e0`, `e1`, `e2` | \(e_0^{\mathrm{alg}},e_1^{\mathrm{alg}},e_2^{\mathrm{alg}}\) |
| free-coordinate `v` and `s` in endpoint output | \(v_{\mathrm{end}}=e_2^{\mathrm{alg}}-f_2,\ s=x+z\) |
| residuals `B`, `F` | \(\mathcal B_{\mathrm{end}},\mathcal F_{\mathrm{end}}\) |
