# SageMath calculations for Rank3

This directory contains the calculations for rank three in the accompanying paper.
The fixed environment is SageMath 10.9, Python 3.12.14 and SymPy 1.14.0
on Linux x86_64. Windows users can run the calculations in WSL.

Create the shared environment using the instructions in `../README.md`.
From this directory, run:

```sh
../sage-environment/bin/python run.py --work-dir ../rebuild-rank3 --workers 2
```

The work directory must not already exist for a new run. Results and logs
are written below its `run/` directory. To continue an interrupted run,
repeat the same command with `--resume`. Use `--workers 1` to reduce
concurrent memory use.

`PROOF_MAP.md` lists the mathematical parts and manuscript locations.
`PROOF_MAP.json` records their exact algebraic statements, coefficient
rings, localizations, source files and complete dependencies. `SYMBOLS.md`
relates manuscript notation to the local Sage names. List the available
parts with `../sage-environment/bin/python run.py --list-proofs`; select
one part and its prerequisites with `--proof ID` and a new work directory.

`sources/` contains the mathematical programs and their imported modules.
`inputs/` contains the required exact coefficient data and witnesses.
`runtime/` contains the execution support. `tasks.json` fixes the dependency
graph. `MANIFEST.json` lists the required files, hashes and consuming tasks.
The `.tex` files in `inputs/paper/` are formula excerpts; the complete
manuscript and supplement are in `../paper/`.

The calculations use exact integer, rational or algebraic arithmetic.
Provide at least 16 GiB of RAM for the complete calculation.

## The relative conic endpoint in Appendix C.9

This calculation is included in the complete run. To reproduce it separately
with all its prerequisites, run from the repository root:

```sh
./sage-environment/bin/python Rank3/run.py --proof C9_ENDPOINT --work-dir ./rebuild-endpoint --workers 2
```

The selection has 42 tasks: 23 mathematical executions and 19 copies of
freshly generated inputs. Nine mathematical checks and nine input copies
were added for this proof. The preceding source calculations are shared
with the rest of Rank3.

`results/conic_endpoint/` contains the recorded exact outputs: the
1444-term relative derivative polynomial through its exact ideal decomposition, the
source denominator census, the graph identities, the endpoint normal
matrix and Jacobian coefficients, and the open conic normal matrix.
Their producing tasks and hashes appear in `PROOF_MAP.json`. A new run
creates its own outputs under `reproduction/conic_endpoint/` in the work
directory. The local ideal argument is in Appendix C.9 of the manuscript.
The polynomial coefficients are stored once in
`generic_G_ideal_decomposition_sparse.json`, using
`G = t*Gt + kappa*Gkappa`. The program also writes the direct polynomial
coefficient list when rerun.

The columns of `endpoint_normal_matrix.txt` are
`(x,z,w,e0,e1,e2,delta,t,kappa)`. Its first 28 rows are
`R0,...,R5,Q0,H0,...,H17,Q1,Q5,Q6`. The last nine rows are grouped
by gauge indices `1,2,4` and, within each group, `Q1,Q5,Q6`.
The six pivot rows, numbered from zero, are `2,3,6,25,27,33`;
row 33 is the complete `G(p21,Q6)`. The endpoint root equation,
source hashes, pivot determinants, tangent matrix and residual
coefficients are in `endpoint_unit_and_obstructions.json`.
The open conic matrix in `conic_normal_matrix.json` carries its
row and column labels.
