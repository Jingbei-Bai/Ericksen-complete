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
