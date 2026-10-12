# Ericksen's classification problem: manuscript and SageMath calculations

Release date: 12 October 2026.

Clone the repository or extract the release archive. Its contents are ready to use:

```text
README.md
paper/
    main.pdf
    main.tex
    supplement.tex
    figures/
Rank012/
Rank3/
VERIFICATION.json
SHA256SUMS.txt
```

`paper/main.pdf` is the complete preprint (210 pages). Supplementary
calculations S1-S9 follow Appendix E, beginning on page 168. The document
has one abstract and one table of contents, and all cross-references stay
within the same PDF. Its complete LaTeX sources and three vector figure
backgrounds are in the same directory. `supplement.tex` is included by
`main.tex`; it is not a separate document. The revised Section 8 is
integrated into `main.tex`.

`Rank012/` contains the calculations for ranks zero, one and two.
`Rank3/` contains the calculations for rank three. Each directory includes
all required programs, inputs, execution support, the task graph,
the environment lock and the manuscript index.

## Reproduce the calculations

Run the commands below from the extracted release directory on Linux
x86_64, or in WSL on Windows. Use Conda or Miniforge to install the fixed
environment. Installation requires access to the package repositories.

```sh
conda create --prefix ./sage-environment --file Rank012/environment-linux-64.lock -y
./sage-environment/bin/python Rank012/run.py --work-dir ./rebuild-rank012 --workers 2
./sage-environment/bin/python Rank3/run.py --work-dir ./rebuild-rank3 --workers 2
```

Both directories contain the same environment lock: SageMath 10.9,
Python 3.12.14 and SymPy 1.14.0, with all dependency versions and builds
fixed. Provide at least 16 GiB of RAM; use `--workers 1` to reduce concurrent
memory use.

Use a new work directory for each new run. Outputs and logs are created
under its `run/` subdirectory. Continue an interrupted run by adding
`--resume` to the original command. The runner resolves dependencies and
creates the required generated inputs.

To inspect or reproduce one mathematical part:

```sh
./sage-environment/bin/python Rank3/run.py --list-proofs
./sage-environment/bin/python Rank3/run.py --proof ID --work-dir ./rebuild-part --workers 1
```

Replace `ID` with an entry returned by `--list-proofs`.

## Correspondence with the paper

Each calculation directory contains:

- `PROOF_MAP.md`: mathematical parts and their manuscript locations.
- `PROOF_MAP.json`: exact algebraic statements, coefficient rings,
  localizations, manuscript excerpts and dependencies.
- `SYMBOLS.md`: manuscript notation and local Sage names.
- `tasks.json`: the registered tasks and dependency graph.
- `MANIFEST.json`: file hashes and the tasks requiring each source or input.

The `.tex` files under `inputs/paper/` are formula excerpts used by the
programs. The complete sources are `paper/main.tex` and
`paper/supplement.tex`. The supplement records the important matrices,
determinants, local ideals, resultants and intersection orders, including
all 425 auxiliary expressions and 20 entries of the moment matrix.
One required 468-byte coefficient input in `Rank3/inputs/` has a `.log`
extension.

Check the indices against the full sources with:

```sh
./sage-environment/bin/python Rank012/paper_index.py --paper paper/main.tex --supplement paper/supplement.tex
./sage-environment/bin/python Rank3/paper_index.py --paper paper/main.tex --supplement paper/supplement.tex
```

## Completed runs and file integrity

| Calculation | Mathematical tasks covered by the recorded runs | Generated-data copies |
|---|---:|---:|
| Ranks 0, 1 and 2 | 63 / 63 passed | 0 |
| Existing Rank 3 calculations | 337 / 337 passed | 129 / 129 completed |
| Added conic endpoint calculations | 9 / 9 passed | 9 / 9 completed |

The existing packages were executed in clean work directories. Their
mathematical sources, inputs, arguments and execution programs are unchanged.
Rank 3 includes full exact integer interpolation at 2039 and 965 nodes,
with both complete independent readers passing.

The added Appendix C.9 calculations were run from another empty work
directory together with their entire dependency closure: 23 mathematical
executions and 19 generated-data copies, all completed. The current Rank3
task graph contains 484 tasks. The verification combines the existing
complete run and this clean endpoint dependency run.

To reproduce the new part alone:

```sh
./sage-environment/bin/python Rank3/run.py --proof C9_ENDPOINT --work-dir ./rebuild-endpoint --workers 2
```

The exact outputs and both coefficient matrices are in
`Rank3/results/conic_endpoint/`. The manuscript proves the local ideal
implication from these identities in Appendix C.9.

`VERIFICATION.json` contains the completed-run evidence, environment
versions, manuscript hashes and directory hashes. The directory hash is
SHA-256 of the sorted sequence of relative POSIX paths, a NUL separator,
each file's SHA-256 and a newline. This revision integrates the expanded Section 8 and adds the
recomputed relative conic endpoint calculations.

Verify the extracted files from the release directory:

```sh
sha256sum --check SHA256SUMS.txt
```

## Compile the complete preprint

Use pdfLaTeX with TeX Live 2025. From `paper/`, run the following command
three times:

```sh
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

This compiles the paper and supplementary calculations together. All
formulas retain their original section and equation numbers.
