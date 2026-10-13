"""Compile the proof diagram with TeX and render its README image.

Requirements: Python 3, pdflatex (with standalone and TikZ), pdftoppm.
Run from any directory: python assets/build_proof_figure.py
Executable paths can be supplied through PDFLATEX and PDFTOPPM.
"""
from pathlib import Path
import os
import shutil
import subprocess

ASSETS = Path(__file__).resolve().parent
BUILD = ASSETS.parent / ".figure-build"


def executable(variable, default):
    value = os.environ.get(variable) or shutil.which(default)
    if not value:
        raise SystemExit(f"Set {variable} to the path of {default}, or add it to PATH.")
    return value


def main():
    BUILD.mkdir(exist_ok=True)
    subprocess.run(
        [executable("PDFLATEX", "pdflatex"), "-interaction=nonstopmode",
         "-halt-on-error", f"-output-directory={BUILD}", "proof_structure.tex"],
        cwd=ASSETS, check=True, stdout=subprocess.DEVNULL,
    )
    shutil.copy2(BUILD / "proof_structure.pdf", ASSETS / "proof_structure.pdf")
    subprocess.run(
        [executable("PDFTOPPM", "pdftoppm"), "-singlefile", "-png", "-r", "170",
         str(ASSETS / "proof_structure.pdf"), str(ASSETS / "proof_structure")],
        check=True,
    )
    print("Created assets/proof_structure.pdf and assets/proof_structure.png")


if __name__ == "__main__":
    main()
