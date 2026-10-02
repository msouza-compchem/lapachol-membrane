#!/usr/bin/env python3
"""
Phase A -- LogP and passive-permeability descriptors for lapachol.

WHY THIS IS THE FIRST STEP
----------------------------
Before building a membrane and running expensive umbrella sampling, ask the
cheap question first: does this molecule even LOOK like something that crosses
a lipid bilayer passively? LogP, polar surface area and hydrogen-bond count are
the classic screen -- the same properties behind Lipinski's Rule of Five for
oral drug absorption.

None of this replaces the simulation. It sets the expectation the simulation
should either confirm or contradict -- and if it contradicts, THAT is
interesting and worth investigating, not a reason to distrust the calculation.

Usage:
    python 01_logp_screen.py
    python 01_logp_screen.py --smiles "CC(C)=CCC1=C(O)C(=O)c2ccccc2C1=O" --name lapachol
"""

from __future__ import annotations

import argparse

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski

RDLogger.DisableLog("rdApp.*")

LAPACHOL_SMILES = "CC(C)=CCC1=C(O)C(=O)c2ccccc2C1=O"


def describe(smiles: str, name: str) -> None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise SystemExit(f"Could not parse SMILES: {smiles}")

    # Crippen LogP: an atom-contribution method, the same family of estimate
    # used throughout medicinal chemistry when an experimental value is not
    # in hand. Fast, and good enough to place a molecule on the right side of
    # the permeability question.
    logp = Crippen.MolLogP(mol)

    # Topological polar surface area: the sum of surface contributed by polar
    # atoms (mainly O, N and their attached H). Above ~140 A^2, passive
    # membrane permeability drops sharply -- this threshold is empirical,
    # from analysis of orally absorbed drugs.
    tpsa = Descriptors.TPSA(mol)

    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    mw = Descriptors.MolWt(mol)
    rotatable = Descriptors.NumRotatableBonds(mol)

    print(f"\n{'=' * 60}")
    print(f"{name}")
    print(f"{'=' * 60}")
    print(f"  SMILES                 : {smiles}")
    print(f"  Molecular weight       : {mw:.1f} Da")
    print(f"  LogP (Crippen)         : {logp:.2f}")
    print(f"  TPSA                   : {tpsa:.1f} A^2")
    print(f"  H-bond donors          : {hbd}")
    print(f"  H-bond acceptors       : {hba}")
    print(f"  Rotatable bonds        : {rotatable}")

    print(f"\n  Interpretation:")
    if logp < 0:
        print(f"    LogP < 0: unusually hydrophilic for passive diffusion.")
    elif logp <= 5:
        print(f"    LogP in the 0-5 range: consistent with passive membrane "
              f"permeability.")
    else:
        print(f"    LogP > 5: unusually lipophilic; may partition INTO the "
              f"membrane rather than crossing it (a permeability liability, "
              f"not an advantage).")

    if tpsa > 140:
        print(f"    TPSA > 140 A^2: this alone would predict POOR passive "
              f"permeability, regardless of LogP.")
    else:
        print(f"    TPSA <= 140 A^2: does not, by itself, predict a "
              f"permeability problem.")

    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    print(f"\n  Lipinski Rule of Five violations: {violations} of 4")
    print(f"  (For context, not a verdict -- Ro5 was built for oral-drug-like")
    print(f"  molecules in general, not quinones specifically, and predicts")
    print(f"  ABSORPTION LIKELIHOOD, not passive bilayer crossing directly.)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--smiles", default=LAPACHOL_SMILES)
    ap.add_argument("--name", default="lapachol")
    args = ap.parse_args()
    describe(args.smiles, args.name)

    print(f"\n{'-' * 60}")
    print("For comparison, three reference points from pharmacology:")
    print(f"{'-' * 60}")
    for name, smi in [
        ("aspirin (crosses membranes well)", "CC(=O)Oc1ccccc1C(=O)O"),
        ("glucose (does NOT cross passively -- needs a transporter)",
         "OCC1OC(O)C(O)C(O)C1O"),
        ("cholesterol (very lipophilic, membrane-embedded)",
         "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"),
    ]:
        describe(smi, name)


if __name__ == "__main__":
    main()
