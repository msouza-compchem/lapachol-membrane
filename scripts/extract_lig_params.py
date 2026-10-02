#!/usr/bin/env python3
"""
Extract the CGenFF parameters specific to one ligand (LIG) from a merged
CHARMM-GUI forcefield.itp, so they can be added to a DIFFERENT topology
(the membrane's) without duplicating anything already there (water, ions,
lipid).

WHY THIS IS NECESSARY
----------------------
CHARMM-GUI's Solution Builder merges everything -- protein, water, ions, and
ligand -- into a single forcefield.itp when a ligand is present. That file is
NOT a clean "ligand-only" parameter set: it also repeats every water and ion
parameter GROMACS already knows from the membrane's own forcefield.itp.
Copying the whole file into the membrane project would raise "already
defined" errors from grompp, or worse, silently accept the first definition
it hits, so it does not fix the problem.

The two force-field files describe two different molecular universes and
share nothing except a *few* CGenFF atom types unique to the ligand
(prefixes like CG, OG, HG here) plus the bonded and pairtypes terms that
involve ONLY those types. Because bonded parameters (bonds, pairs, angles,
dihedrals) are, by construction, always specified between atoms belonging
to the same molecule, restricting to lines where every real (non-'X')
atom-type column belongs to the ligand's own type set is both necessary and
sufficient -- no cross terms with water or ions are ever needed, because no
bond/angle/dihedral in the ligand touches a water or ion atom.

WHAT THIS SCRIPT DOES
-----------------------
1. Reads the ligand's own atom types directly from LIG.itp's [ atoms ]
   section -- not from a hand-typed list, which is exactly the kind of list
   a human silently gets wrong (as happened with the first attempt here,
   which grepped for 7 types and missed 5).
2. Walks the merged forcefield.itp section by section, keeping only the
   lines in [atomtypes], [bondtypes], [pairtypes], [angletypes], and
   [dihedraltypes] whose real (non-'X') atom-type columns are a subset of
   the ligand's own types.
3. Writes a single, self-contained .itp with those sections, ready to be
   #include'd in the target topol.top BEFORE the ligand's own moleculetype
   file (LIG.itp) -- exactly the position CHARMM-GUI itself uses.

Usage:
    python extract_lig_params.py \\
        --lig  path/to/LIG.itp \\
        --ff   path/to/merged/forcefield.itp \\
        --out  lig_cgenff_params.itp
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

SECTION_RE = re.compile(r"^\s*\[\s*(\w+)\s*\]")


def read_ligand_types(lig_itp: Path) -> set[str]:
    """
    Pull the set of unique atom types straight from LIG.itp's [ atoms ]
    table. Column layout: nr type resnr residue atom cgnr charge mass ...
    so the type is the second whitespace-separated token.
    """
    types: set[str] = set()
    in_atoms = False
    for raw in lig_itp.read_text().splitlines():
        line = raw.split(";", 1)[0].strip()  # drop end-of-line comments
        m = SECTION_RE.match(raw)
        if m:
            in_atoms = m.group(1).lower() == "atoms"
            continue
        if in_atoms and line:
            parts = line.split()
            if len(parts) >= 2:
                types.add(parts[1])
    return types


def relevant(tokens: list[str], n_type_cols: int, ligtypes: set[str]) -> bool:
    """
    A bonded-parameter line is ligand-specific if every REAL (non-'X')
    atom-type column is one of the ligand's own types, and at least one
    column is a real type (never all wildcards).
    """
    cols = tokens[:n_type_cols]
    real = [c for c in cols if c != "X"]
    if not real:
        return False
    return all(c in ligtypes for c in real)


# section name -> number of leading columns that are atom-type names
TYPE_COLS = {
    "atomtypes": 1,
    "bondtypes": 2,
    "pairtypes": 2,
    "angletypes": 3,
    "dihedraltypes": 4,
}


def extract(ff_itp: Path, ligtypes: set[str]) -> dict[str, list[str]]:
    """Walk forcefield.itp once, bucketing the lines that survive the filter."""
    buckets: dict[str, list[str]] = {k: [] for k in TYPE_COLS}
    section = None

    for raw in ff_itp.read_text().splitlines():
        m = SECTION_RE.match(raw)
        if m:
            name = m.group(1).lower()
            section = name if name in TYPE_COLS else None
            continue
        if section is None:
            continue

        stripped = raw.split(";", 1)[0].strip()
        if not stripped:
            continue

        tokens = stripped.split()
        if len(tokens) < TYPE_COLS[section]:
            continue

        if relevant(tokens, TYPE_COLS[section], ligtypes):
            buckets[section].append(raw.rstrip())

    return buckets


HEADER = """\
;; Ligand-specific CGenFF parameters, extracted automatically by
;; extract_lig_params.py from a merged CHARMM-GUI forcefield.itp.
;;
;; Contains ONLY lines whose atom-type columns belong to the ligand's own
;; CGenFF type set -- no water, ion, or lipid parameters are duplicated
;; here, since every line was required to be self-contained to the ligand.
;;
;; Include this file BEFORE the ligand's own moleculetype (e.g. LIG.itp).
"""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lig", required=True, type=Path, help="path to LIG.itp")
    ap.add_argument("--ff", required=True, type=Path,
                    help="path to the merged forcefield.itp that has the ligand's parameters")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    ligtypes = read_ligand_types(args.lig)
    print(f"Ligand atom types found in {args.lig.name} ({len(ligtypes)}):")
    print("  " + ", ".join(sorted(ligtypes)))

    buckets = extract(args.ff, ligtypes)

    print("\nLines extracted per section:")
    for name, lines in buckets.items():
        print(f"  {name:14s}: {len(lines)}")

    if not buckets["atomtypes"]:
        print("\n! No atomtypes matched. Check that --lig and --ff correspond "
              "to the same ligand, and that LIG.itp really has an [ atoms ] "
              "section with the type in column 2.")

    with args.out.open("w") as fh:
        fh.write(HEADER + "\n")
        # Emit in the order GROMACS expects a topology's non-bonded/bonded
        # parameter file to be structured.
        order = ["atomtypes", "bondtypes", "pairtypes", "angletypes", "dihedraltypes"]
        for name in order:
            lines = buckets[name]
            if not lines:
                continue
            fh.write(f"[ {name} ]\n")
            fh.write("\n".join(lines) + "\n\n")

    print(f"\nWrote {args.out}")
    print("\nNext steps:")
    print(f"  1. Copy {args.out.name} into the membrane project's toppar/ folder.")
    print(f"  2. In the membrane's topol.top, add this line BEFORE the LIG.itp include:")
    print(f'       #include "toppar/{args.out.name}"')
    print(f"  3. Copy LIG.itp itself into the same toppar/ folder.")
    print(f"  4. Add another include line for LIG.itp right after it.")
    print(f"  5. Add 'LIG   1' (or however many copies) to the [ molecules ] "
          f"section at the end of topol.top -- but only AFTER you've actually "
          f"inserted the ligand's coordinates into the .gro with gmx insert-molecules.")


if __name__ == "__main__":
    main()
