# Lapachol → NQO1: a reproducible DFT / docking / MD pipeline

An end-to-end, laptop-scale computational study of **lapachol**, a naphthoquinone from
the Brazilian *ipê* tree (*Handroanthus* spp.), and its interaction with **human
NAD(P)H:quinone oxidoreductase 1 (NQO1)**.

---

## 1. Scientific question

NQO1 catalyses the two-electron reduction of quinones. For most quinones this is a
detoxification route, but for some — notably β-lapachone, the cyclised isomer of
lapachol — the resulting hydroquinone is unstable and re-oxidises spontaneously,
driving a futile redox cycle that depletes NAD(P)H and generates reactive oxygen
species. Because NQO1 is over-expressed in several solid tumours, this turns the
enzyme into a bioactivation switch rather than a protective one.

**Question addressed here:** can inexpensive, reproducible calculations rationalise
the redox behaviour of lapachol, and does the molecule occupy the NQO1 catalytic site
in a geometry compatible with hydride transfer from the FAD cofactor?

**Approach:** three stages, each independently reproducible.

| Stage | Method | Software | Output |
|-------|--------|----------|--------|
| 1 | DFT geometry optimisation, frequencies, frontier orbitals, adiabatic electron affinity | ORCA 6 | Electronic descriptors |
| 2 | Molecular docking into the NQO1 catalytic site | AutoDock Vina 1.2 | Binding pose + score |
| 3 | Molecular dynamics of the complex, 100 ns | GROMACS 2024 | Pose stability, contacts |

---

## 2. System

- **Ligand:** lapachol, C15H14O4, 19 heavy atoms, 2 rotatable bonds.
- **Receptor:** human NQO1, PDB **1D4A**, 1.7 Å resolution, homodimer with one FAD
  per monomer. The catalytic site is formed at the dimer interface; Tyr128 and
  Phe232 gate the site and are known to be conformationally mobile.
- **Reference ligand for validation:** duroquinone, resolved in PDB **1DXO**, used
  here as a redocking control.

---

## 3. Reproducing this work

```bash
git clone https://github.com/<user>/lapachol-nqo1.git
cd lapachol-nqo1
conda env create -f environment.yml
conda activate lapachol
```

Then follow the stages in order:

```bash
bash 01_dft/run_orca.sh            # ~3 h on 4 cores
bash 02_docking/run_docking.sh     # ~5 min
bash 03_md/build_system.sh         # ~20 min setup; production runs on GPU
python 04_analysis/orca_descriptors.py
python 04_analysis/md_analysis.py
```

Each stage writes its figures to `figures/`.

---

## 4. Computational details

**DFT.** Geometries optimised at B3LYP-D4/def2-SVP with the RIJCOSX approximation;
harmonic frequencies computed at the same level to confirm every structure is a true
minimum (no imaginary frequencies). Single-point energies at B3LYP-D4/def2-TZVP with
the CPCM implicit water model. The adiabatic electron affinity is obtained from
separately optimised neutral and radical-anion structures.

**Docking.** Receptor prepared from 1D4A with FAD retained as part of the rigid
receptor, since the isoalloxazine ring forms one wall of the substrate site. Search
box centred on the duroquinone position from the aligned 1DXO structure. Protocol
validated by redocking duroquinone and checking the RMSD against the crystal pose
(target: below 2.0 Å).

**MD.** CHARMM36m protein force field, CGenFF ligand parameters, TIP3P water,
0.15 M NaCl, cubic box with 1.0 nm padding. Steepest-descent minimisation, 100 ps NVT
at 300 K (V-rescale), 100 ps NPT at 1 bar (C-rescale), then 100 ns production with a
2 fs time step and LINCS constraints on bonds to hydrogen.

---

## 5. Hardware and why it is split this way

Stages 1 and 2 were run on a consumer laptop (4-core i5, 8 GB RAM). Stage 3 was run
on a free-tier cloud GPU: the solvated complex contains roughly 50,000 atoms, which
on CPU alone delivers only 1–3 ns/day, but 30–80 ns/day on a single T4. Trajectory
files are not versioned here; only inputs, parameters and analysis code are, which is
sufficient to regenerate every result.

---

## 6. Results

*(Fill in as the calculations finish. Keep this section to three or four short
paragraphs plus the figures — a reader should be able to see the outcome without
opening any other file.)*

| Descriptor | Value | Unit |
|------------|-------|------|
| E(HOMO) | | eV |
| E(LUMO) | | eV |
| HOMO–LUMO gap | | eV |
| Adiabatic electron affinity | | eV |
| Electrophilicity index ω | | eV |
| Vina score (best pose) | | kcal/mol |
| Backbone RMSD, last 50 ns | | nm |
| Ligand RMSD, last 50 ns | | nm |

---

## 7. Limitations

Implicit solvation and a single conformer per species; Vina scores are not free
energies and are used only to rank poses; the MD samples a single 100 ns replica,
which is enough to test pose stability but not to converge binding thermodynamics.
Hydride transfer itself is not modelled — that would require QM/MM.

---

## 8. References

See [`docs/REFERENCES.md`](docs/REFERENCES.md).

## 9. Licence

MIT — see [`LICENSE`](LICENSE).
