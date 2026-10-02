# lapachol-membrane

**Free-energy profile of lapachol crossing a POPC bilayer (GROMACS, CHARMM36, umbrella sampling).**
Project 3 of my computational chemistry portfolio ([msouza-compchem](https://github.com/msouza-compchem)). Work in progress, documented as an open notebook.

*Resumo em português: calculo o perfil de energia livre ΔG(z) da permeação do lapachol (naftoquinona do ipê) em uma bicamada de POPC, com dinâmica molecular e umbrella sampling. O projeto está em andamento e este repositório registra cada etapa, com parâmetros, scripts e limitações. O caderno visual está em `docs/` (GitHub Pages).*

## Scientific question

Lapachol is a naphthoquinone from *Handroanthus* (ipê) studied for its biological activity. In [`lapachol-nqo1`](https://github.com/msouza-compchem/lapachol-nqo1) I docked it to intracellular human NQO1; to get there, the molecule has to cross a membrane.

> What is the free-energy profile, ΔG(z), of lapachol crossing a lipid bilayer, and at which depth does it prefer to sit?

Expected outcomes (to be tested, not results): a free-energy minimum near the headgroup/tail interface and a moderate barrier at the bilayer center. A comparison with the computed LogP is planned.

**What this project will not answer:** ΔG(z) alone does not give permeability (that also needs the diffusion profile D(z), via the solubility-diffusion model); only one protonation state of the weak acid is treated at first; a single lipid composition with a classical force field is not a real cell membrane.

## Status

- [x] Bilayer built with CHARMM-GUI (120 POPC, 13 K⁺, 13 Cl⁻, TIP3P water)
- [x] 100 ps test run with production parameters (stable; see below)
- [x] Speed benchmark on a Colab L4 GPU
- [ ] 100 ns production run of the pure bilayer (area per lipid, thickness and density profiles to check equilibration)
- [ ] Steered MD along z with lapachol
- [ ] Umbrella sampling: 5 ns per window to screen, then 10–20 ns for production
- [ ] WHAM with bootstrap, histogram overlap and block convergence

## System and parameters

| Item | Value |
|---|---|
| Software | GROMACS 2023.3 (local, CPU) and 2024.5 conda-forge build with CUDA (Colab) |
| Builder | CHARMM-GUI |
| Force field | CHARMM36 (force-switch 1.0–1.2 nm, PME) |
| Bilayer | 120 POPC (60 per leaflet), 13 K⁺, 13 Cl⁻, 32,267 atoms |
| Box (after equilibration) | about 6.31 × 6.31 × 7.75 nm |
| Water | TIP3P (CHARMM) |
| Temperature / pressure | 303.15 K (v-rescale) / 1 bar, semi-isotropic (C-rescale) |
| Time step | 2 fs, h-bonds constrained (LINCS) |

## Results so far

100 ps test run, CHARMM-GUI production parameters, pure bilayer:

| Quantity | Result |
|---|---|
| Temperature | 303.05 K (target 303.15) |
| Pressure | 8.6 ± 9.4 bar (large fluctuations are expected for this system size) |
| Box X = Y | 6.310 nm, negligible drift |
| Box Z | 7.754 nm, drift 0.015 nm |
| Density | 1017.7 kg/m³ |
| Area per lipid | about 66 Å² (not yet evidence of equilibrium; see the 100 ns run) |
| LINCS warnings | none |

Speed: 4.8 ns/day on 4 CPU cores (WSL) and **270 ns/day on a Colab L4 GPU** (100 ps benchmark).

## Repository layout

```
inputs/bilayer/   starting structure, topology, index and production .mdp
results/          small outputs (energies, logs) from tests
scripts/          helper scripts (frame selection, analysis)
notebooks/        Colab notebooks
docs/             open notebook (GitHub Pages)
```

Trajectories (`.xtc`, `.trr`), checkpoints and binary run input files are not tracked because of size. They are kept on Google Drive and can be shared on request.

## How to reproduce the bilayer run

```bash
cd inputs/bilayer
gmx grompp -f prod_100ns.mdp -c step6.6_equilibration.gro -p topol_bilayer.top -n index.ndx -o prod_100ns.tpr
gmx mdrun -deffnm prod_100ns -cpi prod_100ns.cpt -cpt 10 -maxh 11 -nb gpu -pme gpu -bonded gpu -update gpu
```

The `.mdp` file is the CHARMM-GUI `step7_production.mdp` with `nsteps = 50000000` (100 ns).

## Planned umbrella sampling

Pulling along z with `pull-coord1-geometry = direction`, 43 windows spaced 0.15 nm (z from −3.15 to +3.15 nm), force constant 1000 kJ/mol/nm² (to be adjusted if histograms do not overlap). Convergence will be assessed by WHAM over time blocks and bootstrap errors, following common practice for permeation PMFs: Cordeiro (2018) is a close methodological reference.

## References

1. Torrie GM, Valleau JP. *J Comput Phys* 1977, 23, 187.
2. Kumar S, et al. *J Comput Chem* 1992, 13, 1011.
3. Hub JS, de Groot BL, van der Spoel D. *J Chem Theory Comput* 2010, 6, 3713.
4. Marrink SJ, Berendsen HJC. *J Phys Chem* 1994, 98, 4155.
5. Jo S, et al. *J Comput Chem* 2008, 29, 1859 (CHARMM-GUI).
6. Cordeiro RM. *J Phys Chem B* 2018, 122, 8211.

## Author

Marcelo Souza. Chemistry and biology teacher building a computational chemistry portfolio. Code under MIT license (add a `LICENSE` file).
