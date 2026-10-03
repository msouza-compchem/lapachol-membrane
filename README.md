# lapachol-membrane

**Free-energy profile of lapachol crossing a POPC bilayer (GROMACS, CHARMM36, umbrella sampling).**
Project 3 of my computational chemistry portfolio ([msouza-compchem](https://github.com/msouza-compchem)). Work in progress, documented as an open notebook: https://msouza-compchem.github.io/lapachol-membrane/

*Resumo em português: calculo o perfil de energia livre ΔG(z) da permeação do lapachol (naftoquinona do ipê) em uma bicamada de POPC, com dinâmica molecular e umbrella sampling. O projeto está em andamento. Este repositório registra cada etapa, com parâmetros, scripts, resultados e limitações. O caderno visual está em `docs/`.*

## Scientific question

Lapachol is a naphthoquinone from *Handroanthus* (ipê) studied for its biological activity. In [`lapachol-nqo1`](https://github.com/msouza-compchem/lapachol-nqo1) I docked it to the intracellular enzyme human NQO1; to get there, the molecule has to cross a membrane.

> What is the free-energy profile, ΔG(z), of lapachol crossing a lipid bilayer, and at which depth does it prefer to sit?

Expected outcomes (to be tested, not results): a free-energy minimum near the headgroup/tail interface and a moderate barrier at the bilayer center. A comparison with the computed LogP is planned.

**What this project will not answer:** ΔG(z) alone does not give permeability (that also needs the diffusion profile D(z), via the solubility-diffusion model); only the neutral form of the weak acid is treated at first; one lipid composition with a classical force field is not a real cell membrane.

## Status

- [x] Bilayer built with CHARMM-GUI (120 POPC, 13 K⁺, 13 Cl⁻, TIP3P water)
- [x] 100 ns production of the pure bilayer (area per lipid 64.4 ± 0.3 Å²)
- [x] Lapachol placed in water at 3.6 nm, 1 ns equilibration, 1.8 ns steered MD to the bilayer center
- [x] 37 window frames extracted (every 0.1 nm, 3.6 nm to 0)
- [ ] Stage 1: 37 windows × 5 ns (screening) — started 2026-10-02
- [ ] WHAM with bootstrap, histogram overlap and block convergence
- [ ] Stage 2: extend windows to 10–20 ns where needed

## System and parameters

| Item | Value |
|---|---|
| Software | GROMACS 2023.3 (local, CPU) and 2024.5 conda-forge build with CUDA (Colab L4) |
| Builder | CHARMM-GUI |
| Force field | CHARMM36 (force-switch 1.0–1.2 nm, PME); lapachol with CGenFF |
| Bilayer | 120 POPC (60 per leaflet), 13 K⁺, 13 Cl⁻, 32,267 atoms |
| System with ligand | 32,179 atoms (5,347 TIP3P waters after removing 40 clashing with lapachol) |
| Box after 100 ns | about 6.2 × 6.2 × 8.0 nm |
| Temperature / pressure | 303.15 K (v-rescale) / 1 bar, semi-isotropic (C-rescale) |
| Time step | 2 fs, h-bonds constrained (LINCS) |
| Reaction coordinate | z distance between lapachol (center of mass) and the POPC hydrophobic core (9,532 atoms) |
| Windows | 37, every 0.1 nm from 3.6 nm to 0 (one leaflet; the bilayer is symmetric), 1000 kJ/mol/nm² |

The ligand is the neutral, open-chain lapachol (C15H14O3, 32 atoms, net charge 0.000).

## Results so far

Pure bilayer, 100 ns:

| Quantity | Result |
|---|---|
| Temperature | 303.20 K (target 303.15) |
| Pressure | 0.89 bar (fluctuations of about ±230 bar, expected for this size) |
| Density | 1016.4 kg/m³ |
| Area per lipid | 64.4 ± 0.3 Å² (mean of ten 10 ns blocks, standard error); oscillates about ±2 Å² with 10–20 ns periods |
| P–P thickness | 3.86 nm in the frame used to start the pulling (68.7 ns) |

The linear change of the area in the last 50 ns was −2.2 Å², larger than the 0.5 Å² criterion I had set. I judged the criterion too strict for 60 lipids per leaflet, since slow oscillations dominate the slope; the second half of the run is about 1 Å² above the first, which 100 ns cannot separate from fluctuation.

Lapachol in water and steered MD:

| Quantity | Result |
|---|---|
| Position during 1 ns equilibration | 3.598 ± 0.047 nm (expected for the spring: about 0.05 nm) |
| Steered MD (1.8 ns, 0.002 nm/ps) | z from 3.614 to 0.110 nm; spring force from −208 to 127 kJ/mol/nm |

Speed on a Colab L4: 310 ns/day (pure bilayer) and 344 ns/day (with lapachol). On 4 CPU cores: 4.8 ns/day.

The free-energy profile does not exist yet; it will be added after stage 1 and the overlap check.

## Repository layout

```
inputs/bilayer/    starting structure, topology, index and production .mdp of the pure bilayer
inputs/ligand/     lapachol coordinates named as in LIG.itp
inputs/system/     topology, index and .mdp files of the system with lapachol, windows.dat
inputs/umbrella/   .mdp of each window (stage 1)
results/           figures, logs and small outputs
scripts/           workflow helpers and analysis
docs/              open notebook (GitHub Pages)
```

Trajectories, checkpoints and binary run inputs are not tracked because of size. They are kept on Google Drive and can be shared on request.

## How to reproduce

All gmx steps used GROMACS 2024.5. `scripts/umbrella_workflow.py` generates the inputs; the gmx commands go between its steps.

```bash
# 1. pure bilayer, 100 ns
gmx grompp -f prod_100ns.mdp -c step6.6_equilibration.gro -p topol_bilayer.top -n index.ndx -o prod_100ns.tpr
gmx mdrun -deffnm prod_100ns -cpi prod_100ns.cpt -cpt 10 -maxh 11 -nb gpu -pme gpu -bonded gpu -update gpu
printf "Temperature\nPressure\nBox-X\nBox-Y\nBox-Z\nDensity\n0\n" | gmx energy -f prod_100ns.edr -o prod_energy.xvg
python scripts/analysis_bilayer.py prod_energy.xvg

# 2. pick a frame near the mean area and box height, make molecules whole and center the bilayer
gmx trjconv -f prod_100ns.xtc -s prod_100ns.tpr -dump 68700 -o frame_pulling.gro
printf "POPC\nSystem\n" | gmx trjconv -f frame_pulling.gro -s prod_100ns.tpr -pbc mol -center -o frame_whole.gro

# 3. lapachol in water, minimization, equilibration
python scripts/umbrella_workflow.py insert
gmx grompp -f em.mdp -c system_with_lapachol.gro -p topol_system.top -n index_system.ndx -o em.tpr && gmx mdrun -deffnm em
gmx grompp -f eq_lig.mdp -c em.gro -p topol_system.top -n index_system.ndx -o eq_lig.tpr && gmx mdrun -deffnm eq_lig

# 4. steered MD and window frames
python scripts/umbrella_workflow.py pull-mdp
gmx grompp -f pull.mdp -c eq_lig.gro -t eq_lig.cpt -p topol_system.top -n index_system.ndx -o pull.tpr && gmx mdrun -deffnm pull
python scripts/umbrella_workflow.py frames

# 5. umbrella windows (stage 1) and WHAM
python scripts/umbrella_workflow.py umbrella --ns 5
python scripts/umbrella_workflow.py run
bash scripts/wham.sh        # not run yet
```

## References

1. Torrie GM, Valleau JP. *J Comput Phys* 1977, 23, 187.
2. Kumar S, et al. *J Comput Chem* 1992, 13, 1011.
3. Hub JS, de Groot BL, van der Spoel D. *J Chem Theory Comput* 2010, 6, 3713.
4. Marrink SJ, Berendsen HJC. *J Phys Chem* 1994, 98, 4155.
5. Jo S, et al. *J Comput Chem* 2008, 29, 1859 (CHARMM-GUI).
6. Cordeiro RM. *J Phys Chem B* 2018, 122, 8211.

## Author

Marcelo Souza. Chemistry and biology teacher building a computational chemistry portfolio.
