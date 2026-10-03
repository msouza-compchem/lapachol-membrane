#!/usr/bin/env python3
"""Helpers for the lapachol-membrane workflow (GROMACS 2024.5 + MDAnalysis).

Run each step from the folder that holds the GROMACS files (on Colab:
/content/drive/MyDrive/lapachol-membrane/bilayer_prod). The gmx commands that
go between the steps are listed in the README.

  insert    put lapachol in the water, build topol_system.top, index and .mdp files
  pull-mdp  make pull.mdp (steered MD from eq_lig.mdp)
  frames    pick one frame per window from the pulling trajectory
  umbrella  make the .mdp and .tpr of every window
  run       run the windows, skipping the finished ones
"""
import argparse, os, re, subprocess, time
import numpy as np

EM_MDP = """integrator = steep
emtol = 1000.0
emstep = 0.001
nsteps = 5000
cutoff-scheme = Verlet
nstlist = 10
pbc = xyz
coulombtype = PME
rcoulomb = 1.2
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw_switch = 1.0
rvdw = 1.2
DispCorr = no
constraints = h-bonds
constraint_algorithm = LINCS
"""

EQ_MDP = """integrator = md
dt = 0.002
nsteps = 500000
nstxout-compressed = 5000
nstcalcenergy = 100
nstenergy = 1000
nstlog = 1000
cutoff-scheme = Verlet
nstlist = 20
rlist = 1.2
vdwtype = Cut-off
vdw-modifier = Force-switch
rvdw_switch = 1.0
rvdw = 1.2
coulombtype = PME
rcoulomb = 1.2
DispCorr = no
tcoupl = v-rescale
tc_grps = MEMB SOLV
tau_t = 1.0 1.0
ref_t = 303.15 303.15
pcoupl = C-rescale
pcoupltype = semiisotropic
tau_p = 5.0
compressibility = 4.5e-5 4.5e-5
ref_p = 1.0 1.0
constraints = h-bonds
constraint_algorithm = LINCS
continuation = no
gen_vel = yes
gen_temp = 303.15
gen_seed = -1
nstcomm = 100
comm_mode = linear
comm_grps = MEMB SOLV
pull = yes
pull-ncoords = 1
pull-ngroups = 2
pull-group1-name = MEMCORE
pull-group1-pbcatom = {pbcatom}
pull-group2-name = LIG
pull-coord1-type = umbrella
pull-coord1-geometry = direction
pull-coord1-vec = 0 0 1
pull-coord1-dim = N N Y
pull-coord1-groups = 1 2
pull-coord1-start = yes
pull-coord1-rate = 0
pull-coord1-k = 1000
pull-nstxout = 500
pull-nstfout = 500
"""


def setkey(text, key, value):
    """Set `key = value` in an .mdp text, adding the line if it is missing."""
    pat = rf"^{re.escape(key)}\s*=.*$"
    if re.search(pat, text, flags=re.M):
        return re.sub(pat, f"{key} = {value}", text, flags=re.M)
    return text + f"{key} = {value}\n"


def write_group(f, name, atoms):
    idx = atoms.indices + 1
    f.write(f"[ {name} ]\n")
    for i in range(0, len(idx), 15):
        f.write(" ".join(map(str, idx[i:i + 15])) + "\n")


def cmd_insert(a):
    import MDAnalysis as mda
    u = mda.Universe(a.bilayer)
    lig = mda.Universe(a.ligand)
    box = u.dimensions.copy()
    popc = u.select_atoms("resname POPC")
    zc = popc.center_of_mass()[2]
    target = np.array([box[0] / 2, box[1] / 2, zc + a.z * 10.0])      # Å
    lig.atoms.translate(target - lig.atoms.center_of_mass())
    lig.residues.resids = [int(u.residues.resids.max()) + 1]
    mm = mda.Merge(u.atoms, lig.atoms)
    mm.dimensions = box
    clash = mm.select_atoms("same residue as (resname TIP3 and around 3.0 resname LIG)")
    keep = mm.atoms.difference(clash)
    keep.write("system_with_lapachol.gro")
    nw = keep.select_atoms("resname TIP3").n_residues
    print(f"waters removed: {clash.n_residues} | TIP3 left: {nw} | atoms: {keep.n_atoms}")

    top = open(a.topology).read()
    top = re.sub(r"^TIP3\s+\d+", f"TIP3           {nw}", top, flags=re.M)
    if not top.endswith("\n"):
        top += "\n"
    open("topol_system.top", "w").write(top + "LIG              1\n")

    s = mda.Universe("system_with_lapachol.gro")
    spopc = s.select_atoms("resname POPC")
    zc = spopc.center_of_mass()[2]
    core = spopc.select_atoms(f"prop z > {zc - 12} and prop z < {zc + 12}")   # hydrophobic core
    pbcatom = int(core.atoms[np.argmin(np.abs(core.positions[:, 2] - zc))].index) + 1
    with open("index_system.ndx", "w") as f:
        write_group(f, "System", s.atoms)
        write_group(f, "MEMB", spopc)
        write_group(f, "SOLV", s.select_atoms("not resname POPC"))
        write_group(f, "LIG", s.select_atoms("resname LIG"))
        write_group(f, "MEMCORE", core)
    open("em.mdp", "w").write(EM_MDP)
    open("eq_lig.mdp", "w").write(EQ_MDP.format(pbcatom=pbcatom))
    print(f"MEMCORE: {len(core)} atoms | reference atom {pbcatom}")


def cmd_pull_mdp(a):
    t = open("eq_lig.mdp").read()
    for k, v in [("nsteps", "900000"), ("nstxout-compressed", "500"), ("continuation", "yes"),
                 ("gen_vel", "no"), ("pull-coord1-rate", f"-{a.rate}"),
                 ("pull-pbc-ref-prev-step-com", "yes"), ("pull-group1-pbcatom", "-1")]:
        t = setkey(t, k, v)
    open("pull.mdp", "w").write(t)
    print("pull.mdp written")


def cmd_frames(a):
    x = np.loadtxt("pull_pullx.xvg", comments=("#", "@"))
    assert x[:, 1].min() < 0.3, "the pulling did not reach the bilayer center"
    targets = np.arange(a.zmax, -0.0001, -a.step)
    rows = []
    for i, z0 in enumerate(targets):
        k = np.argmin(np.abs(x[:, 1] - z0))
        tp = x[k, 0]
        r = subprocess.run(["gmx", "trjconv", "-f", "pull.xtc", "-s", "pull.tpr", "-dump", str(tp),
                            "-o", f"conf{i:02d}.gro"], input="0\n", text=True, capture_output=True)
        if r.returncode != 0:
            print("trjconv failed for window", i, r.stderr[-300:])
        rows.append((i, z0, tp, x[k, 1]))
    np.savetxt("windows.dat", np.array(rows), fmt="%d %.2f %.1f %.3f",
               header="window z_target(nm) t(ps) z_real(nm)")
    print(len(rows), "windows | largest |z_target - z_real|:", f"{max(abs(r[1] - r[3]) for r in rows):.3f} nm")


def cmd_umbrella(a):
    base = open("pull.mdp").read()
    W = np.loadtxt("windows.dat")
    os.makedirs("umbrella", exist_ok=True)
    steps = int(a.ns * 1000 / 0.002)
    for i, z in enumerate(W[:, 1]):
        t = base
        for k, v in [("nsteps", steps), ("nstxout-compressed", 5000), ("continuation", "no"),
                     ("gen_vel", "yes"), ("gen_temp", "303.15"), ("gen_seed", "-1"),
                     ("pull-coord1-start", "no"), ("pull-coord1-init", f"{z:.3f}"),
                     ("pull-coord1-rate", "0"), ("pull-nstxout", 500), ("pull-nstfout", 500)]:
            t = setkey(t, k, v)
        open(f"umbrella/umb_{i:02d}.mdp", "w").write(t)
        r = subprocess.run(["gmx", "grompp", "-f", f"umbrella/umb_{i:02d}.mdp", "-c", f"conf{i:02d}.gro",
                            "-p", "topol_system.top", "-n", "index_system.ndx",
                            "-o", f"umbrella/umb_{i:02d}.tpr", "-po", f"umbrella/mdout_{i:02d}.mdp"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("grompp failed for window", i, r.stderr[-400:])
    print(len(W), "windows prepared,", a.ns, "ns each")


def cmd_run(a):
    n = len(np.loadtxt("windows.dat"))
    t0 = time.time()
    for i in range(n):
        d = f"umbrella/umb_{i:02d}"
        if os.path.exists(d + ".gro"):
            continue
        if (time.time() - t0) / 3600 > a.hours - 0.4:
            print("session budget reached; run again to continue")
            break
        r = subprocess.run(["gmx", "mdrun", "-deffnm", d, "-cpi", d + ".cpt", "-cpt", "10",
                            "-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu"],
                           capture_output=True, text=True)
        if r.returncode != 0:
            print("mdrun failed for window", i, r.stderr[-500:])
            break
        px = np.loadtxt(d + "_pullx.xvg", comments=("#", "@"))
        s = px[px[:, 0] >= 500]
        print(f"window {i:02d} | z mean {s[:, 1].mean():.3f} +- {s[:, 1].std():.3f} nm")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("insert")
    s.add_argument("--bilayer", default="frame_whole.gro")
    s.add_argument("--ligand", default="lig_named.pdb")
    s.add_argument("--topology", default="topol_bilayer.top")
    s.add_argument("--z", type=float, default=3.6, help="start distance from the bilayer center (nm)")
    s.set_defaults(f=cmd_insert)
    s = sub.add_parser("pull-mdp")
    s.add_argument("--rate", type=float, default=0.002, help="pulling rate (nm/ps)")
    s.set_defaults(f=cmd_pull_mdp)
    s = sub.add_parser("frames")
    s.add_argument("--zmax", type=float, default=3.6)
    s.add_argument("--step", type=float, default=0.1)
    s.set_defaults(f=cmd_frames)
    s = sub.add_parser("umbrella")
    s.add_argument("--ns", type=float, default=5.0)
    s.set_defaults(f=cmd_umbrella)
    s = sub.add_parser("run")
    s.add_argument("--hours", type=float, default=10.0)
    s.set_defaults(f=cmd_run)
    a = p.parse_args()
    a.f(a)
