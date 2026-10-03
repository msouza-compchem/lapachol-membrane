#!/usr/bin/env python3
"""Bilayer equilibration check: area per lipid in 10 ns blocks and a 4-panel figure.

Input: prod_energy.xvg from
  printf "Temperature\nPressure\nBox-X\nBox-Y\nBox-Z\nDensity\n0\n" | gmx energy -f prod_100ns.edr -o prod_energy.xvg
(in GROMACS 2024.5 these terms were 14, 15, 17, 18, 19 and 21).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fname = sys.argv[1] if len(sys.argv) > 1 else "prod_energy.xvg"
n_lipids_per_leaflet = 60

d = np.loadtxt(fname, comments=("#", "@"))
t = d[:, 0] / 1000                                   # ns
T, P, bx, by, bz, rho = d[:, 1:7].T
apl = bx * by / n_lipids_per_leaflet * 100           # Å² per lipid

blocks = []
print("block (ns)   APL mean ± sd (Å²)   Box-Z (nm)")
for i in range(int(t.max() // 10)):
    s = (t >= i * 10) & (t < (i + 1) * 10)
    blocks.append(apl[s].mean())
    print(f"{i*10:3d}-{(i+1)*10:3d}      {apl[s].mean():6.2f} ± {apl[s].std():4.2f}        {bz[s].mean():6.3f}")
b = np.array(blocks)
print(f"\nmean of blocks {b.mean():.2f} Å², standard error {b.std(ddof=1)/np.sqrt(len(b)):.2f} Å²")
s = t > t.max() / 2
slope = np.polyfit(t[s], apl[s], 1)[0]
print(f"second half: mean {apl[s].mean():.2f} Å², slope {slope:+.4f} Å²/ns, total change {slope*t.max()/2:+.2f} Å²")

k = np.ones(1000) / 1000                             # 2 ns running mean (points every 2 ps)
fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
for a, y, name in zip(ax.flat, [apl, bz, rho, T],
                      ["Area per lipid (Å²)", "Box-Z (nm)", "Density (kg/m³)", "Temperature (K)"]):
    a.plot(t, y, alpha=.25, lw=.6)
    a.plot(t[999:], np.convolve(y, k, "valid"), lw=1.6)
    a.set_xlabel("time (ns)")
    a.set_ylabel(name)
fig.savefig("fig_bilayer_equilibration.png", dpi=200)
