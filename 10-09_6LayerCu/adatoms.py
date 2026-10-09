"""Count adatoms and surface roughness over time for every 50 ps trajectory.

Layers are defined from the first frame: 4 Cu layers of 48 atoms, the top one is layer 0.
Each Cu is assigned to the nearest layer by its height above the initial top layer:
    layer 0  = top layer, layer +1 = adatom, layer >= +2 = further out (dissolving)
The selected (biased) Cu, ASE index 118, is left out of all counts and reported separately.

Usage:  python adatoms.py [path/to/trajectories_50ps]
Writes adatoms.csv and adatoms.png to the current directory.
"""
import sys, glob, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import read

TRAJ_DIR = sys.argv[1] if len(sys.argv) > 1 else "cu_data_trajectories_2026-10-07/trajectories_50ps"
SELECTED = 118
N_PER_LAYER = 48

rows = []
for path in sorted(glob.glob(os.path.join(TRAJ_DIR, "*_every50ps.extxyz.gz"))):
    run = os.path.basename(path).split("_every")[0]          # e.g. wCO_w0
    frames = read(path, index=":")                           # ASE reading all frames
    print(run, len(frames), "frames")

    # reference layers from the first frame
    cu0 = np.where(frames[0].numbers == 29)[0]
    z0 = np.sort(frames[0].positions[cu0, 2])
    top = z0[-N_PER_LAYER:].mean()
    second = z0[-2 * N_PER_LAYER:-N_PER_LAYER].mean()
    d = top - second                                          # layer spacing

    for atoms in frames:
        cu = np.where(atoms.numbers == 29)[0]
        cu = cu[cu != SELECTED]
        z = atoms.positions[cu, 2]
        layer = np.round((z - top) / d).astype(int)
        surface = z[layer >= 0]                               # top layer + everything above it
        rows.append(dict(
            run=run,
            system=run.split("_")[0],
            time_ns=atoms.info["campaign_time_ns"],
            n_adatoms=int((layer == 1).sum()),
            n_above_adatom=int((layer >= 2).sum()),
            n_top_layer=int((layer == 0).sum()),              # 47 when intact (48 minus selected Cu)
            roughness_A=surface.std(),                        # std of heights of top layer + adatoms
            selected_cu_height_A=atoms.positions[SELECTED, 2] - top,
        ))

df = pd.DataFrame(rows)
df.to_csv("adatoms.csv", index=False)

fig, ax = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
for j, system in enumerate(["nCO", "wCO"]):
    for run, g in df[df.system == system].groupby("run"):
        ax[0, j].plot(g.time_ns, g.n_adatoms, label=run)
        ax[1, j].plot(g.time_ns, g.roughness_A, label=run)
    ax[0, j].set_title(system)
    ax[0, j].set_ylabel("adatoms (excl. selected Cu)")
    ax[1, j].set_ylabel("roughness (A)")
    ax[1, j].set_xlabel("time (ns)")
    ax[0, j].legend(fontsize=7)
fig.tight_layout()
fig.savefig("adatoms.png", dpi=150)
print(df.groupby("run")[["n_adatoms", "roughness_A"]].agg(["mean", "max"]).round(2))
