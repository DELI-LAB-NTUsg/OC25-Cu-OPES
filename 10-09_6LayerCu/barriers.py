#!/usr/bin/env python
# coding: utf-8

# # Barrier Calculation and Plotting - OC25
# Free-energy barriers along height for lattice, adatom and solvated Cu, using PLUMED's FES_from_State.py.
# 
# For each walker, FES_from_State.py turns the final OPES State.data into a 2D FES F(height, cno).
# That is projected onto height:  F(h) = -kT log sum_cno exp(-F(h, cno) / kT).
# Barrier A->B = highest F(h) between the two basins minus the lowest F(h) in basin A.
# 
# Needs FES_from_State.py (github.com/Miro-Astore/plumed_scripts or plumed2/plugins/opes) next to this file.
# Usage:  python barriers.py [path/to/walker_inputs_and_states]
# Writes fes_<run>.dat, barriers.csv, barriers.png and barriers_summary.png to the current directory.

# In[14]:


# SETUP
import sys, glob, os, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

STATE_DIR = "cu_data_trajectories_2026-10-07/walker_inputs_and_states"
FES_SCRIPT = "FES_from_State.py"
KT = 8.617333e-5 * 300                     # eV, 300 K (PLUMED energies are in eV here)
BASINS = {"lattice": (-0.75, 0.75), "adatom": (1.25, 2.25), "solvated": (5.0, 7.5)}   # height ranges (A)

states = sorted(glob.glob(os.path.join(STATE_DIR, "*", "State.data")))
if not states:
    sys.exit(f"No */State.data found in {STATE_DIR} (pass the walker_inputs_and_states folder)")

rows, profiles, grids = [], {}, {}


# In[15]:


"""
for state in states:
    run = os.path.basename(os.path.dirname(state))              # e.g. wCO_w0
    out = f"fes_{run}.dat"
    subprocess.run([sys.executable, FES_SCRIPT, "-f", state, "-o", out, "--kt", str(KT),
                    "--min=-3,0", "--max=9,8", "--bin=120,80"], check=True)  # "=" so "-3,0" isn't read as a flag
"""


# Log-Sum-Exp trick, going from $F(h,cno)\rightarrow F(h)$ $F(h)=-k_BT\cdot\ln\left(\sum\limits_\text{cno}e^{-(F(h,cno)/k_BT)}\right)$
# 
# For each height, all `cno` values are combined into one free energy. 
# Fmin is the minimum value of all results for one walker, ensuring the exonentials stay between 0 and 1. Then the Fh.min() i s subtracted normalises the plot to have the lowest point at 0 eV.
# 
# Barriers obtained from finding the lowest F inside each basin's height range: bottom of basin, finding highest F between basins as the top that needs to be crossed.

# In[16]:


for state in states:
    run = os.path.basename(os.path.dirname(state))              # e.g. wCO_w0
    out = f"fes_{run}.dat"

    # Read free energy from CSV
    fes = pd.read_csv(out, sep=r"\s+", comment="#", names=["height", "cno", "F"])

    # MAP free energy onto a grid of height vs cno
    grid = fes.pivot(index="height", columns="cno", values="F")  # rows = height, columns = cno

    # Compute free energy profile along height by integrating out cno
    h = grid.index.values
    Fmin = grid.values.min(axis=1)

    # Use the log-sum-exp trick to compute the free energy profile along height
    Fh = Fmin - KT * np.log(np.exp(-(grid.values - Fmin[:, None]) / KT).sum(axis=1))
    Fh -= Fh.min()
    profiles[run] = (h, Fh)
    grids[run] = grid

    row = dict(run=run, system=run.split("_")[0])
    for a, b in [("lattice", "adatom"), ("adatom", "solvated"), ("lattice", "solvated")]:
        in_a = (h > BASINS[a][0]) & (h < BASINS[a][1])
        in_b = (h > BASINS[b][0]) & (h < BASINS[b][1])
        top = Fh[(h >= BASINS[a][0]) & (h <= BASINS[b][1])].max()
        row[f"{a}->{b}"] = top - Fh[in_a].min()
        row[f"{b}->{a}"] = top - Fh[in_b].min()
    rows.append(row)

res = pd.DataFrame(rows)
res.to_csv("barriers.csv", index=False)
print(res.drop(columns="run").groupby("system").agg(["mean", "sem"]).round(2).T)


# In[17]:


# F(height) profiles per walker, basins shaded
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for j, system in enumerate(["nCO", "wCO"]):
    for run, (h, Fh) in profiles.items():
        if run.startswith(system):
            ax[j].plot(h, Fh, label=run)
    for lo, hi in BASINS.values():
        ax[j].axvspan(lo, hi, color="grey", alpha=0.15)
    ax[j].set_title(system)
    ax[j].set_xlabel("height (A)")
    ax[j].legend(fontsize=7)
ax[0].set_ylabel("F (eV)")
fig.tight_layout()
fig.savefig("barriers.png", dpi=150)


# In[18]:


# barrier comparison: bar = mean over walkers, error bar = sem, dots = individual walkers
cols = [c for c in res.columns if "->" in c]
x = np.arange(len(cols))
fig, ax = plt.subplots(figsize=(10, 4))
for k, system in enumerate(["nCO", "wCO"]):
    g = res[res.system == system][cols]
    ax.bar(x + (k - 0.5) * 0.4, g.mean(), 0.4, yerr=g.sem(), capsize=3, label=system, alpha=0.7)
    for i, c in enumerate(cols):
        ax.plot(np.full(len(g), x[i] + (k - 0.5) * 0.4), g[c], "k.", ms=4)
ax.set_xticks(x)
ax.set_xticklabels(cols, rotation=20)
ax.set_ylabel("barrier (eV)")
ax.legend()
fig.tight_layout()
fig.savefig("barriers_summary.png", dpi=150)

# Print averages
print("Barrier summary (mean ± sem, eV):")
for system in ["nCO", "wCO"]:
    g = res[res.system == system][cols]
    print(f"  {system}")
    for c in cols:
        print(f"    {c:20s} {g[c].mean():.2f} ± {g[c].sem():.2f}  (n={g[c].count()})")


# In[19]:


# 2D FES per walker from FES_from_State.py, basin height ranges marked
fig, axes = plt.subplots(2, 6, figsize=(20, 6), sharex=True, sharey=True)
for ax, (run, grid) in zip(axes.flat, sorted(grids.items())):
    h, c = grid.index.values, grid.columns.values
    im = ax.imshow(grid.values.T, origin="lower", aspect="auto", vmin=0, vmax=3, cmap="Blues",
                   extent=(h[0], h[-1], c[0], c[-1]))
    for lo, hi in BASINS.values():
        ax.axvline(lo, color="w", lw=0.8)
        ax.axvline(hi, color="w", lw=0.8)
    ax.set_title(run)
for ax in axes[1]:
    ax.set_xlabel("height (A)")
for ax in axes[:, 0]:
    ax.set_ylabel("cno")
fig.colorbar(im, ax=axes, label="F (eV)")
fig.savefig("fes_2d.png", dpi=120)


# In[ ]:




