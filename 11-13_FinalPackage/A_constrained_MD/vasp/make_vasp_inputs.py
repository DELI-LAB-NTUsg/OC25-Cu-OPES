"""eSEN window frames -> VASP constrained-MD folders (one per window), run in A_constrained_MD/vasp/:

    python3 make_vasp_inputs.py <name> [scan_dir]      e.g. 4x4ad_nCO_00   (scan_dir default ../esen)

Reads <scan_dir>/win_<name>_k00.xyz ... k12.xyz (last frame of each 100 ps eSEN window: positions + velocities),
<scan_dir>/inputs/EQ_<name>.xyz (fixed bottom layers) and <scan_dir>/inputs/target_<name>.txt, and the window heights
from <scan_dir>/scan_<name>.csv. Writes <name>/k<kk>_<up|down>_h<h>/ with
  POSCAR    species order Cu, C, O, H, Cs (those present); selective dynamics: target T T F, bottom 2 Cu layers F F F,
            all else T T T; velocities in A/fs
  INCAR     OC25 level + MD: MDALGO 3 (Langevin 5/ps per species), 300 K, POTIM 0.5, NSW 40000 (20 ps)
  KPOINTS   4x4x1 (OC25 rule round(40/a) for the 10.1 A cell)
  ICONST    wCO only: "R C target 4" = Fermi-type wall on the C-target distance (FBIAS_R0 2.3 A, A 3 eV, D 100)
  index_map.txt  POSCAR line of the target (and C), window height, eSEN mean F_z of that window
"""
import os
import sys

import numpy as np
import pandas as pd
from ase import units
from ase.io import read

name = sys.argv[1]
src = sys.argv[2] if len(sys.argv) > 2 else "../esen"
SPECIES = ["Cu", "C", "O", "H", "Cs"]
NSW = 40000

INCAR = """SYSTEM = {name} {win}, constrained-height MD, OC25 level
GGA    = RP
IVDW   = 11
ENCUT  = 400
PREC   = Normal
ALGO   = Normal
ISMEAR = 0
SIGMA  = 0.1
ISPIN  = 1
LASPH  = .TRUE.
LREAL  = Auto
LDIPOL = .TRUE.
IDIPOL = 3
DIPOL  = 0.5 0.5 0.5
EDIFF  = 1E-5
NELM   = 100
ISYM   = 0
IBRION = 0
MDALGO = 3
ISIF   = 2
TEBEG  = 300
POTIM  = 0.5
NSW    = {nsw}
LANGEVIN_GAMMA = {gamma}
LCHARG = .FALSE.
LWAVE  = .FALSE.
"""
WALL = "FBIAS_R0 = 2.3\nFBIAS_A  = 3.0\nFBIAS_D  = 100\n"

eq = read(f"{src}/inputs/EQ_{name}.xyz")
fixed = {i for c in eq.constraints for i in c.index}
d = {k: [int(x) for x in v] for k, *v in (l.split() for l in open(f"{src}/inputs/target_{name}.txt"))}
t, c = d["target_ase"][0], d.get("C_ase", [None])[0]
top = [i for i in d["top_layer_ase"] if i != t]
scan = pd.read_csv(f"{src}/scan_{name}.csv")
scan = scan[scan.k >= 0]
assert len(fixed) == 32, "expected 32 fixed Cu (bottom 2 layers)"

n_win = scan.k.max() + 1
for k in range(n_win):
    w = scan[scan.k.eq(k)]
    h_set, direc = w.h_set.iloc[0], w.dir.iloc[0]
    a = read(f"{src}/win_{name}_k{k:02d}.xyz")
    assert len(a) == len(eq) and (a.numbers == eq.numbers).all()
    sym = np.array(a.get_chemical_symbols())
    order = [i for s in SPECIES for i in np.where(sym == s)[0]]
    assert len(order) == len(a), f"unexpected element in {name}"
    present = [s for s in SPECIES if (sym == s).any()]
    v = a.get_velocities() * units.fs                               # A/fs
    assert np.abs(v).max() > 0, "no velocities in the window frame"
    h = a.positions[t, 2] - a.positions[top, 2].mean()
    win = f"k{k:02d}_{direc}_h{h_set:.2f}"
    out = os.path.join(name, win)
    os.makedirs(out, exist_ok=False)
    lines = [f"{name} {win}", "1.0"] + [" ".join(f"{x:.10f}" for x in row) for row in a.cell.array]
    lines += [" ".join(present), " ".join(str((sym == s).sum()) for s in present), "Selective dynamics", "Cartesian"]
    for i in order:
        flag = "F F F" if i in fixed else ("T T F" if i == t else "T T T")
        lines.append(" ".join(f"{x:.10f}" for x in a.positions[i]) + " " + flag)
    lines += ["", "Cartesian"] + [" ".join(f"{x:.10f}" for x in v[i]) for i in order]
    open(os.path.join(out, "POSCAR"), "w").write("\n".join(lines) + "\n")
    incar = INCAR.format(name=name, win=win, nsw=NSW, gamma=" ".join(["5.0"] * len(present)))
    tp = order.index(t) + 1
    if c is not None:
        cp = order.index(c) + 1
        open(os.path.join(out, "ICONST"), "w").write(f"R {cp} {tp} 4\n")
        incar += WALL
    open(os.path.join(out, "INCAR"), "w").write(incar)
    open(os.path.join(out, "KPOINTS"), "w").write("OC25 rule round(40/a), round(40/b), 1\n0\nMonkhorst-Pack\n4 4 1\n0 0 0\n")
    with open(os.path.join(out, "index_map.txt"), "w") as f:
        f.write(f"species {' '.join(present)}\ntarget_poscar {tp}\n")
        if c is not None:
            f.write(f"C_poscar {cp}\n")
        f.write(f"h_set {h_set}\nh_start {h:.4f}\neSEN_mean_F_z {w[w.time_ps > 1.0].F_z.mean():.4f}\n")
    print(f"{out}: {len(a)} atoms, target POSCAR line {tp}, h {h:.3f} A")
