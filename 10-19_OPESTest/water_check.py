"""Water integrity check: count broken water in a structure or trajectory, or stop an MD run when it breaks.

Counts per frame: free H (no water O within 1.25 A), H2, OH-, H3O+, bare O, and free H sitting on Cu.
The CO oxygen (O within 1.4 A of a C) is not counted as water.
Expected: 0 for all of them. A free H / OH- / H3O+ for a frame or two is a proton hopping between waters.

1) On files:   python water_check.py [--every N] file [file ...]   any ASE-readable file, .gz ok
               --every N checks every Nth frame. Exit code 1 if any checked frame is not clean.

2) During MD (ASE), stop the run when water stays broken:
       from water_check import WaterWatchdog
       dyn.attach(WaterWatchdog(atoms, dyn), interval=2000)     # check every 1 ps at 0.5 fs
   Raises RuntimeError on any H2, or when water is broken in 10 consecutive checks (10 ps),
   and prints the MD step at which the water was first seen broken.
"""
import os
import sys
import numpy as np
from ase.io import iread
from scipy.spatial import cKDTree

R_OH, R_HH, R_CUH, R_CO = 1.25, 1.0, 2.0, 1.4


def census(atoms):
    sym = np.array(atoms.get_chemical_symbols())
    box = np.where(atoms.pbc, atoms.cell.lengths(), 1e6)       # periodic only where pbc is set
    pos = atoms.positions % box
    iH, iO, iC, iCu = (np.where(sym == s)[0] for s in ("H", "O", "C", "Cu"))
    if len(iC):                                                 # drop CO oxygen(s)
        d = cKDTree(pos[iC], boxsize=box).query(pos[iO])[0]
        iO = iO[d > R_CO]
    dHO, jO = cKDTree(pos[iO], boxsize=box).query(pos[iH])
    nH = np.bincount(jO[dHO < R_OH], minlength=len(iO))
    free = iH[dHO >= R_OH]
    h2 = [(a, b) for a, b in cKDTree(pos[free], boxsize=box).query_pairs(R_HH)] if len(free) > 1 else []
    on_cu = free[cKDTree(pos[iCu], boxsize=box).query(pos[free])[0] < R_CUH] if len(free) and len(iCu) else []
    return dict(free_H=len(free), H2=len(h2), OH=int((nH == 1).sum()), H3O=int((nH == 3).sum()),
                O_bare=int((nH == 0).sum()), H_on_Cu=len(on_cu),
                free_H_idx=free.tolist(), OH_O_idx=iO[nH == 1].tolist())


def broken(c):
    return any(c[x] for x in ("free_H", "H2", "OH", "H3O", "O_bare"))


def describe(c):
    return (f"free H {c['free_H']} H2 {c['H2']} OH- {c['OH']} H3O+ {c['H3O']} bare O {c['O_bare']} "
            f"H on Cu {c['H_on_Cu']} | free H idx (0-based) {c['free_H_idx']} OH- O idx {c['OH_O_idx']}")


class WaterWatchdog:
    """Attach to an ASE dynamics object; stops the run (RuntimeError) on any H2 or on
    broken water in `patience` consecutive checks. Pass `dyn` to get the MD step at which
    the water was first seen broken (roll back to before that step)."""

    def __init__(self, atoms, dyn=None, patience=10):
        self.atoms, self.dyn, self.patience = atoms, dyn, patience
        self.count, self.first_step = 0, None

    def __call__(self):
        c = census(self.atoms)
        if broken(c):
            self.count += 1
            if self.first_step is None and self.dyn is not None:
                self.first_step = self.dyn.nsteps
        else:
            self.count, self.first_step = 0, None
        if c["H2"] or self.count >= self.patience:
            raise RuntimeError(f"water_check: broken water, run stopped. First seen broken at MD step "
                               f"{self.first_step} of this segment ({self.count} checks). {describe(c)}")


def main(args):
    every = 1
    if args and args[0] == "--every":
        every, args = int(args[1]), args[2:]
    bad = 0
    for path in args:
        n = 0
        for k, atoms in zip(range(0, 10**9, every), iread(path, index=f"::{every}")):
            c = census(atoms)
            n += 1
            if broken(c):
                bad += 1
                t = atoms.info.get("campaign_time_ps", atoms.info.get("time", ""))
                print(f"{os.path.basename(path)} frame {k} {t}: {describe(c)}")
        print(f"{os.path.basename(path)}: {n} frames checked")
    print("clean" if bad == 0 else f"{bad} frames with broken water")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
