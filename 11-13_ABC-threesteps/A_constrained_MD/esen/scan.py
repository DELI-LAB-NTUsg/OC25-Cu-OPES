"""Constrained-z scan with eSEN-OC25, done the way it would be done in VASP.

    python -u scan.py <name> <ps> <h0> <mode> [step] [top] [k0 r]   mode = up | updown | eq (inputs/target_<name>.txt)
    adatom equilibration: python -u scan.py 4x4adcs_nCO_00 100 1.75 eq   BUILT_<name>.xyz -> FIRE (target z fixed)
                          + 100 ps at h0 -> inputs/EQ_<name>.xyz (log eq_<name>.csv), then stops
    this package:         python -u scan.py 4x4ad_nCO_00 100 1.75 updown 0.75 6.25
                          (1.75, 2.50, ..., 6.25 A up = windows k00-k06, then 5.50 ... 1.75 A down = k07-k12)
Default step 0.25 A, top 6.0 A. Scans read inputs/EQ_<name>.xyz and run 10 ps at the first height before the first
window (k = -1, dir = eq).
Restart after a water break (roll back):  ... updown 0.75 6.25 <k0> <r>   reruns from window k0 on, starting from
win_<name>_k<k0-1>.xyz (k0 = 0: end of the run-in; k0 = -1: from inputs/EQ_<name>.xyz with the run-in), with new
velocities and thermostat seeds (restart number r = 1, 2, 3); log scan_<name>_r<r>.csv, window files k0... rewritten.

wCO: the CO is kept on the target by an upper wall on the C-target distance (2.3 A, kappa 2000 eV/A^2, PLUMED), in
the equilibration and every window; F_z includes the wall's force (zero while the CO is bound). VASP: ICONST R on
C-target with a Fermi-type upper wall (status 4).
VASP equivalent: selective dynamics with the target `T T F` (z frozen, x/y free), bottom layers `F F F`; MDALGO = 3
(Langevin, LANGEVIN_GAMMA = 5 ps^-1), POTIM = 0.5 fs, normal H mass, 300 K. One window = one MD job; the next
window starts from the end of the previous one with the target moved by one step; h is measured from the
mean z of the other top-layer Cu (in the starting structure).
Moving between windows: the target (wCO: together with its CO, C and O) goes to the next height in 0.05 A increments
with 20 fs of MD after each (0.3 ps per 0.75 A step, not recorded), so the Cu-C bond is never compressed or stretched
and the water can follow. Each window's time_ps starts once the target is at its height.

Output scan_<name>.csv every 2 fs: k, h_set, time_ps, h, F_z (force on the frozen z, eV/A; VASP prints it in OUTCAR),
cno, cn, C-Cu, E_pot. win_<name>_k<kk>.xyz = last frame of each window. dF/dh = -<F_z>  (ti.py).
"""
import sys

import numpy as np
from ase import units
from ase.constraints import FixAtoms, FixCartesian
from ase.io import read
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import FIRE

from cv import CVs, calculator, read_target
from water_check import WaterWatchdog

RAMP_DZ, RAMP_STEPS = 0.05, 40      # between windows: move 0.05 A at a time, 40 steps (20 fs) of MD after each

name = sys.argv[1]
ps, h0, mode = float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
step = float(sys.argv[5]) if len(sys.argv) > 5 else 0.25         # window spacing (A)
top = float(sys.argv[6]) if len(sys.argv) > 6 else 6.0           # highest window (A)
k0 = int(sys.argv[7]) if len(sys.argv) > 7 else -1               # restart: first window to run again
rs = int(sys.argv[8]) if len(sys.argv) > 8 else 0                # restart number (new seeds, new csv)
up = np.round(np.arange(h0, top + 1e-4, step), 2)
grid = np.concatenate([up, up[-2::-1]]) if mode == "updown" else up

atoms = read(f"inputs/{'BUILT' if mode == 'eq' else 'EQ'}_{name}.xyz")
fixed = [c for c in atoms.constraints if isinstance(c, FixAtoms)]
atoms.pbc = True
target, carbon, top = read_target(f"inputs/target_{name}.txt")
atoms.set_constraint(atoms.constraints + [FixCartesian(target, mask=(False, False, True))])
atoms.calc = calculator()
if carbon is not None:                       # wCO: keep the CO on the target (as the OPES runs), whole run
    from ase.calculators.plumed import Plumed
    atoms.calc = Plumed(calc=atoms.calc, atoms=atoms, timestep=0.5 * units.fs, kT=units.kB * 300,
                        log=f"plumed_{name}.log",
                        input=["UNITS LENGTH=A TIME=fs ENERGY=eV",
                               f"c_cu: DISTANCE ATOMS={carbon + 1},{target + 1}",
                               "uwall: UPPER_WALLS ARG=c_cu AT=2.3 KAPPA=2000.0"])
cvs = CVs(atoms, target, carbon, top)
co = []                                      # wCO: C and O of the CO, moved together with the target
if carbon is not None:
    o_all = np.where(atoms.numbers == 8)[0]
    co = [carbon, int(o_all[np.argmin(atoms.get_distances(carbon, o_all, mic=True))])]
z_top0 = atoms.positions[top, 2].mean()                # from the starting structure, also on a restart
if k0 >= 0:                                             # restart: continue from the end of window k0 - 1
    atoms.set_positions(read(f"win_{name}_k{k0 - 1:02d}.xyz").positions, apply_constraint=False)
else:
    atoms.positions[target, 2] = z_top0 + grid[0]
if mode == "eq":
    FIRE(atoms, logfile=f"relax_{name}.log").run(fmax=0.1, steps=500)
MaxwellBoltzmannDistribution(atoms, temperature_K=300)

log = open(f"{'eq' if mode == 'eq' else 'scan'}_{name}{f'_r{rs}' if rs else ''}.csv", "w")
log.write("k,dir,h_set,time_ps,h,F_z,cno,cn,c_cu,epot\n")
eq_ps = ps if mode == "eq" else 10
for k, h_set in [(-1, grid[0])] + list(enumerate(grid)):
    if k < k0:                                            # restart: windows before k0 are kept from the first run
        continue
    d = "eq" if k < 0 else "up" if k < len(up) else "down"
    dz = z_top0 + h_set - atoms.positions[target, 2]
    n_ramp = int(round(abs(dz) / RAMP_DZ))
    if n_ramp:                                            # go to the new height gradually (target + CO)
        ramp = Langevin(atoms, timestep=0.5 * units.fs, temperature_K=300, friction=0.005 / units.fs,
                        rng=np.random.default_rng(1000 * (k + 1) + len(name) + 500 + 100000 * rs))
        for _ in range(n_ramp):
            atoms.positions[[target] + co, 2] += dz / n_ramp
            ramp.run(RAMP_STEPS)
    atoms.positions[target, 2] = z_top0 + h_set          # exact window height
    dyn = Langevin(atoms, timestep=0.5 * units.fs, temperature_K=300, friction=0.005 / units.fs,
                   rng=np.random.default_rng(1000 * (k + 1) + len(name) + 100000 * rs))

    def record():
        if dyn.nsteps % 4 == 0:                       # every 2 fs
            h, cno, cn, dc = cvs(atoms)
            f_z = atoms.calc.results["forces"][target, 2]   # raw force; the constraint zeroes only the copy
            log.write(f"{k},{d},{h_set},{dyn.nsteps / 2000:.3f},{h:.4f},{f_z:.5f},{cno:.3f},{cn:.3f},{dc:.3f},"
                      f"{atoms.calc.results['energy']:.5f}\n")

    dyn.attach(record, interval=1)
    dyn.attach(WaterWatchdog(atoms, dyn), interval=2000)
    dyn.run(int((eq_ps if k < 0 else ps) * 2000))     # 0.5 fs steps
    log.flush()
    if mode == "eq":                                      # equilibrated cell for the scan; keep FixAtoms only
        a = atoms.copy()
        a.set_constraint(fixed)
        a.write(f"inputs/EQ_{name}.xyz")
        print("equilibrated ->", f"inputs/EQ_{name}.xyz", flush=True)
        break
    atoms.write(f"win_{name}_k{k:02d}.xyz")
    print(f"window {k}  h_set {h_set:.2f}  done", flush=True)
log.close()
